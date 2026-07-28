"""idl_random.py — bit-exact Python reimplementation of IDL 8.6 RANDOMU / RANDOMN.

Established by controlled experiment against IDL 8.6.0 (linux x86_64), 2026-07-27.
Oracle dumps: probes/rng_ctrl/idl_rng.txt, idl_rng3.txt

Findings
--------
IDL 8.6's *default* random generator is Mersenne Twister MT19937, seeded and advanced
identically to numpy's legacy ``np.random.RandomState``. Given the same scalar seed the
two produce the *same* 32-bit word stream, so the whole IDL stream is reproducible.

Per-call semantics (all verified bit-exact):

==================  ==========================================================
IDL call            equivalent
==================  ==========================================================
randomu(s,n,/doub)  ``rs.random_sample(n)``            — 2 MT words per value
randomu(s,n)        ``float32(word / 2**32)``          — 1 MT word  per value
randomu(s,n,/long)  ``word >> 1``
randomn(s,n,/doub)  ``rs.standard_normal(n)``          — polar over 53-bit doubles
randomn(s,n)        polar over 1-word doubles, float32 rounding (see ``_randomn_float``)
==================  ==========================================================

Spare-deviate rule (the one place IDL and numpy genuinely disagree)
------------------------------------------------------------------
The polar (Marsaglia) method makes deviates in *pairs*. When a call asks for an odd
count, one deviate is left over, and the two precisions handle it differently — verified
independently by four probes each (blocks T/V/X/ZA for float, U/W/Y/N for double):

* ``randomn(s,n)``        (float)  — the spare **persists** across call boundaries, and
  survives an intervening ``randomu``. Six ``randomn(s,1)`` calls give exactly the same
  six numbers as one ``randomn(s,6)``.
* ``randomn(s,n,/double)`` — the spare is **discarded** at the end of each call. Six
  ``randomn(s,1,/double)`` calls return elements 0,2,4,6,8,10 of a single ``randomn(s,6+)``
  call; ``randomn(s,3,/double)`` twice returns elements 0,1,2 then 4,5,6.

numpy caches in both cases, so the double path needs its cache cleared per call.

Not covered: the ``/RAN1`` keyword selects a different (Numerical-Recipes-derived)
generator with its own 36-long state; ``mcmc_dem.pro`` does not use it. The interaction
of a pending *float* spare with an intervening ``/DOUBLE`` normal call is untested.
"""
from __future__ import annotations

import numpy as np

f32 = np.float32

__all__ = ["IDLRandom"]


class IDLRandom:
    """Reproduces IDL 8.6 ``RANDOMU``/``RANDOMN`` bit-for-bit.

    Parameters
    ----------
    seed : int
        The scalar value handed to IDL as the ``Seed`` argument. IDL treats INT, LONG
        and DOUBLE seeds of equal value identically (verified), so any integer works.
    """

    def __init__(self, seed: int):
        self._rs = np.random.RandomState(int(seed))
        self._spare_f = None      # pending float32 deviate, persists across calls

    # ------------------------------------------------------------------ internals
    def _words(self, n: int) -> np.ndarray:
        """Draw n raw MT19937 32-bit words (as float64, exactly representable)."""
        return self._rs.randint(0, 2**32, n, dtype=np.uint32).astype(np.float64)

    def _drop_spare(self) -> None:
        """Clear numpy's cached second gaussian — IDL drops it at each call boundary."""
        st = list(self._rs.get_state())
        st[3], st[4] = 0, 0.0
        self._rs.set_state(tuple(st))

    def _randomn_float(self, n: int) -> np.ndarray:
        """Polar method over 1-word uniforms with IDL's float32 rounding placement.

        Rounding placement (brute-forced against the oracle, unique exact solution):
        uniforms stay float64; v = float32(2u-1); rsq accumulated in float32;
        the log/sqrt transform runs in float64 and is rounded to float32 once.
        """
        out: list[float] = []
        while len(out) < n:
            if self._spare_f is not None:
                out.append(self._spare_f)
                self._spare_f = None
                continue
            while True:
                u1, u2 = self._words(2) / 2.0**32
                v1, v2 = f32(2.0 * u1 - 1.0), f32(2.0 * u2 - 1.0)
                rsq = v1 * v1 + v2 * v2          # float32 arithmetic
                if 0.0 < rsq < 1.0:
                    break
            r = np.float64(rsq)
            fac = f32(np.sqrt(-2.0 * np.log(r) / r))
            out.append(f32(v2 * fac))
            self._spare_f = f32(v1 * fac)        # kept for the next call (IDL float rule)
        return np.array(out[:n], dtype=np.float32)

    # -------------------------------------------------------------------- public
    def randomu(self, n: int | None = None, double: bool = False):
        """IDL ``RANDOMU(seed, n)``; ``double=True`` for the ``/DOUBLE`` keyword."""
        k = 1 if n is None else int(n)
        if double:
            out = self._rs.random_sample(k)
        else:
            out = f32(self._words(k) / 2.0**32)
        return out[0] if n is None else out

    def randomn(self, n: int | None = None, double: bool = False):
        """IDL ``RANDOMN(seed, n)``; ``double=True`` for the ``/DOUBLE`` keyword."""
        k = 1 if n is None else int(n)
        if double:
            out = self._rs.standard_normal(k)
            self._drop_spare()
        else:
            out = self._randomn_float(k)
        return out[0] if n is None else out

    def randomu_long(self, n: int | None = None):
        """IDL ``RANDOMU(seed, n, /LONG)`` — uniform integers in [0, 2**31)."""
        k = 1 if n is None else int(n)
        out = (self._words(k).astype(np.uint64) >> np.uint64(1)).astype(np.int64)
        return out[0] if n is None else out
