"""parity_lib.py — run a parity set and report per-checkpoint verdicts.

Shared by the CLI (run_parity.py) and the regression suite (tests/test_parity.py).
"""
from __future__ import annotations

import os

import numpy as np
from scipy.io import readsav

from mcmc_dem import mcmc_dem
from poa_core import lineflx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # work/pintofale_mcmc_parity

SETS = {
    # name: oracle probe dir + the mcmc_dem configuration the IDL driver used
    "set1": dict(dir="probes/idl_set1", nsim=10, nbatch=5, nburn=2,
                 nosrch=True, abrng=False,
                 desc="/nosrch, frozen abundances (npar=12)"),
    "set2": dict(dir="probes/oracle_set2c", nsim=10, nbatch=5, nburn=2,
                 nosrch=False, abrng=False,
                 desc="grid-search start, frozen abundances (npar=12)"),
    "set3": dict(dir="probes/oracle_set3", nsim=15, nbatch=5, nburn=3,
                 nosrch=True, abrng=True,
                 desc="/nosrch, abundances thawed via ABRNG (npar=42)"),
}


def _probe(oracle, n):
    c = readsav(os.path.join(oracle, f"probe_{n}.sav"), verbose=False)["chk"]
    return {f.lower(): np.asarray(c[f][0])
            for f in c.dtype.names if f.lower() != "probe_id"}


def run_set(name, oracle=None):
    """Run one parity set. Returns (results, npass) where results are
    (label, verdict, detail) tuples and verdict is EXACT / DIFFER / SHAPE."""
    cfg = SETS[name]
    oracle = oracle or os.path.join(ROOT, cfg["dir"])
    results = []

    def probe(n):
        return _probe(oracle, n)

    def check(label, idl, py):
        idl = np.atleast_1d(np.asarray(idl, dtype=np.float64))
        py = np.atleast_1d(np.asarray(py, dtype=np.float64))
        if idl.shape != py.shape:
            results.append((label, "SHAPE", f"idl{idl.shape} vs py{py.shape}"))
            return
        both_nan = np.isnan(idl) & np.isnan(py)
        eq = (idl == py) | both_nan
        if eq.all():
            results.append((label, "EXACT", f"n={idl.size}"))
            return
        d = np.abs(idl - py)
        d[both_nan] = 0
        fin = np.isfinite(d)
        rel = (np.max(d[fin] / np.maximum(np.abs(idl[fin]), 1e-300))
               if fin.any() else np.inf)
        results.append((label, "DIFFER",
                        f"n={idl.size} bad={int((~eq).sum())} maxrel={rel:.3e}"))

    # inputs: EMIS/TRUEDEM are exp()-built test data, taken from the oracle so the
    # comparison isolates the converted algorithm from libm exp differences
    p00 = probe("00_input")
    logt, wvl, zz = p00["logt"], p00["wvl"], p00["zz"]
    emis, truedem, ab = p00["emis"].T.copy(), p00["truedem"].copy(), p00["ab"]
    dem0 = probe("00b_input")["dem0"]

    flx = lineflx(emis, logt, wvl, zz, truedem, ab)
    check("00 flx (converted lineflx)", p00["flx"], flx)
    fsig = 0.1 * flx
    check("00 fsig", p00["fsig"], fsig)

    lscal = probe("01_setup")["lscal"]
    abrng = np.stack([ab / 3.0, ab * 3.0], axis=1) if cfg["abrng"] else None

    tr = {}
    res = mcmc_dem(wvl, flx, emis, zz, logt, dem0, ab, fsig, lscal,
                   nsim=cfg["nsim"], nbatch=cfg["nbatch"], nburn=cfg["nburn"],
                   nosrch=cfg["nosrch"], abrng=abrng, seed=42, trace=tr)

    p01, p01b = probe("01_setup"), probe("01b_setup")
    s = tr["setup"]
    for k in ("allpar", "ipar", "opar", "aparmx", "aparmn", "sigpar", "emenv"):
        check(f"01 {k}", p01[k], s[k])
    for k in ("cenv", "xdem", "xab"):
        check(f"01b {k}", p01b[k], s[k])
    check("01b rngd", p01b["rngd"].T, s["rngd"])
    check("01b rnga", p01b["rnga"].T, s["rnga"])

    p02 = probe("02_init")
    for k in ("dem", "abnd", "ff", "oldprb", "bestprb", "allpar", "sigpar", "bestpar"):
        check(f"02 {k}", p02[k], tr["init"][k])

    for i, b in enumerate(tr["batch"]):
        if not os.path.exists(os.path.join(oracle, f"probe_10_batch_{i:03d}.sav")):
            break
        pb = probe(f"10_batch_{i:03d}")
        for k in ("rn", "ru", "allpar", "sigpar"):
            check(f"10_{i:03d} {k}", pb[k], b[k])

    for i, st in enumerate(tr["step"]):
        if not os.path.exists(os.path.join(oracle, f"probe_11_step_{i:03d}.sav")):
            break
        ps = probe(f"11_step_{i:03d}")
        for k in ("tmpar", "tmpdem", "tmpabnd", "tmpff", "prb", "accept_fn",
                  "metro_check", "allpar"):
            check(f"11_{i:03d} {k}", ps[k], st[k])

    p99, p99b = probe("99_final"), probe("99b_final")
    check("99 dem", p99["dem"], res.dem)
    check("99 demerr", p99["demerr"].T, res.demerr)
    check("99 abund", p99["abund"], res.abund)
    check("99 aberr", p99["aberr"].T, res.aberr)
    check("99 simprb", p99["simprb"], res.simprb)
    for k in ("simdem", "simabn", "simflx"):
        check(f"99 {k}", p99[k].T, getattr(res, k))
    for k in ("storpar", "storidx", "stordem", "bestpar"):
        check(f"99b {k}", p99b[k], getattr(res, k))

    npass = sum(1 for r in results if r[1] == "EXACT")
    return results, npass
