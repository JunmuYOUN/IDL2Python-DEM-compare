"""idl_compat.py — IDL language/library primitives needed for bit-exact parity.

Only the semantics actually exercised by the pintofale mcmc_dem path are
implemented. Each function documents the IDL behaviour it reproduces.
"""
from __future__ import annotations

import ctypes
import ctypes.util

import numpy as np

f32 = np.float32

__all__ = ["total", "interpol_n", "interpol_xu", "histogram", "idl_sort",
           "gt_", "lt_", "where", "alog10_f32", "alog_f32", "exp_f32", "pow10"]

# ---------------------------------------------------------------------------
# IDL's ALOG/ALOG10 on a FLOAT argument call the platform's single-precision
# libm (log10f/logf), which is not correctly rounded. numpy's float32 log10 uses
# its own kernel and disagrees by 1 ulp on some inputs, so libm is called
# directly. Verified: libm.log10f reproduces IDL 30/30 on the abundance vector,
# where f32(np.log10(float64)) missed 3.
# ---------------------------------------------------------------------------
_libm = ctypes.CDLL(ctypes.util.find_library("m") or "libm.so.6")
for _n in ("log10f", "logf"):
    _f = getattr(_libm, _n)
    _f.restype = ctypes.c_float
    _f.argtypes = [ctypes.c_float]


def _apply_f32(fn, a):
    arr = np.asarray(a, dtype=np.float32)
    scalar = arr.ndim == 0
    flat = np.atleast_1d(arr).ravel()
    out = np.array([fn(float(v)) for v in flat], dtype=np.float32)
    return out[0] if scalar else out.reshape(arr.shape)


def alog10_f32(a):
    """IDL ALOG10 on FLOAT input — platform log10f."""
    return _apply_f32(_libm.log10f, a)


def alog_f32(a):
    """IDL ALOG on FLOAT input — platform logf."""
    return _apply_f32(_libm.logf, a)


def exp_f32(a):
    """IDL EXP on FLOAT input — computed in double, then rounded to float32.

    Note the asymmetry with ALOG/ALOG10 above: those match the platform's
    *single*-precision libm, but EXP does not. Verified against the grid-search
    probability profiles: float32(exp(float64(x))) matched 202/202 values, while
    libm expf missed one in each profile and numpy's float32 exp missed ~28%.
    """
    arr = np.asarray(a, dtype=np.float64)
    return np.float32(np.exp(arr))


_libm.pow.restype = ctypes.c_double
_libm.pow.argtypes = [ctypes.c_double, ctypes.c_double]


def pow10(a):
    """IDL `10.D^x` — platform pow(10, x).

    numpy's ** for float64 disagrees with libm pow by 1 ulp on ~1 value in 4
    (verified 7/30 on the abundance vector), so libm is called directly.
    """
    arr = np.asarray(a, dtype=np.float64)
    scalar = arr.ndim == 0
    flat = np.atleast_1d(arr).ravel()
    out = np.array([_libm.pow(10.0, float(v)) for v in flat], dtype=np.float64)
    return out[0] if scalar else out.reshape(arr.shape)


def total(a, nan=False):
    """IDL TOTAL — *sequential* accumulation, unlike numpy's pairwise sum.

    numpy's np.sum uses pairwise summation, which can differ from IDL's simple
    running total in the last ulp. Arrays here are small (<=~100), so an explicit
    loop costs nothing and removes the discrepancy. The accumulator keeps the
    input dtype, matching IDL (TOTAL returns float for float input).
    """
    a = np.asarray(a)
    flat = a.ravel(order="F")
    if nan:
        flat = flat[np.isfinite(flat)]
    if flat.size == 0:
        return a.dtype.type(0)
    acc = flat[0]
    for v in flat[1:]:
        acc = acc + v
    return acc


def interpol_n(v, n):
    """IDL INTERPOL(V, N) — linear interpolation of V onto N regular points.

    IDL builds the abscissa as FINDGEN(N)*(m-1)/(N-1) in float, then does
    linear interpolation on the original integer grid.
    """
    v = np.asarray(v)
    m = v.size
    if n == m:
        return v.copy()
    x = f32(np.arange(n, dtype=np.float32)) * f32(m - 1) / f32(n - 1)
    return interpol_xu(v, np.arange(m, dtype=np.float64), x)


def interpol_xu(v, x, u):
    """IDL INTERPOL(V, X, U) — linear interpolation, extrapolating at the ends.

    IDL locates each U in X via VALUE_LOCATE (clamped to a valid interval) and
    applies  v[i] + (u-x[i]) * (v[i+1]-v[i])/(x[i+1]-x[i]).  Unlike np.interp it
    extrapolates rather than clamping, so it is written out explicitly.

    X may be monotonically *decreasing* — mcmc_dem's lower confidence bound calls
    this with REVERSE(FINDGEN(n)). VALUE_LOCATE supports that, np.searchsorted
    does not, so the descending case is handled separately.
    """
    v = np.asarray(v, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    u = np.atleast_1d(np.asarray(u, dtype=np.float64))
    m = v.size
    descending = m > 1 and x[0] > x[-1]
    if descending:
        # VALUE_LOCATE on a decreasing vector: largest j with x[j] >= u
        idx = m - 1 - np.searchsorted(x[::-1], u, side="left")
    else:
        idx = np.searchsorted(x, u, side="right") - 1
    idx = np.clip(idx, 0, m - 2)
    x0, x1 = x[idx], x[idx + 1]
    v0, v1 = v[idx], v[idx + 1]
    return v0 + (u - x0) * (v1 - v0) / (x1 - x0)


def histogram(a, mn, mx, binsize):
    """IDL HISTOGRAM(A, MIN=, MAX=, BINSIZE=).

    Bin index is FLOOR((a-min)/binsize); values outside [min,max] are dropped.
    The bin count is FLOOR((max-min)/binsize)+1, and the final bin absorbs
    anything landing exactly on MAX.
    """
    a = np.asarray(a, dtype=np.float64)
    nbin = int(np.floor((mx - mn) / binsize)) + 1
    out = np.zeros(nbin, dtype=np.int64)
    sel = a[(a >= mn) & (a <= mx)]
    if sel.size:
        k = np.floor((sel - mn) / binsize).astype(np.int64)
        k = np.clip(k, 0, nbin - 1)
        np.add.at(out, k, 1)
    return out


def idl_sort(a):
    """IDL SORT — ascending; ties broken by original position (stable)."""
    return np.argsort(np.asarray(a), kind="stable")


def gt_(a, b):
    """IDL `>` operator (elementwise maximum)."""
    return np.maximum(a, b)


def lt_(a, b):
    """IDL `<` operator (elementwise minimum)."""
    return np.minimum(a, b)


def where(cond):
    """IDL WHERE — returns (indices, count)."""
    idx = np.nonzero(np.asarray(cond).ravel(order="F"))[0]
    return idx, idx.size
