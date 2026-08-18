"""dem_sites — Python port of dem_sites.pro (SITES DEM, Morgan & Pickering 2019).

Iterative deconvolution DEM inversion for a single multi-channel pixel. Faithful 1:1
port; verified against the IDL oracle by the IDL2Python-Parity harness.

IDL idioms preserved:
  - element-wise ``>`` are IDL max operators -> np.maximum (NOT comparisons).
  - ``rebin(reform(v,1,nwl),nt,nwl)`` row-broadcast -> v[None,:]; ``rebin(v_nt,nt,nwl)``
    col-broadcast -> v[:,None].
  - ``total(x,2)`` (sum over IDL dim 2 = nwl) -> np.sum(axis=1); ``total(x,1)`` -> axis=0.
  - ``convol(a,ker,/edge_zero)`` -> np.convolve(a,ker,'same') (symmetric kernel, zero pad).
  - ``gaussian_function(sigma,/norm)`` reimplemented from the IDL 8.6 library
    source: width = 2*((ceil(3*sigma)) OR 1)+1, float32 chain, glibc expf.
    (2026-08-18 fix: the odd-forcing ``OR 1`` step was previously missing.)
  - arrays are [nt,nwl] (IDL logical), demmain float; internal math float64.
"""
import numpy as np

from chk_dump import chk_dump


try:                                # IDL float32 exp == glibc expf (1 ulp vs
    import ctypes                   # numpy's SIMD exp) — same approach as the
    _libm = ctypes.CDLL("libm.so.6")  # pintofale_mcmc port's idl_compat.
    _libm.expf.restype = ctypes.c_float
    _libm.expf.argtypes = [ctypes.c_float]

    def _expf(arr):
        return np.array([_libm.expf(np.float32(v)) for v in arr],
                        dtype=np.float32)
except (OSError, AttributeError):   # non-glibc fallback: 1-ulp differences
    def _expf(arr):
        return np.exp(np.asarray(arr, dtype=np.float32)).astype(np.float32)


def gaussian_function(sigma):
    """IDL 8.6 lib/gaussian_function.pro, 1-D, /NORMALIZE — float32 faithful.

    Width rule (verified against the IDL library source): w = CEIL(3*sigma);
    w OR= 1 (force odd — this step was missing in the previous port and made
    the kernel 21 instead of 23 whenever ceil(3*sigma) is even); size = 2*w+1.
    Values: exp(-x^2/(2*sigma^2)) in float32, normalized by the float32 total
    (TOTAL(..., /PRESERVE_TYPE)).
    """
    f32 = np.float32
    s = f32(max(float(sigma), 0.00001))
    w = int(np.ceil(f32(s * f32(3))))
    w |= 1
    n = 2 * w + 1
    x = np.arange(n, dtype=f32) - f32(n // 2)
    a1 = ((x * x) / (f32(2) * (s * s))).astype(f32)
    g = _expf(-a1)
    t = f32(0.0)
    for v in g:
        t = f32(t + v)
    return (g / t).astype(f32)


def dem_sites(obs_in, err_in, response, response_err, delta_temp,
              convergence=1e-2, ker=None):
    """Return (demmain, demerr, obsmod, irep). response is [nt,nwl]."""
    obs_in = np.asarray(obs_in, dtype=np.float64)
    err_in = np.asarray(err_in, dtype=np.float64)
    response = np.asarray(response, dtype=np.float64)
    response_err = np.asarray(response_err, dtype=np.float64)
    delta_temp = np.asarray(delta_temp, dtype=np.float64)
    chk_dump("00_input", obs_in=obs_in, err_in=err_in, response=response,
             response_err=response_err, delta_temp=delta_temp)

    nt, nwl = response.shape
    wt = 1.0 / np.sqrt((err_in / obs_in) ** 2 + response_err ** 2)
    num = response * wt[None, :]
    gres = num / np.sum(num, axis=1)[:, None]

    obs = obs_in.copy()
    # IDL: nker=(0.08*nt)>0.5 — float32 literal chain preserved
    nker = max(np.float32(np.float32(0.08) * np.float32(nt)), np.float32(0.5))
    if ker is None:
        ker = gaussian_function(nker)
    ker = np.asarray(ker)

    irep = 0
    demmain = np.zeros(nt)                       # IDL fltarr(nt)
    res2 = response * delta_temp[:, None]
    totres2 = np.sum(res2 ** 2, axis=0)
    res3 = res2 * gres
    chk_dump("01_setup", wt=wt, gres=gres, ker=ker, res2=res2, totres2=totres2, res3=res3)

    e = (err_in / obs_in) ** 2 + response_err ** 2
    demerr = np.sqrt(np.sum(e[None, :] * gres, axis=1))
    totwt = np.sum(wt)
    obswt = wt / obs_in

    while True:
        dem = (obs / totres2)[None, :] * res3
        conv = np.convolve(np.sum(dem, axis=1), ker, mode="same")
        demmain = np.maximum(demmain + conv, 0.0)
        obs = obs_in - np.sum(demmain[:, None] * res2, axis=0)
        convobs = np.sum(np.abs(obs * obswt)) / totwt
        irep += 1
        if irep == 1:
            chk_dump("02_iter0", dem=dem, demmain=demmain, obs=obs,
                     convobs=np.float64(convobs))
        if irep > 300 or convobs < convergence:
            break

    obsmod = np.sum((demmain * delta_temp)[:, None] * response, axis=0)
    demerr = demerr * demmain
    chk_dump("99_final", demmain=demmain, demerr=demerr, obsmod=obsmod, irep=np.int32(irep))
    return demmain, demerr, obsmod, irep
