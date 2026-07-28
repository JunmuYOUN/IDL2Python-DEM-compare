"""dem_sites — Python port of dem_sites.pro (SITES DEM, Morgan & Pickering 2019).

Iterative deconvolution DEM inversion for a single multi-channel pixel. Faithful 1:1
port; verified against the IDL oracle by the IDL2Python-Parity harness.

IDL idioms preserved:
  - element-wise ``>`` are IDL max operators -> np.maximum (NOT comparisons).
  - ``rebin(reform(v,1,nwl),nt,nwl)`` row-broadcast -> v[None,:]; ``rebin(v_nt,nt,nwl)``
    col-broadcast -> v[:,None].
  - ``total(x,2)`` (sum over IDL dim 2 = nwl) -> np.sum(axis=1); ``total(x,1)`` -> axis=0.
  - ``convol(a,ker,/edge_zero)`` -> np.convolve(a,ker,'same') (symmetric kernel, zero pad).
  - ``gaussian_function(sigma,/norm)`` reimplemented: size 2*ceil(3*sigma)+1, float32.
  - arrays are [nt,nwl] (IDL logical), demmain float; internal math float64.
"""
import numpy as np

from chk_dump import chk_dump


def gaussian_function(sigma):
    """IDL gaussian_function(sigma,/normalize): centered Gaussian, size 2*ceil(3*sigma)+1."""
    n = int(2 * np.ceil(3 * sigma) + 1)
    x = np.arange(n) - (n - 1) / 2.0
    g = np.exp(-x ** 2 / (2.0 * sigma ** 2))
    return (g / np.sum(g)).astype(np.float32)


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
    nker = max(0.08 * nt, 0.5)
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
