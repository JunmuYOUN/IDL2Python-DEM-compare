"""aia_teem — Python port of the aschwanden aia_teem DEM core (Aschwanden 2011/2013).

Two pieces, ported verbatim from aia_teem_table.pro and the fit loop of aia_teem_map.pro
(lines 124-158, nx<4096 all-pixels branch):
  - build_flux_table: AIA single-Gaussian flux lookup table (forward model).
  - fit: per-pixel brute-force (Te, sigma) grid search minimising reduced chi (deterministic
    — the source uses nested loops, NOT an optimizer, so parity is exact modulo float32).

Excluded (I/O scaffolding, not the DEM algorithm): FITS read, aia_prep, make_map, plotting.
Verified against the IDL oracle by the IDL2Python-Parity harness.
"""
import numpy as np

from chk_dump import chk_dump


def build_flux_table(resp, telog, dte, tsig):
    """flux[nte,nsig,nwave] for unity-EM single-Gaussian DEMs. resp is [nte,nwave]."""
    resp = np.asarray(resp, dtype=np.float64)
    telog = np.asarray(telog, dtype=np.float64)
    dte = np.asarray(dte, dtype=np.float64)
    tsig = np.asarray(tsig, dtype=np.float64)
    nte, nwave = resp.shape
    nsig = tsig.size
    flux = np.zeros((nte, nsig, nwave))
    for i in range(nte):
        for j in range(nsig):
            em_kelvin = np.exp(-(telog - telog[i]) ** 2 / (2.0 * tsig[j] ** 2))
            for iw in range(nwave):
                flux[i, j, iw] = np.sum(resp[:, iw] * em_kelvin * dte)
    return flux


def fit(images, texp_, flux, telog, tsig):
    """Per-pixel single-Gaussian T/EM fit. Returns (te_map, em_map, sig_map, chi_map).

    Kept in float32 (IDL fltarr): the reduced-chi residual (flux_obs - flux_dem) is a small
    difference of large numbers near the best fit, so float64 would diverge from IDL's float32
    by ~1e-3 via cancellation. float32 reproduces IDL exactly.
    """
    images = np.asarray(images, dtype=np.float32)
    texp_ = np.asarray(texp_, dtype=np.float32)
    flux = np.asarray(flux, dtype=np.float32)
    telog = np.asarray(telog, dtype=np.float32)
    tsig = np.asarray(tsig, dtype=np.float32)
    nte, nsig, nwave = flux.shape
    nxx, nyy = images.shape[0], images.shape[1]
    nfree = 3
    te_map = np.zeros((nxx, nyy), dtype=np.float32); em_map = np.zeros((nxx, nyy), dtype=np.float32)
    sig_map = np.zeros((nxx, nyy), dtype=np.float32); chi_map = np.zeros((nxx, nyy), dtype=np.float32)
    for j in range(nyy):
        for i in range(nxx):
            flux_obs = images[i, j, :]
            counts = flux_obs * texp_
            noise = np.sqrt(counts) / texp_
            chimin = 9999.0
            em_best = te_best = sig_best = 0.0
            for k in range(nte):
                for l in range(nsig):
                    flux_dem1 = flux[k, l, :]
                    em1 = np.sum(flux_obs) / np.sum(flux_dem1)
                    flux_dem = flux_dem1 * em1
                    chi = np.sqrt(np.sum((flux_obs - flux_dem) ** 2 / noise ** 2) / (nwave - nfree))
                    if chi <= chimin:              # IDL 'le': ties keep the later (k,l)
                        chimin = chi
                        em_best = np.log10(em1)
                        te_best = telog[k]
                        sig_best = tsig[l]
            em_map[i, j] = em_best; te_map[i, j] = te_best
            sig_map[i, j] = sig_best; chi_map[i, j] = chimin
            if i == 0 and j == 0:
                chk_dump("02_pixel", flux_obs=flux_obs, noise=noise, te_best=te_best,
                         sig_best=sig_best, em_best=em_best, chimin=chimin)
    return te_map, em_map, sig_map, chi_map
