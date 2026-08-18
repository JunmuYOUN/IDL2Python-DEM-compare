"""aia_sparse_em — Python port of aia_sparse_em_init.pro (Cheung et al. 2015 sparse-EM DEM).

Basis-pursuit DEM: minimise L1(EM coefficients) s.t. the reconstructed DN lies within a
tolerance band of the observed DN, coefficients >= 0. The IDL core uses IDL's built-in
`simplex` (Numerical Recipes simplx); this port uses scipy.optimize.linprog (HiGHS).

Faithful port of the deterministic setup (basis functions, dictionary D=K.B). The LP solve
is mathematically equivalent but the *vertex* selected among degenerate optima can differ
between simplex implementations — see the parity certificate.
"""
import numpy as np
from scipy.optimize import linprog

from chk_dump import chk_dump


def aia_sparse_em_init(tresp, lgtaxis, bases_sigmas=(0.0, 0.1, 0.2), dictfac=1e26):
    """Build basis_funcs [(nsig*nt), nt] and Dict [nchan, nbasis] = (K.B)*dictfac."""
    tresp = np.asarray(tresp, dtype=np.float64)          # (nchan, nt)
    lgtaxis = np.asarray(lgtaxis, dtype=np.float64)
    nchannels, ntemp = tresp.shape
    nsigmas = len(bases_sigmas)
    basis_funcs = np.zeros((nsigmas * ntemp, ntemp))
    for s, sig in enumerate(bases_sigmas):
        if sig == 0.0:
            for m in range(ntemp):
                basis_funcs[ntemp * s + m, m] = 1.0      # Dirac delta bases
        else:
            for m in range(ntemp):
                line = np.exp(-((lgtaxis - lgtaxis[m]) / sig) ** 2.0)
                line[line < 0.04] = 0.0
                basis_funcs[ntemp * s + m, :] = line
    Dict = (tresp @ basis_funcs.T) * dictfac             # (nchan, nbasis)
    chk_dump("01_init", lgtaxis=lgtaxis, basis_funcs=basis_funcs, tresp=tresp,
             dict=Dict, dictfac=np.float32(dictfac))
    return Dict, basis_funcs


def aia_sparse_em_solve(image, Dict, basis_funcs, tolfac=1.4, eps=1e-3, symmbuff=1.0):
    """Per-pixel sparse-EM LP. Returns (coeffs, oem, zmax, status). image is [nx,ny,nchan]."""
    image = np.asarray(image, dtype=np.float64)
    nx, ny, nchannels = image.shape
    nbasis = Dict.shape[1]
    nt = basis_funcs.shape[1]
    coeffs = np.zeros((nx, ny, nbasis))
    zmax = np.zeros((nx, ny))
    status = np.zeros((nx, ny))
    chk_dump("00_input", image=image, eps=np.float32(eps), tolfac=np.float32(tolfac))

    c = np.ones(nbasis)                                  # minimise sum(coeffs) = L1
    A_ub = np.vstack([Dict, -Dict])
    for i in range(nx):
        for j in range(ny):
            y = np.maximum(image[i, j, :], 0.0)
            tol = tolfac * np.sqrt(y)                     # tolfunc = 'sqrt(y)'
            b_ub = np.concatenate([y + tol, -np.maximum(y - symmbuff * tol, 0.0)])
            res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
            if res.success:
                x = np.asarray(res.x)
                s = 0
            else:
                x = np.zeros(nbasis)
                s = 1
            if x.min() < 0.0:                            # simplex-parity: negative -> zero, status 10
                x = np.zeros(nbasis)
                s = 10
            coeffs[i, j, :] = x
            zmax[i, j] = -np.sum(x)                       # IDL r[0] = max(-sum) = -sum(coeffs)
            status[i, j] = s
            if i == 0 and j == 0:
                chk_dump("02_lp", y=y, tol=tol, r=np.concatenate([[zmax[i, j]], x]))

    oem = np.zeros((nx, ny, nt))
    for i in range(nx):
        for j in range(ny):
            oem[i, j, :] = coeffs[i, j, :] @ basis_funcs
    chk_dump("99_solve", coeffs=coeffs, oem=oem,
             zmax=np.atleast_1d(zmax.ravel()), status=np.atleast_1d(status.ravel()))
    return coeffs, oem, zmax, status
