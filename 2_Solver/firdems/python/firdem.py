"""firdem — Python port of the verifiable deterministic core of firdem.pro
(Plowman 2012 fast/iterative/regularized DEM).

Scope (oracle-verified): the first-pass REGULARIZED inversion —
  * firdem_regularize_data : iterative Tikhonov-style data regularization
  * first_pass_dem         : DEM coefficients = basis22 # (a_inv # (data/normfac))
  * firdem_triangle_basis  : higher-resolution triangle basis (float32-grid limited)

Injected library boundary (not re-derived here): a_struc (response spline interpolation
via spl_interp, basis normalisation and a_array via int_tabulated, SVD via svdc/column ->
a_inv), basis22 (spl_interp), a2_array (int_tabulated), chisqr_cvf. The 2000-iteration
negative-EM removal (firdem_iterate) is path-dependent and OUT of strict-parity scope.

IDL idioms: element-wise ``>`` -> np.maximum; ``M # v`` (2D#1D) -> M @ v;
``a # b`` (1D#1D) -> np.outer; invert(/double) -> np.linalg.inv; findgen -> np.arange.
"""
import numpy as np


def firdem_triangle(x, x0, b, h):
    sep = np.abs(x - x0)
    y = h * (0.5 * b - sep) / (0.5 * b)
    y = np.where(y < 0, 0.0, y)
    return y


def firdem_triangle_basis(nb, logt):
    """nb triangle basis functions on grid logt. Returns (basis[nt,nb], t0a[nb])."""
    logt = np.asarray(logt, dtype=np.float64)
    logtlo = np.min(logt)
    logthi = np.max(logt)
    logdt = (logthi - logtlo) / (nb + 1)
    nt = logt.size
    t0a = logdt + logtlo + logdt * np.arange(nb, dtype=np.float64)
    basis = np.zeros((nt, nb))
    for j in range(nb):
        basis[:, j] = firdem_triangle(logt, t0a[j], 2.0 * logdt, 1.0)
    return basis, t0a


def firdem_regularize_data(datavec0, sigmas, tr_norm, a_inv,
                           chi2_end, chi2_tol=0.05, niter_max=50):
    """Port of firdem_regularize_data. a_inv, tr_norm from the (injected) a_struc.
    tr_norm here already includes exptimes (a_struc.tr_norm*a_struc.exptimes)."""
    datavec0 = np.asarray(datavec0, dtype=np.float64)
    sigmas = np.asarray(sigmas, dtype=np.float64)
    tr_norm = np.asarray(tr_norm, dtype=np.float64)
    a_inv = np.asarray(a_inv, dtype=np.float64)
    nchan = datavec0.size

    if np.sum((datavec0 / sigmas) ** 2) < chi2_end:
        return np.zeros(nchan)
    if chi2_end / nchan < 1.0e-4:
        return datavec0

    a_inv_scaled = a_inv / np.outer(tr_norm, tr_norm)
    sigs2_diag = np.diag((1.0 / sigmas) ** 2)

    alpha_low = 0.0
    alpha_high = 0.0
    bisect_start = 0
    alpha = 0.0
    datavec = datavec0 / tr_norm
    for k in range(niter_max):
        datavec = datavec0 / tr_norm
        data2vec = (a_inv @ datavec) / tr_norm
        if k == 0:
            alpha = np.sqrt(9.0 * chi2_end / np.sum(sigmas ** 2 * data2vec ** 2))
        b_inv = np.linalg.inv(alpha * a_inv_scaled + sigs2_diag)
        datavec = datavec - (alpha * (b_inv @ data2vec)) / tr_norm
        chi2 = np.sum(((datavec0 - datavec * tr_norm) / sigmas) ** 2)

        if bisect_start == 1:
            if chi2 < chi2_end:
                alpha_low = alpha
            if chi2 > chi2_end:
                alpha_high = alpha
            alpha = alpha_low + 0.5 * (alpha_high - alpha_low)
        if chi2 < chi2_end and bisect_start == 0:
            alpha = 5 * alpha
        if chi2 > chi2_end and bisect_start == 0:
            bisect_start = 1
            alpha_high = alpha
            alpha = alpha_low + 0.5 * (alpha_high - alpha_low)

        if abs(chi2 - chi2_end) / chi2_end < chi2_tol:
            break

    return datavec * tr_norm


def first_pass_dem(data_out, a_inv, normfac, basis22):
    """DEM coefficients in the high-res basis: basis22 # (a_inv # (data_out/normfac))."""
    a_inv = np.asarray(a_inv, dtype=np.float64)
    normfac = np.asarray(normfac, dtype=np.float64)
    basis22 = np.asarray(basis22, dtype=np.float64)
    return basis22 @ (a_inv @ (np.asarray(data_out, dtype=np.float64) / normfac))
