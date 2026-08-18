"""aia_sparse_em — Python port of aia_sparse_em_init.pro (Cheung et al. 2015 sparse-EM DEM).

Basis-pursuit DEM: minimise L1(EM coefficients) s.t. the reconstructed DN lies within a
tolerance band of the observed DN, coefficients >= 0.

2026-08-18 engine change: the default LP engine is now a faithful float32
reimplementation of IDL's built-in `simplex` (Numerical Recipes simplx, per the
IDL documentation), including its EPS semantics and non-convergence behaviour
(status=3 with an all-zero result). This reproduces the IDL original's numbers
on real data — including its failures. The previous scipy.optimize.linprog
(HiGHS) engine is kept as `engine="highs"`; it is more robust (it solves pixels
IDL's simplex cannot) but is NOT numerically parity-equivalent to the IDL code.

Status values (IDL SIMPLEX convention):
  0 success · 1 unbounded · 2 infeasible · 3 did not converge ·
  10 negative coefficient zeroed (AIA_SPARSE_EM_SOLVE) · 11 no counts.

Faithful-arithmetic notes:
  - basis_funcs is float32 (IDL fltarr); Dict = double(tresp) @ double(f32 basis).
  - the constraint tableau is float32 (IDL `constraint = fltarr(m,ntemp+1)`).
  - tolfac/eps are rounded through float32 before use, because the IDL callers
    pass single-precision literals (1.4, 1e-3).
"""
import numpy as np

from chk_dump import chk_dump

_F = np.float32


def aia_sparse_em_init(tresp, lgtaxis, bases_sigmas=(0.0, 0.1, 0.2), dictfac=1e26):
    """Build basis_funcs [(nsig*nt), nt] (float32) and Dict [nchan, nbasis]."""
    tresp = np.asarray(tresp, dtype=np.float64)          # (nchan, nt)
    lgtaxis = np.asarray(lgtaxis, dtype=np.float64)
    nchannels, ntemp = tresp.shape
    nsigmas = len(bases_sigmas)
    basis_funcs = np.zeros((nsigmas * ntemp, ntemp), dtype=_F)   # IDL fltarr
    for s, sig in enumerate(bases_sigmas):
        if sig == 0.0:
            for m in range(ntemp):
                basis_funcs[ntemp * s + m, m] = 1.0      # Dirac delta bases
        else:
            for m in range(ntemp):
                line = np.exp(-((lgtaxis - lgtaxis[m]) / sig) ** 2.0)
                line[line < 0.04] = 0.0
                basis_funcs[ntemp * s + m, :] = line     # f32 store, like IDL
    Dict = (tresp @ basis_funcs.astype(np.float64).T) * dictfac  # (nchan, nbasis)
    chk_dump("01_init", lgtaxis=lgtaxis, basis_funcs=basis_funcs, tresp=tresp,
             dict=Dict, dictfac=np.float32(dictfac))
    return Dict, basis_funcs


# ---------------------------------------------------------------------------
# IDL SIMPLEX (Numerical Recipes in C, 2nd ed, §10.8 simplx) — float32.
# ---------------------------------------------------------------------------

def _idl_simplex(zeq, rows, m1, m2, m3, eps, itmax=5000):
    """NR simplx. rows = (m, n+1) float32 [b_i, -A_i*] with b_i >= 0,
    ordered <=, >=, ==. Returns (r[(n+1,)], status, npiv); on non-convergence
    returns zeros with status 3, matching IDL empirically."""
    eps = _F(eps)
    zero = _F(0.0)
    m = m1 + m2 + m3
    n = zeq.shape[0]
    a = np.zeros((m + 3, n + 2), dtype=_F)               # NR 1-based tableau
    a[1, 2:n + 2] = zeq
    a[2:m + 2, 1:n + 2] = rows

    izrov = np.arange(0, n + 1)
    iposv = np.zeros(m + 1, dtype=np.int64)
    for i in range(1, m + 1):
        if a[i + 1, 1] < zero:
            raise ValueError("bad input tableau: negative b")
        iposv[i] = n + i
    l1 = list(range(1, n + 1))
    npiv = 0

    def simp1(mm, iabf):
        if not l1:
            return 0, zero
        cols = np.asarray(l1, dtype=np.int64)
        vals = a[mm + 1, cols + 1]
        if iabf == 0:
            idx = int(np.argmax(vals))
            return l1[idx], vals[idx]
        idx = int(np.argmax(np.abs(vals)))
        return l1[idx], vals[idx]

    def simp2(kp):
        i = 1
        while i <= m:
            if a[i + 1, kp + 1] < -eps:
                break
            i += 1
        if i > m:
            return 0
        q1 = _F(-a[i + 1, 1] / a[i + 1, kp + 1])
        ip = i
        i += 1
        while i <= m:
            if a[i + 1, kp + 1] < -eps:
                q = _F(-a[i + 1, 1] / a[i + 1, kp + 1])
                if q < q1:
                    ip = i
                    q1 = q
                elif q == q1:                            # degenerate tie-break
                    q0 = qp = zero
                    for k in range(1, n + 1):
                        qp = _F(-a[ip + 1, k + 1] / a[ip + 1, kp + 1])
                        q0 = _F(-a[i + 1, k + 1] / a[i + 1, kp + 1])
                        if q0 != qp:
                            break
                    if q0 < qp:
                        ip = i
            i += 1
        return ip

    def simp3(i1, k1, ip, kp):
        piv = _F(_F(1.0) / a[ip + 1, kp + 1])
        for ii in range(1, i1 + 2):
            if ii - 1 != ip:
                a[ii, kp + 1] = _F(a[ii, kp + 1] * piv)
                coef = a[ii, kp + 1]
                save = coef
                a[ii, 1:k1 + 2] = (a[ii, 1:k1 + 2]
                                   - a[ip + 1, 1:k1 + 2] * coef).astype(_F)
                a[ii, kp + 1] = save
        a[ip + 1, 1:k1 + 2] = (a[ip + 1, 1:k1 + 2] * (-piv)).astype(_F)
        a[ip + 1, kp + 1] = piv

    def fail(code):
        return np.zeros(n + 1, dtype=_F), code, npiv

    if m2 + m3:
        l3 = [1] * (m2 + 1)
        for k in range(1, n + 2):
            q1 = zero
            for i in range(m1 + 1, m + 1):
                q1 = _F(q1 + a[i + 1, k])
            a[m + 2, k] = -q1
        while True:
            kp, bmax = simp1(m + 1, 0)
            if bmax <= eps and a[m + 2, 1] < -eps:
                return fail(2)                           # infeasible
            if bmax <= eps and a[m + 2, 1] <= eps:
                skip_to_one = False
                for ip_ in range(m1 + m2 + 1, m + 1):
                    if iposv[ip_] == ip_ + n:
                        kp, bmax = simp1(ip_, 1)
                        if bmax > eps:
                            ip = ip_
                            skip_to_one = True
                            break
                if not skip_to_one:
                    for i in range(m1 + 1, m1 + m2 + 1):
                        if l3[i - m1] == 1:
                            a[i + 1, 1:n + 2] = -a[i + 1, 1:n + 2]
                    break
            else:
                ip = simp2(kp)
                if ip == 0:
                    return fail(2)                       # infeasible
            simp3(m + 1, n, ip, kp)
            npiv += 1
            if npiv > itmax:
                return fail(3)                           # did not converge
            if iposv[ip] >= n + m1 + m2 + 1:
                l1.remove(kp)
            else:
                kh = iposv[ip] - m1 - n
                if kh >= 1 and l3[kh]:
                    l3[kh] = 0
                    a[m + 2, kp + 1] = _F(a[m + 2, kp + 1] + _F(1.0))
                    a[1:m + 3, kp + 1] = -a[1:m + 3, kp + 1]
            izrov[kp], iposv[ip] = iposv[ip], izrov[kp]

    while True:
        kp, bmax = simp1(0, 0)
        if bmax <= eps:
            break                                        # optimum reached
        ip = simp2(kp)
        if ip == 0:
            return fail(1)                               # unbounded
        simp3(m, n, ip, kp)
        npiv += 1
        if npiv > itmax:
            return fail(3)
        izrov[kp], iposv[ip] = iposv[ip], izrov[kp]

    r = np.zeros(n + 1, dtype=_F)
    r[0] = a[1, 1]
    for i in range(1, m + 1):
        if iposv[i] <= n:
            r[iposv[i]] = a[i + 1, 1]
    return r, 0, npiv


def aia_sparse_em_solve(image, Dict, basis_funcs, tolfac=1.4, eps=1e-3,
                        symmbuff=1.0, engine="idl", itmax=5000):
    """Per-pixel sparse-EM LP. Returns (coeffs, oem, zmax, status). image is [nx,ny,nchan]."""
    image = np.asarray(image, dtype=np.float64)
    nx, ny, nchannels = image.shape
    nbasis = Dict.shape[1]
    nt = basis_funcs.shape[1]
    coeffs = np.zeros((nx, ny, nbasis), dtype=_F)        # IDL fltarr
    zmax = np.zeros((nx, ny), dtype=_F)
    status = np.zeros((nx, ny), dtype=_F)
    chk_dump("00_input", image=image, eps=np.float32(eps), tolfac=np.float32(tolfac))

    # IDL callers pass single-precision literals
    tolfac32 = float(np.float32(tolfac))
    eps32 = float(np.float32(eps))
    e8 = float(np.float32(8e-4))
    nocounts = image.sum(axis=2) < 10.0 * eps32          # IDL NOCOUNTS

    if engine == "idl":
        zeq = np.full(nbasis, -1.0, dtype=_F)            # maximize -(L1 norm)
        negdict32 = (-Dict).astype(_F)                   # constraint cols, set once
        rows = np.empty((2 * nchannels, nbasis + 1), dtype=_F)
        rows[:nchannels, 1:] = negdict32
        rows[nchannels:, 1:] = negdict32
        for i in range(nx):
            for j in range(ny):
                y = np.maximum(image[i, j, :], 0.0)
                tol = tolfac32 * np.sqrt(y)              # double, like IDL
                rows[:nchannels, 0] = (y + tol).astype(_F)
                rows[nchannels:, 0] = np.maximum(y - symmbuff * tol, 0.0).astype(_F)
                r, s, _ = _idl_simplex(zeq, rows, nchannels, nchannels, 0,
                                       _F(eps32 * y.max() * e8), itmax=itmax)
                x = r[1:]
                if x.min() < 0.0:                        # IDL: negative -> zero, 10
                    x = np.zeros(nbasis, dtype=_F)
                    s = 10
                coeffs[i, j, :] = x
                zmax[i, j] = r[0]                        # IDL r[0]
                status[i, j] = s
                if i == 0 and j == 0:
                    chk_dump("02_lp", y=y, tol=tol,
                             r=np.concatenate([[zmax[i, j]], x]))
    elif engine == "highs":
        from scipy.optimize import linprog
        c = np.ones(nbasis)
        A_ub = np.vstack([Dict, -Dict])
        for i in range(nx):
            for j in range(ny):
                y = np.maximum(image[i, j, :], 0.0)
                tol = tolfac * np.sqrt(y)
                b_ub = np.concatenate([y + tol,
                                       -np.maximum(y - symmbuff * tol, 0.0)])
                res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None),
                              method="highs")
                if res.success:
                    x = np.asarray(res.x)
                    s = 0
                else:
                    x = np.zeros(nbasis)
                    s = 1
                if x.min() < 0.0:
                    x = np.zeros(nbasis)
                    s = 10
                coeffs[i, j, :] = x
                zmax[i, j] = -np.sum(x)
                status[i, j] = s
    else:
        raise ValueError("engine must be 'idl' or 'highs'")

    if np.any(nocounts):                                 # IDL: status 11
        status[nocounts] = 11.0
        coeffs[nocounts, :] = 0.0

    oem = np.zeros((nx, ny, nt), dtype=_F)               # IDL fltarr, f32 '#'
    for i in range(nx):
        for j in range(ny):
            oem[i, j, :] = coeffs[i, j, :] @ basis_funcs
    chk_dump("99_solve", coeffs=coeffs, oem=oem,
             zmax=np.atleast_1d(zmax.ravel()), status=np.atleast_1d(status.ravel()))
    return coeffs, oem, zmax, status
