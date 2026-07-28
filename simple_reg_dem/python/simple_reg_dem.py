"""simple_reg_dem — Python port of solarsoft_dem/simple_reg_dem/idl/simple_reg_dem.pro.

Regularized DEM inversion (Plowman 2020, arXiv:2006.06828). Faithful 1:1 port of the
IDL solver; verified against the IDL oracle by the IDL2Python-Parity harness.

IDL idioms preserved:
  - element-wise ``>`` / ``<`` are IDL max/min operators -> np.maximum / np.minimum
    (NOT Python comparisons); see analysis/02_construct_report.md items 4,7,9.
  - IDL ``#`` matrix operator: 2D#2D / 2D#1D -> ``@`` ; 1D#1D -> np.outer
    (confirmed by controlled experiment against the oracle).
  - ``shift(diag_matrix(v), s)`` -> ``np.roll(np.diag(v), s, axis=0)``.
  - la_choldc/la_cholsol -> scipy cho_factor/cho_solve; la_choldc status!=0
    (not positive-definite) -> break the pixel loop (LinAlgError).
  - output dems/chi2 are float32 (IDL fltarr); internal math is float64 (IDL promotes
    via the 2.0/6.0 literals and double intermediates).

chk_dump twin-probes mirror staging/simple_reg_dem_probed.pro (no-op unless CHK_DIR set).
"""
import os

import numpy as np
from scipy.linalg import cho_factor, cho_solve, LinAlgError

from chk_dump import chk_dump


def simple_reg_dem(data, errors, exptimes, logt, tresps, *,
                   kmax=100, kcon=5, steps=(0.1, 0.5), drv_con=8.0,
                   chi2_th=1.0, tol=0.1):
    """Return (dems, chi2). dems: (nx,ny,nt) float32; chi2: (nx,ny) float32.

    data, errors: (nx,ny,nchan). exptimes: (nchan,). logt: (nt,). tresps: (nt,nchan).
    """
    i0 = int(os.environ.get("PROBE_I0", "0") or 0)
    j0 = int(os.environ.get("PROBE_J0", "0") or 0)

    data = np.asarray(data)
    errors = np.asarray(errors)
    exptimes = np.asarray(exptimes, dtype=np.float64)
    logt = np.asarray(logt, dtype=np.float64)
    tresps = np.asarray(tresps, dtype=np.float64)
    steps = np.asarray(steps, dtype=np.float64)

    chk_dump("00_input", data=data, errors=errors, exptimes=exptimes,
             logt=logt, tresps=tresps)

    nt = logt.size
    nx = data.shape[0]
    ny = data.shape[1]
    dT = logt[1:nt] - logt[0:nt - 1]
    d0 = np.concatenate(([0.0], dT))          # [0, dT]
    d1 = np.concatenate((dT, [0.0]))          # [dT, 0]
    Bij = (np.diag(d0 + d1) * 2.0
           + np.roll(np.diag(d0), -1, axis=0)
           + np.roll(np.diag(d1), 1, axis=0)) / 6.0
    W = tresps * exptimes[None, :]            # (nt, nchan)
    Rij = W.T @ Bij                           # (nchan, nt)
    e0 = np.concatenate(([0.0], 1.0 / dT))
    e1 = np.concatenate((1.0 / dT, [0.0]))
    Dij = (np.diag(e0 + e1)
           - np.roll(np.diag(e0), -1, axis=0)
           - np.roll(np.diag(e1), 1, axis=0))
    regmat = Dij * exptimes.size / (drv_con ** 2 * (logt[nt - 1] - logt[0]))
    rvec = np.sum(Rij, axis=1)                # (nchan,)   IDL total(Rij,2)

    chk_dump("01_matrices", bij=Bij, rij=Rij, dij=Dij, regmat=regmat, rvec=rvec)

    dems = np.zeros((nx, ny, nt), dtype=np.float32)
    chi2 = np.zeros((nx, ny), dtype=np.float32) - 1.0
    for i in range(nx):
        for j in range(ny):
            err = np.asarray(errors[i, j, :], dtype=np.float64)
            dat0 = np.maximum(np.asarray(data[i, j, :], dtype=np.float64), 0.0)
            num = np.sum(rvec * (np.maximum(dat0, 1.0e-2) / err ** 2))
            den = np.sum((rvec / err) ** 2)
            s = np.log(np.full(nt, num / den))
            if i == i0 and j == j0:
                chk_dump("02a_init", s=s)
            for k in range(kmax):
                dat = (dat0 - Rij @ ((1.0 - s) * np.exp(s))) / err
                mmat = Rij * np.outer(1.0 / err, np.exp(s))        # (nchan, nt)
                amat = mmat.T @ mmat + regmat
                amat_asm = amat.copy() if (i == i0 and j == j0 and k == 0) else None
                try:
                    cf = cho_factor(amat)
                    stat = 0
                except LinAlgError:
                    stat = 1
                if stat == 0:
                    c2p = np.mean((dat0 - Rij @ np.exp(s)) ** 2.0 / err ** 2)
                    deltas = cho_solve(cf, mmat.T @ dat) - s
                    maxad = np.max(np.abs(deltas))
                    deltas = deltas * (np.minimum(maxad, 0.5 / steps[0]) / maxad)
                    ds = 1 - 2 * (1 if c2p < chi2_th else 0)
                    c20 = np.mean((dat0 - Rij @ np.exp(s + deltas * ds * steps[0])) ** 2.0 / err ** 2)
                    c21 = np.mean((dat0 - Rij @ np.exp(s + deltas * ds * steps[1])) ** 2.0 / err ** 2)
                    interp_step = ((steps[0] * (c21 - chi2_th) + steps[1] * (chi2_th - c20)) / (c21 - c20))
                    s = s + deltas * ds * np.minimum(np.maximum(interp_step, steps[0]), steps[1])
                    chi2[i, j] = np.mean((dat0 - Rij @ np.exp(s)) ** 2.0 / err ** 2)
                    if i == i0 and j == j0 and k == 0:
                        chi2_ij = chi2[i, j]      # float32, mirrors IDL chi2[i,j]
                        chk_dump("02b_iter0a", dat=dat, mmat=mmat, amat_asm=amat_asm,
                                 c2p=c2p, ds=ds)
                        chk_dump("02c_iter0b", deltas=deltas, c20=c20, c21=c21,
                                 interp_step=interp_step, s=s, chi2_ij=chi2_ij)
                else:
                    break
                cond1 = (ds * (c2p - c20) / steps[0] < tol) and (k > kcon)
                cond2 = abs(chi2[i, j] - chi2_th) < tol
                if cond1 or cond2:
                    break
            dems[i, j, :] = np.exp(s).astype(np.float32)
            if i == i0 and j == j0:
                chi2_ij = chi2[i, j]              # float32, mirrors IDL chi2[i,j]
                chk_dump("03_final", s=s, chi2_ij=chi2_ij)
    chk_dump("99_final", dems=dems, chi2=chi2)
    return dems, chi2
