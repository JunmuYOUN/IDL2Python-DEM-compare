"""chianti_dem_fit — Python port of chianti_dem/dem_fit.pro (classic CHIANTI spline-DEM
Levenberg-Marquardt fit, Del Zanna/Dere). Faithful float32 port; verified vs IDL oracle.
Includes an exact port of IDL SPLINE (spline-under-tension, sigma=1.0). Plotting removed."""
import numpy as np
from chk_dump import chk_dump

f32 = np.float32


def idl_spline(x, y, t, sigma=1.0):
    x = np.asarray(x, f32); y = np.asarray(y, f32); t = np.asarray(t, f32)
    n = min(len(x), len(y)); nm1 = n - 1
    sigma = f32(max(sigma, 0.001)); yp = np.zeros(2 * n, f32)
    delx1 = x[1]-x[0]; dx1 = (y[1]-y[0])/delx1
    delx2 = x[2]-x[1]; delx12 = x[2]-x[0]
    c1 = -(delx12+delx1)/delx12/delx1; c2 = delx12/delx1/delx2; c3 = -delx1/delx12/delx2
    slpp1 = c1*y[0]+c2*y[1]+c3*y[2]
    deln = x[nm1]-x[nm1-1]; delnm1 = x[nm1-1]-x[nm1-2]; delnn = x[nm1]-x[nm1-2]
    c1 = (delnn+deln)/delnn/deln; c2 = -delnn/deln/delnm1; c3 = deln/delnn/delnm1
    slppn = c3*y[nm1-2]+c2*y[nm1-1]+c1*y[nm1]
    sigmap = f32(sigma*nm1/(x[nm1]-x[0]))
    dels = sigmap*delx1; exps = np.exp(dels); sinhs = f32(0.5)*(exps-1./exps)
    sinhin = 1./(delx1*sinhs); diag1 = sinhin*(dels*f32(0.5)*(exps+1./exps)-sinhs)
    diagin = 1./diag1; yp[0] = diagin*(dx1-slpp1)
    spdiag = sinhin*(sinhs-dels); yp[n] = diagin*spdiag
    dxa = x[1:]-x[:-1]; dxv = (y[1:]-y[:-1])/dxa
    dv = sigmap*dxa; ev = np.exp(dv); sv = f32(0.5)*(ev-1./ev); si = 1./(dxa*sv)
    diag2 = si*(dv*(f32(0.5)*(ev+1./ev))-sv)
    diag2 = np.concatenate(([f32(0)], (diag2[:-1]+diag2[1:]))).astype(f32)
    dx2nm1 = dxv[nm1-1]
    dx2 = np.concatenate(([f32(0)], (dxv[1:]-dxv[:-1]))).astype(f32)
    spd = si*(sv-dv)
    for i in range(1, nm1):
        gi = 1./(diag2[i]-spd[i-1]*yp[i+n-1])
        yp[i] = gi*(dx2[i]-spd[i-1]*yp[i-1]); yp[i+n] = gi*spd[i]
    gi = 1./(diag1-spd[nm1-1]*yp[n+nm1-1])
    yp[nm1] = gi*(slppn-dx2nm1-spd[nm1-1]*yp[nm1-1])
    for i in range(n-2, -1, -1):
        yp[i] = yp[i]-yp[i+n]*yp[i+1]
    subs = np.clip(np.searchsorted(x, t, side='right'), 1, nm1).astype(int); subs1 = subs-1
    d1 = t-x[subs1]; d2 = x[subs]-t; ds = x[subs]-x[subs1]
    e1 = np.exp(sigmap*d1); sh1 = f32(0.5)*(e1-1./e1)
    e = np.exp(sigmap*d2); sh2 = f32(0.5)*(e-1./e); e = e1*e; shs = f32(0.5)*(e-1./e)
    spl = (yp[subs]*sh1+yp[subs1]*sh2)/shs + ((y[subs]-yp[subs])*d1+(y[subs1]-yp[subs1])*d2)/ds
    return spl.astype(f32)


class DemFit:
    """State object standing in for the IDL common blocks obs/dem/contr."""
    def __init__(s, obs_int, obs_sig, ch_tot_contr, log_t_mesh, log_dem_temp,
                 dem_temp, d_dem_temp, log_dem_mesh):
        s.obs_int = np.asarray(obs_int, f32); s.obs_sig = np.asarray(obs_sig, f32)
        s.n_obs = s.obs_int.size
        s.ch_tot_contr = np.asarray(ch_tot_contr, f32)        # (n_dem_temp, n_obs)
        s.log_t_mesh = np.asarray(log_t_mesh, f32)
        s.log_dem_temp = np.asarray(log_dem_temp, f32)
        s.dem_temp = np.asarray(dem_temp, f32)
        s.d_dem_temp = f32(np.asarray(d_dem_temp).reshape(-1)[0])
        s.log_dem_mesh = np.asarray(log_dem_mesh, f32).copy()

    def functn(s):
        dem = idl_spline(s.log_t_mesh, s.log_dem_mesh, s.log_dem_temp)
        dem = f32(10.0) ** dem
        dlnt = np.log(f32(10.0) ** s.d_dem_temp)
        out = np.zeros(s.n_obs, f32)
        for i in range(s.n_obs):
            out[i] = np.sum(s.ch_tot_contr[:, i] * dem * s.dem_temp) * dlnt
        return out

    def deriv(s):
        nt = s.log_dem_mesh.size
        d = np.zeros((s.n_obs, nt), f32)
        sav = s.log_dem_mesh.copy()
        da = np.maximum(f32(0.01) * s.log_dem_mesh, f32(1e-6)).astype(f32)
        yfit = s.functn()
        for j in range(nt):
            s.log_dem_mesh[j] = sav[j] + da[j]
            d[:, j] = (s.functn() - yfit) / da[j]
            s.log_dem_mesh[j] = sav[j]
        return d

    def chisqr(s):
        dy = s.obs_int - s.functn()
        return np.float32(np.sum(dy ** 2 / s.obs_sig ** 2))

    def fit(s, flambda=10.0, scale=1.0, niter=20, dchisq_min=1e-5):
        flambda = f32(flambda); scale = f32(scale)
        chk_dump("00_input", obs_int=s.obs_int, obs_sig=s.obs_sig,
                 ch_tot_contr=s.ch_tot_contr, log_t_mesh=s.log_t_mesh, log_dem_temp=s.log_dem_temp)
        chk_dump("00b_grids", dem_temp=s.dem_temp, d_dem_temp=s.d_dem_temp, log_dem_mesh=s.log_dem_mesh)
        yinit = s.functn(); chk_dump("01_forward", yinit=yinit)
        chisqr = s.chisqr(); chisq1 = chisqr; dchisq = f32(1.0)
        y = s.functn(); weight = 1.0 / np.maximum(s.obs_sig, f32(1.0))
        nterms = s.log_t_mesh.size; nfree = s.n_obs - nterms
        failed = 0
        if nfree <= 0:
            failed = 1
        it = 0
        while (it <= niter) and (dchisq >= dchisq_min) and not failed:
            alpha = np.zeros((nterms, nterms), f32); beeta = np.zeros(nterms, f32)
            b = s.log_dem_mesh.copy()
            deriv1 = s.deriv()
            dy = s.obs_int - s.functn()
            for i in range(s.n_obs):
                for j in range(nterms):
                    beeta[j] += weight[i] * dy[i] * deriv1[i, j]
                    for k in range(j + 1):
                        alpha[j, k] += weight[i] * deriv1[i, j] * deriv1[i, k]
            alpha = alpha.T.copy()
            chisq1 = s.chisqr()
            while True:  # start_this_iter (flambda retry)
                array = np.zeros((nterms, nterms), f32)
                for j in range(nterms):
                    for k in range(nterms):
                        array[j, k] = alpha[j, k] / np.sqrt(alpha[j, j] * alpha[k, k])
                    array[j, j] = f32(1.0) + flambda
                array = np.linalg.inv(array).astype(f32)
                for j in range(nterms):
                    s.log_dem_mesh[j] = b[j]
                    for k in range(nterms):
                        s.log_dem_mesh[j] = s.log_dem_mesh[j] + beeta[k] * array[j, k] / np.sqrt(alpha[j, j] * alpha[k, k]) * scale
                chisqr = s.chisqr()
                if np.isnan(chisqr):
                    if flambda > 1e10:
                        failed = 1; break
                    flambda = f32(10.0) * flambda
                    if scale > 0.1: scale = scale / f32(2.0)
                    continue
                if chisqr < chisq1:
                    if flambda == 0:
                        failed = 1; break
                    flambda = flambda / f32(10.0)
                    break
                else:
                    if flambda == 0 or flambda > 1e10:
                        failed = 1; break
                    flambda = f32(10.0) * flambda
                    if scale > 0.1: scale = scale / f32(2.0)
                    continue
            if failed:
                break
            chisqr = s.chisqr()
            dchisq = chisq1 - chisqr; chisq1 = chisqr; it += 1
            if it == 1:
                chk_dump("02_iter1", alpha=alpha, beeta=beeta, log_dem_mesh=s.log_dem_mesh, chisqr=chisqr)
        chisqr = s.chisqr(); y = s.functn()
        chk_dump("99_final", log_dem_mesh=s.log_dem_mesh, y=y, chisqr=chisqr, iter=np.int32(it))
        return s.log_dem_mesh, y, chisqr, it
