"""mcmc_dem.py — PINTofALE MCMC DEM reconstruction, converted from IDL.

Source: PINTofALE current release, pro/fitting/mcmc_dem.pro + mcmc_dem_only.pro
        (PoA_pro_current.tar.gz, sha256 d205c2b8...df166)
Method: Kashyap & Drake 1998 — Metropolis MCMC over log10(DEM) and log10(abundance).

Scope of this conversion
------------------------
Implemented: the default configuration used for parity — /nosrch or grid search,
/chi2 likelihood, varsmooth DEM smoothing, frozen or thawed abundances.
Not implemented (each raises): LOOPY, SPLINY, ONLYRAT/ratios, MIXIE deblending,
SYSDEV emissivity jitter, SOFTLIM, SAMPENV, upper limits, savfil.

Injected rather than converted (established project pattern — support-library
outputs come from the oracle, the algorithm itself is converted):
  * `abnd`  <- GETABUND('anders & grevesse')
  * `lscal` <- FINDSCALE(line)  (which calls WVLT_SCALE)

Faithful reproduction of upstream defects
-----------------------------------------
1. `xfrac = 0*wvl*1.0` in mcmc_dem (mcmc_dem_only uses `0*wvl+1`), so the
   pre-loop flux is divided by zero: ff = +Inf and oldprb = bestprb = +Inf.
   Consequence: accept_fn = min(Inf - prb, 0) = 0 at the first step, so the
   first proposal is always accepted. Reproduced exactly.
2. mcmc_dem_only passes `pardem=lscal`, an undefined variable in its scope, so
   the per-step smoothing half-width is a constant 3 (see poa_core).
3. PRED_FLX is called once per outer batch purely to feed a progress printout;
   it has no effect on the chain and is omitted here.
"""
from __future__ import annotations

import numpy as np

from idl_compat import (alog10_f32, alog_f32, exp_f32, histogram, idl_sort,
                        interpol_xu, pow10, total)
from idl_random import IDLRandom
from poa_core import likeli, lineflx, mk_dem_varsmooth

f32 = np.float32

__all__ = ["mcmc_dem", "mcmc_dem_only", "McmcResult"]


class McmcResult:
    """Outputs of mcmc_dem (the IDL return value plus its output keywords)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def mcmc_dem_only(pars, fobs, sig, ulim, line, logT, wvl, Z, nt, rchi=False,
                  xfrac=None):
    """IDL MCMC_DEM_ONLY — decode parameters, predict fluxes, return chi^2/2.

    Returns (ff, prb, dem, abnd).

    `xfrac` is IDL's *by-reference* keyword: mcmc_dem_only unconditionally
    assigns `xfrac=0*wvl+1` (no MIXSSTR), which overwrites the caller's copy.
    mcmc_dem relies on that side effect — see the note in mcmc_dem below — so the
    caller's array is mutated in place here rather than merely read.
    """
    dem = pars[:nt]
    abnd = pars[nt:]
    # pardem=lscal is out of scope in IDL -> constant half-width 3 (see poa_core)
    dem = mk_dem_varsmooth(dem, None, nt)
    dem = pow10(dem)
    abnd = pow10(abnd)
    # adjustie() is a no-op without TIES
    local_xfrac = np.zeros(len(wvl)) + 1.0      # IDL: 0*wvl+1
    if xfrac is not None:
        xfrac[:] = local_xfrac                  # by-reference write-back
    ff = lineflx(line, logT, wvl, Z, dem, abnd) / local_xfrac
    prb = likeli(fobs, ff, sig, ulim=ulim, chi2=True, rchi=rchi)
    return ff, prb, dem, abnd


def mcmc_dem(wvl, flx, emis, Z, logt, diffem, abund, fsigma, lscal,
             nsim=100, nbatch=10, nburn=None, nosrch=True, bound=None,
             abrng=None, seed=42, trace=None):
    """IDL MCMC_DEM. `emis` is (nt, nw) in IDL logical order.

    `abund` and `lscal` are the injected GETABUND / FINDSCALE outputs.
    `trace` — optional dict that collects per-step diagnostics for parity probes.
    """
    rng = IDLRandom(seed)

    # ---- cast variables ----------------------------------------------------
    ww = np.asarray(wvl, dtype=np.float64)
    fx = np.asarray(flx, dtype=np.float64)
    line = np.asarray(emis, dtype=np.float64)
    nw, nf = ww.size, fx.size
    zz = np.asarray(Z).astype(np.int64)
    sig = np.asarray(fsigma, dtype=np.float64).copy()
    ulim = np.zeros(nw, dtype=np.int64)
    logt = np.asarray(logt, dtype=np.float64)
    nt = logt.size

    dem = np.asarray(diffem, dtype=np.float64).copy()
    abnd = np.asarray(abund, dtype=np.float32).copy()      # getabund -> float32
    abmn = 1e-20
    abnd[abnd <= 0] = abmn
    if abnd[0] != 1:
        raise ValueError("H abundance must be 1 (upstream HALTs otherwise)")

    nsim = max(int(nsim), 10)
    nbatch = max(int(nbatch), 5)
    nburn = int(nsim / 10.0 + 1) if nburn is None else int(nburn)

    # ---- parameter ranges ---------------------------------------------------
    rngD = np.concatenate([dem / 1e5, dem * 1e5]).astype(np.float64)
    rngD = rngD.reshape(2, nt).T                            # IDL reform(nt,2)
    nab = abnd.size
    # IDL: rngA starts as the FLOAT vector [abnd,abnd]; when ABRNG is supplied its
    # (double) values are assigned *into* that float array, so they are rounded to
    # float32 before REFORM. Reproduce the narrowing.
    rngA_flat = np.concatenate([abnd, abnd])
    if abrng is not None:
        tmp = np.asarray(abrng, dtype=np.float64).reshape(nab, 2).ravel(order="F")
        rngA_flat[:] = tmp[:2 * nab].astype(np.float32)
    rngA = rngA_flat.reshape(2, nab).T

    rngD[:, 0] = np.where(rngD[:, 0] <= 0, 1.0, rngD[:, 0])
    rngD[:, 1] = np.where(rngD[:, 1] <= 0, 1.0, rngD[:, 1])
    swap = rngD[:, 0] > rngD[:, 1]
    rngD[swap, 0], rngD[swap, 1] = rngD[swap, 1].copy(), rngD[swap, 0].copy()
    rngA[:, 0] = np.where(rngA[:, 0] <= 0, abmn, rngA[:, 0])
    rngA[:, 1] = np.where(rngA[:, 1] <= 0, abmn, rngA[:, 1])
    swap = rngA[:, 0] > rngA[:, 1]
    rngA[swap, 0], rngA[swap, 1] = rngA[swap, 1].copy(), rngA[swap, 0].copy()

    # ---- thawed parameters --------------------------------------------------
    # IDL's ALOG10 of a FLOAT array returns FLOAT: the abundance half of these
    # vectors is float32-rounded before being promoted by the concatenation.
    allpar = np.concatenate([np.log10(dem), alog10_f32(abnd).astype(np.float64)])
    ipar = np.arange(allpar.size, dtype=np.int64)
    aparmx = np.concatenate([np.log10(rngD[:, 1]),
                             alog10_f32(rngA[:, 1]).astype(np.float64)])
    aparmn = np.concatenate([np.log10(rngD[:, 0]),
                             alog10_f32(rngA[:, 0]).astype(np.float64)])
    ipar[aparmx == aparmn] = -1
    opar = np.nonzero(ipar >= 0)[0]
    npar = opar.size
    if npar == 0:
        raise ValueError("no thawed parameters")

    # ---- histogram bins -----------------------------------------------------
    wdem = wab = 0.1
    mindem, maxdem = np.nanmin(rngD), np.nanmax(rngD)
    minab, maxab = np.nanmin(rngA), np.nanmax(rngA)
    mindem, maxdem = np.log10(mindem), np.log10(maxdem)
    minab, maxab = alog10_f32(minab), alog10_f32(maxab)
    n_dem = int((maxdem - mindem) / wdem + 0.5) + 1
    # IDL: FINDGEN(n)*wdem is float32, but adding the double MINDEM promotes the
    # sum to double. numpy would keep float32 when adding a Python scalar, so the
    # promotion is made explicit.
    xdem = (np.arange(n_dem, dtype=np.float32) * f32(wdem)).astype(np.float64) + mindem
    n_ab = int((maxab - minab) / wab + 0.5) + 1
    xab = (np.arange(n_ab, dtype=np.float32) * f32(wab) + f32(minab)).astype(np.float32)
    if xdem.max() < maxdem:
        xdem = np.append(xdem, xdem.max() + wdem)
        n_dem += 1
    if xab.max() < maxab:
        xab = np.append(xab, f32(xab.max() + wab))
        n_ab += 1

    # ---- emissivity envelope ------------------------------------------------
    emenv = np.array([total(line[i, :]) for i in range(nt)], dtype=np.float64)
    frozen = np.nonzero(aparmx[:nt] == aparmn[:nt])[0]
    if frozen.size:
        emenv[frozen] = 0.0
    oo = np.nonzero(aparmx[:nt] != aparmn[:nt])[0]
    moo = oo.size
    cenv = np.zeros(moo, dtype=np.float64) + emenv[0]
    for i in range(1, moo):
        cenv[i] = cenv[i - 1] + emenv[oo[i]]
    cenv = np.maximum(cenv / cenv[moo - 1], 0)

    lscal = np.asarray(lscal)
    # IDL: `xfrac = 0*wvl*1.0` -> zeros. MCMC_DEM_ONLY takes XFRAC as a
    # by-reference keyword and overwrites it with ones, so whether the pre-loop
    # flux below is divided by zero depends on whether the grid search ran.
    xfrac = np.zeros(nw) * 1.0

    # ---- initialise ---------------------------------------------------------
    isim = 0
    ksim = int(nsim) * int(npar)
    istor = 0
    iburn = 0
    # IDL: `ebound=0.9` is a FLOAT literal when BOUND is not supplied, and only
    # becomes double via `double(bound(0))` when it is. EBOUND*MOU is therefore
    # float32 in the default case, which shifts the confidence bounds by ~3e-8.
    ebound = f32(0.9) if not bound else np.float64(min(max(float(bound), 0.1), 1.0))
    sigpar = (aparmx - aparmn) / 5.0

    parcnt = np.zeros(npar, dtype=np.int64)
    storpar = np.zeros(ksim, dtype=np.float64)
    stordem = np.zeros(ksim, dtype=np.float64)
    storidx = np.zeros(ksim, dtype=np.int64) - 1
    bestpar, bestsig = allpar.copy(), sigpar.copy()

    simprb = np.zeros(nsim + 1, dtype=np.float32)
    simdem = np.zeros((nt, nsim + 1), dtype=np.float64)
    simabn = np.zeros((nab, nsim + 1), dtype=np.float64)
    simflx = np.zeros((nf, nsim + 1), dtype=np.float32)

    # the oracle's 01_setup probe is dumped BEFORE the grid search, so snapshot here
    if trace is not None:
        trace["setup"] = dict(allpar=allpar.copy(), ipar=ipar.copy(),
                              opar=opar.copy(), aparmx=aparmx.copy(),
                              aparmn=aparmn.copy(), sigpar=sigpar.copy(),
                              lscal=lscal.copy(), emenv=emenv.copy(),
                              cenv=cenv.copy(), xdem=xdem.copy(), xab=xab.copy(),
                              rngd=rngD.copy(), rnga=rngA.copy())

    if not nosrch:
        # IDL: crude 1-D grid search per thawed parameter to seed the chain.
        # PRBPAR is FLTARR, so the softmax runs in single precision (expf).
        for i in range(npar):
            j = opar[i]
            xx = xdem if j < nt else xab
            nx = xx.size
            prbpar = np.zeros(nx, dtype=np.float32)
            for k in range(nx):
                tmpar = allpar.copy()
                tmpar[j] = xx[k]
                _, prb_k, _, _ = mcmc_dem_only(tmpar, fx, sig, ulim, line, logt,
                                               ww, zz, nt, rchi=True, xfrac=xfrac)
                prbpar[k] = prb_k
            prbpar = -prbpar
            prbpar = np.maximum(np.minimum(prbpar - prbpar.max(), f32(69)), f32(-69))
            prbpar = exp_f32(prbpar)
            prbpar = prbpar / total(prbpar)
            oo = np.nonzero(prbpar == prbpar.max())[0]
            if oo.size > 1:
                bestpar[j] = total(prbpar[oo] * xx[oo]) / total(prbpar[oo])
            elif oo.size == 1:
                bestpar[j] = xx[oo][0]
            minsig = 3 * abs(np.min(xx[1:] - xx[:-1]))
            bestsig[j] = max(total(prbpar * xx * xx) - total(prbpar * xx) ** 2,
                             minsig ** 2)
            bestsig[j] = np.sqrt(bestsig[j])
            allpar[j] = bestpar[j]

    allpar, sigpar = bestpar.copy(), bestsig.copy()
    bestpar, bestsig = allpar.copy(), sigpar.copy()

    # ---- pre-loop state: upstream divides by a zero XFRAC (see module docstring)
    dem = pow10(allpar[:nt])
    abnd = pow10(allpar[nt:])
    with np.errstate(divide="ignore", invalid="ignore"):
        ff = lineflx(line, logt, ww, zz, dem, abnd) / xfrac
    oldprb = likeli(fx, ff, sig, ulim=ulim, chi2=True)
    bestprb = oldprb
    bestdem, bestabnd, bestflx = dem, abnd, ff

    if trace is not None:
        trace.setdefault("batch", [])
        trace.setdefault("step", [])
        trace["init"] = dict(dem=dem.copy(), abnd=abnd.copy(), ff=ff.copy(),
                             oldprb=oldprb, bestprb=bestprb,
                             allpar=allpar.copy(), sigpar=sigpar.copy(),
                             bestpar=bestpar.copy())

    # ---- MCMC ---------------------------------------------------------------
    prddat = None
    while isim < nsim:
        for ibtc in range(nbatch):
            rn = rng.randomn(npar)                       # float32 stream
            ru = rng.randomu(2 * npar + 1)
            hitpar = np.zeros(npar, dtype=np.int64)
            nnpar = npar + 1
            if trace is not None and len(trace["batch"]) < 20:
                trace["batch"].append(dict(rn=rn.copy(), ru=ru.copy(),
                                           allpar=allpar.copy(),
                                           sigpar=sigpar.copy(),
                                           isim=isim, ibtc=ibtc))

            for jpar in range(nnpar):
                tmpar = allpar.copy()
                if jpar < npar and opar[jpar] < nt:
                    kpar = int(np.nonzero(cenv >= ru[jpar])[0][0])
                    hitpar[kpar] = 1
                else:
                    kpar = jpar

                if kpar < npar:
                    tmpar[opar[kpar]] = tmpar[opar[kpar]] + rn[jpar] * sigpar[opar[kpar]]
                else:
                    ono = np.nonzero(hitpar[opar] == 0)[0]
                    if ono.size:
                        tmpar[opar[ono]] = (tmpar[opar[ono]]
                                            + rng.randomn(ono.size) * sigpar[opar[ono]])

                # force compliance with hard limits
                oo = np.nonzero(((tmpar > aparmx) | (tmpar < aparmn)) & (ipar >= 0))[0]
                while oo.size:
                    tmpar[oo] = allpar[oo] + rng.randomn(oo.size) * sigpar[oo]
                    oo = np.nonzero(((tmpar > aparmx) | (tmpar < aparmn)) & (ipar >= 0))[0]

                tmpff, prb, tmpdem, tmpabnd = mcmc_dem_only(
                    tmpar, fx, sig, ulim, line, logt, ww, zz, nt)

                accept_fn = min(oldprb - prb, 0.0)
                # ru is FLOAT: IDL's ALOG here is single precision (logf)
                metro_check = alog_f32(ru[jpar + npar])
                accepted = metro_check < accept_fn
                if accepted:
                    allpar, dem, oldprb, abnd = tmpar, tmpdem, prb, tmpabnd
                    ff, prddat = tmpff, None

                if trace is not None and len(trace["step"]) < 20:
                    trace["step"].append(dict(
                        tmpar=tmpar.copy(), tmpdem=tmpdem.copy(),
                        tmpabnd=tmpabnd.copy(), tmpff=tmpff.copy(), prb=prb,
                        accept_fn=accept_fn, metro_check=metro_check,
                        allpar=allpar.copy()))

                if prb < bestprb:
                    bestprb, bestpar, bestdem, bestabnd = prb, allpar, dem, abnd
                    bestflx = ff
                    simprb[nsim] = bestprb
                    simdem[:, nsim] = bestdem
                    simabn[:, nsim] = bestabnd
                    simflx[:, nsim] = bestflx

                if kpar < npar:
                    storpar[istor] = tmpar[opar[kpar]]
                    storidx[istor] = opar[kpar]
                    if opar[kpar] < nt:
                        ko = opar[kpar]
                        stordem[istor] = min(max(np.log10(abs(dem[ko])), aparmn[ko]),
                                             aparmx[ko])
                    if iburn == nburn:
                        parcnt[kpar] += 1
                        if parcnt[kpar] == int(parcnt[kpar] / nbatch) * nbatch:
                            istor += 1
                    else:
                        istor += 1

        simprb[isim] = oldprb
        simdem[:, isim] = dem
        simabn[:, isim] = abnd
        simflx[:, isim] = ff
        dem = np.log10(simdem[:, isim])
        abnd = pow10(allpar[nt:])
        dem = pow10(dem)

        if iburn < nburn:
            jstor = istor - nbatch * npar
            norm = 0
            tmpsig = sigpar.copy()
            dsig = 0.0
            if jstor >= 1:
                opars, pars_ = storpar[:jstor], storpar[:istor]
                opari, pari = storidx[:jstor], storidx[:istor]
                for i in range(npar):
                    j = opar[i]
                    oo0 = np.nonzero(opari == j)[0]
                    oo1 = np.nonzero(pari == j)[0]
                    if oo1.size < oo0.size:
                        raise RuntimeError("insect infestation!")
                    if oo0.size == 0:
                        continue
                    xx = xdem if j < nt else xab
                    minx, maxx = xx.min(), xx.max()
                    xx = 0.5 * (xx[1:] + xx[:-1])
                    f0, f1 = opars[oo0], pars_[oo1]
                    binx = abs(xx[1] - xx[0])
                    hf0 = histogram(f0, minx, maxx, binx).astype(np.float32)
                    hf1 = histogram(f1, minx, maxx, binx).astype(np.float32)
                    hf0 = hf0 / total(hf0)
                    hf1 = hf1 / total(hf1)
                    dsig = dsig + total(hf0 * hf1) ** 2 / total(hf0 ** 2) / total(hf1 ** 2)
                    norm += 1
                    os_ = idl_sort(f1)
                    f1 = f1[os_]
                    moo1 = oo1.size
                    xf1 = np.arange(moo1, dtype=np.float32) / np.float32(moo1)
                    if moo1 > 3:
                        # IDL writes [0.16,0.84] as FLOAT literals
                        t = interpol_xu(f1, xf1, np.float32([0.16, 0.84]))
                        if t[0] < t[1]:
                            tmpsig[j] = max(3.0 * (t[1] - t[0]), 0.3 * binx)
                        else:
                            tmpsig[j] = 3 * binx
                    else:
                        tmpsig[j] = max(f1.max() - f1.min(), 3.0 * binx)
            dsig = dsig / max(norm, 1)
            iburn = nburn if dsig > 0.98 else iburn - 1
            if abs(iburn) == nburn:
                iburn = nburn
            if istor >= ksim:
                iburn = nburn
            if iburn == nburn:
                istor = 0
                storpar[:] = 0.0
                storidx[:] = -1
                sigpar = tmpsig
        else:
            isim += 1

    # ---- outputs ------------------------------------------------------------
    def _bounds(vals_row, centre, n, dtype):
        """IDL's confidence-bound block. DEMERR is DBLARR but ABERR is FLTARR."""
        err = np.zeros((n, 2), dtype=dtype)
        for i in range(n):
            err[i, :] = centre[i]
            tmp = np.sort(vals_row[i, :])
            if total(tmp) > 0:
                tu = tmp[tmp >= centre[i]]
                tl = tmp[tmp <= centre[i]]
                if tl.size > 1:
                    # IDL: interpol(tmpl, reverse(findgen(mol)), ebound*mol)
                    x = np.arange(tl.size, dtype=np.float32)[::-1].astype(np.float64)
                    # EBOUND*MOL stays in EBOUND's precision (float32 by default), and
                    # rounding that product is what picks the interpolation point.
                    u = ebound * np.asarray(tl.size, dtype=ebound.dtype)
                    err[i, 0] = interpol_xu(tl, x, [u])[0]
                if tu.size > 1:
                    x = np.arange(tu.size, dtype=np.float32).astype(np.float64)
                    u = ebound * np.asarray(tu.size, dtype=ebound.dtype)
                    err[i, 1] = interpol_xu(tu, x, [u])[0]
        return err

    dem_out = simdem[:, nsim].copy()
    demerr = _bounds(simdem[:, :nsim + 1], dem_out, nt, np.float64)
    abund_out = simabn[:, nsim].copy()
    aberr = _bounds(simabn[:, :nsim + 1], abund_out, nab, np.float32)

    return McmcResult(dem=dem_out, demerr=demerr, abund=abund_out, aberr=aberr,
                      simprb=simprb, simdem=simdem, simabn=simabn, simflx=simflx,
                      storpar=storpar, storidx=storidx, stordem=stordem,
                      bestpar=bestpar, allpar=allpar, sigpar=sigpar)
