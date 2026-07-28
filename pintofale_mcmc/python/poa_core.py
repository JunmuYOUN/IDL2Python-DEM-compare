"""poa_core.py — PINTofALE support routines converted from IDL.

Sources (PINTofALE current release, PoA_pro_current.tar.gz
sha256 d205c2b896d10d7f9ea6e0b4c2ccd84bff5925bac5ee825c8a01f81b699df166):
  pro/lineflx.pro        -> lineflx
  pro/stat/likeli.pro    -> likeli          (/chi2, /rchi paths only)
  pro/util/varsmooth.pro -> varsmooth       (boxcar path only)
  pro/mk_dem.pro         -> mk_dem_varsmooth('varsmooth' branch only)

Only the branches reachable from mcmc_dem's default configuration are
implemented; every other branch raises NotImplementedError rather than
silently doing something different from IDL.
"""
from __future__ import annotations

import numpy as np

from idl_compat import interpol_n, total

f32 = np.float32

__all__ = ["lineflx", "likeli", "varsmooth", "mk_dem_varsmooth"]

# IDL evaluates `alog(10.)` in FLOAT — the single-precision value of ln(10) is
# what multiplies the temperature step in lineflx, so it must not be float64.
_ALOG10_F32 = np.log(f32(10.0))
# `h=6.626176e-27 & c=2.9979e10` are float literals; h*c*1e8 is a float32 chain.
_HC1E8_F32 = f32(f32(f32(6.626176e-27) * f32(2.9979e10)) * f32(1e8))


def lineflx(line, logT, wvl, Z, dem, abund, noph=False):
    """IDL LINEFLX — fold emissivities through a DEM to get predicted fluxes.

    Parameters follow the IDL positional/keyword split:
      line  (nt, nw) emissivities [1e-23 erg cm^3/s]   (IDL logical order)
      logT  (nt,)    log10 T
      wvl   (nw,)    wavelengths [Ang]
      Z     (nw,)    atomic numbers
      dem   (nt,)    differential emission measure
      abund (30,)    abundances

    Reproduces IDL's float/double mixing exactly: nh_ne is FLTARR (float32),
    abundances are float32, and ln(10) enters in float32.
    """
    line = np.asarray(line)
    nt, nw = line.shape
    tt = np.asarray(logT, dtype=np.float64)
    ang = np.abs(np.asarray(wvl, dtype=np.float64))
    atom = np.asarray(Z).astype(np.int64)

    # `oz=where(atom eq 1)` — note IDL then does nh_ne(oz) where nh_ne is a
    # per-TEMPERATURE array indexed by WAVELENGTH indices. That is an upstream
    # index-space bug; it is inert whenever no line has Z==1, which is the only
    # case reached here. Refuse the buggy path rather than guess at it.
    if np.any(atom == 1):
        raise NotImplementedError(
            "lineflx: Z==1 reaches an upstream nh_ne indexing bug (nh_ne is "
            "sized by temperature but indexed by wavelength); not converted.")

    diffem = np.asarray(dem, dtype=np.float64)
    if diffem.size != nt:
        raise NotImplementedError("lineflx: DEM regridding branch not converted")

    nh_ne = np.full(nt, f32(0.83), dtype=np.float32)   # IDL: fltarr(nt)+0.83

    # mean dlogT, then converted to a natural-log step with FLOAT ln(10)
    dlogt = tt[1:] - tt[:-1]
    mdt = total(dlogt) / np.float64(nt - 1.0)
    mdt = _ALOG10_F32 * mdt            # float32 * float64 -> float64

    avdem = diffem * 1e-23
    ab = np.asarray(abund)
    flx = np.zeros(nw, dtype=np.float64)
    for iw in range(nw):
        zab = ab[atom[iw] - 1]                       # float32 if abund is float32
        cfnc = nh_ne * line[:, iw]                   # float32 * float64 -> float64
        tmp = zab * avdem * cfnc
        flx[iw] = mdt * total(tmp, nan=True)

    if not noph:
        nrg = _HC1E8_F32 / ang                       # float32 / float64 -> float64
        flx = flx / nrg
    return flx                                        # area == 1 in this configuration


def likeli(data, model, dsigma, ulim=None, chi2=False, rchi=False):
    """IDL LIKELI — the default (chi-square) branch.

    Returns chi^2/2 when /CHI2 is set (mcmc_dem always sets it). The upper-limit
    handling is reproduced: with no SOFTLIM, any model above an upper-limit datum
    short-circuits to 2d30.
    """
    x = np.atleast_1d(np.asarray(data))
    m = np.atleast_1d(np.asarray(model))
    nd = x.size
    s = np.atleast_1d(np.asarray(dsigma)).astype(x.dtype, copy=True)
    uu = np.zeros(nd, dtype=np.int64) if ulim is None else np.asarray(ulim).astype(np.int64)
    if nd != m.size:
        return 0.0

    bad = s <= 0
    if np.any(bad):
        s[bad] = 1.0 + np.sqrt(np.abs(x[bad]) + 0.75)

    ou2 = np.nonzero((uu != 0) & (m > x))[0]
    oo = np.nonzero((uu == 0) & (s > 0))[0]
    chi22 = total(((x[oo] - m[oo]) / s[oo]) ** 2) / 2.0 if oo.size else 0.0
    if ou2.size:                       # no SOFTLIM branch
        return 2e30 if chi2 else 1e-30
    if rchi:
        chi22 = chi22 / max(oo.size, 1)
    if not chi2 and not rchi:
        chi22 = min(chi22, 69)
        return min(np.exp(-chi22), 1)
    return chi22


def varsmooth(func, scale):
    """IDL VARSMOOTH — the plain boxcar branch (no TYPE/STEPS/GAUSS/WEIGHT).

    Variable-half-width running mean: for point i with half-width s[i], average
    over [i-s[i], i+s[i]] clipped to the array, normalised by the number of
    points actually used. Weights are all 1 in this configuration.
    """
    f = np.asarray(func, dtype=np.float64).ravel(order="F")
    nf = f.size
    sc = np.asarray(scale)
    s = np.zeros(nf, dtype=f.dtype)
    ns = sc.size
    if ns > 0:
        s[:] = sc.ravel(order="F")[ns - 1]
        if ns <= nf:
            s[:ns] = sc.ravel(order="F")
        else:
            s = sc.ravel(order="F")[:nf].astype(f.dtype)

    sfunc = np.zeros(nf, dtype=f.dtype)
    for i in range(nf):
        scl = int(s[i])                       # IDL truncates when indexing
        im = max(i - scl, 0)
        ip = min(i + scl, nf - 1)
        dim, dip = i - im, ip - i
        acc = f[i]
        if dim > 0:
            acc = acc + total(f[im:i])
        if dip > 0:
            acc = acc + total(f[i + 1:ip + 1])
        sfunc[i] = acc / (dim + dip + 1.0)
    return sfunc


def mk_dem_varsmooth(indem, pardem, nT):
    """IDL MK_DEM('varsmooth', logT=..., indem=..., pardem=...).

    pardem may be None: mcmc_dem_only passes `pardem=lscal`, but LSCAL is not in
    that routine's scope (its own argument is named SCALE), so IDL sees an
    undefined keyword. That makes np==0, which selects the ss=[3.,3.] default —
    i.e. the per-step smoothing scale is a constant 3, and FINDSCALE's result is
    never used inside the MCMC step. This upstream bug is reproduced deliberately.
    """
    dem0 = np.asarray(indem, dtype=np.float64)
    nD = dem0.size
    if pardem is None:
        ss = np.array([3.0, 3.0])
        np_ = 2
    else:
        ss = np.asarray(pardem, dtype=np.float64)
        np_ = ss.size
        if np_ <= 1:
            ss = np.array([3.0, 3.0])
            np_ = 2
    if nD <= 1:
        dem0 = np.array([1.0, 1.0])
        nD = 2
    if np_ != nT:
        ss = interpol_n(ss, nT)
    if nD != nT:
        dem0 = interpol_n(dem0, nT)
    return varsmooth(dem0, ss)
