"""vdem — Python port of vdem.pro (Velocity DEM; Newton, Emslie & Mariska 1995).

NOTE: this is a VELOCITY DEM (line-of-sight velocity distribution of an emitting ion,
units photons/s/(km/s)) — a DIFFERENT quantity from the temperature DEMs. GCV-regularized
linear inversion of a spectral line profile.

Faithful 1:1 port of the non-CAXIX, non-VERBOSE core. IDL idioms preserved:
  - IDL `#` (2D) -> numpy `@` under logical storage: A=kernel@kernel.T, baseh=ke3@ke3.T,
    data=kernel@inten, retflux=invertvdem@basekernel, covariance=inv@A@inv.
  - v5_invert (LUDC+LUSOL) -> np.linalg.solve(matrix0, data); invert -> np.linalg.inv.
  - `fix` -> trunc-toward-zero; element-wise `>` -> np.maximum.
  - WEIGHT mode (q = weight*initial) is deterministic; the GCV/Brent optimizer
    (gcvcompute) is verified separately (gcvfunction) and is a path-dependent boundary.
chk_dump twin-probes mirror staging/vdem_probed.pro.
"""
import numpy as np

from chk_dump import chk_dump


def _f32(x):
    """IDL `e` literals are FLOAT32; reproduce their exact value in float64 arithmetic."""
    return float(np.float32(x))


K_BOLTZ = _f32(1.38054e-16)          # erg/K   (IDL float32 literal)
C_LIGHT = _f32(2.997925e10)          # cm/s
# sun_sat_dist and its square are pure float32 ops in IDL -> do them in numpy float32
_AU32 = np.float32(1.495979e13)
_SSD32 = np.float32(0.99) * _AU32                 # float32*float32 -> float32
AU = float(_AU32)
SUN_SAT_DIST = float(_SSD32)
_SSD2 = float(_SSD32 * _SSD32)                     # float32^2 -> float32
# IDL !PI is FLOAT32 (3.1415927), not double! (!DPI is the double one)
_PI32 = float(np.float32(np.pi))
_SQRT_2PI = float(np.sqrt(np.float32(2.0) * np.float32(np.pi)))   # IDL sqrt(2.*!pi) in float32
# flux scaling accumulates 4*!pi*ssd^2 in DOUBLE (into the double flux);
# but evdem/retflux compute 4.*!pi*sun_sat_dist^2 STANDALONE -> all-float32 product:
_SCALE_F32 = float((np.float32(4.0) * np.float32(np.pi)) * np.float32(_SSD2))


def vdem(flux_in, eflux_in, wave, mass, L0, tbar, weight=1.0, xi=0.0, width_sumer=0.0):
    """flux_in, eflux_in: (nflux, nspec). wave: (nflux,). Returns (vdem, evdem, velocity, retflux)."""
    flux_in = np.atleast_2d(np.asarray(flux_in, dtype=np.float64))
    eflux_in = np.atleast_2d(np.asarray(eflux_in, dtype=np.float64))
    wave = np.asarray(wave, dtype=np.float64)
    mass = float(mass); L0 = float(L0); tbar = float(tbar)
    if flux_in.shape[0] != wave.size:      # ensure (nflux, nspec)
        flux_in = flux_in.T; eflux_in = eflux_in.T
    chk_dump("00_input", fin=np.squeeze(flux_in), efin=np.squeeze(eflux_in), wavein=wave,
             mass=np.float64(mass), l0=np.float64(L0), tbar=np.float64(tbar))

    m = mass * _f32(1.6713e-24)
    scale = 4.0 * _PI32 * _SSD2                      # 4.*!pi*sun_sat_dist^2 (!pi & ssd^2 float32)
    flux = flux_in * scale                 # (nflux, nspec)
    eflux = eflux_in
    nflux = wave.size
    specnum = flux.shape[1]

    # velocity bins
    diff = wave[1] - wave[0]
    dvel = C_LIGHT * diff / L0
    vmax = -1.0 * C_LIGHT * (np.min(wave) - L0) / L0
    vmin = -1.0 * C_LIGHT * (np.max(wave) - L0) / L0
    nvel = int(np.trunc((vmax - vmin) / dvel))
    velocity = np.zeros(nvel)
    velocity[0] = vmin + dvel / 2.0
    for i in range(1, nvel):
        velocity[i] = velocity[i - 1] + dvel

    # base kernel (Gaussian), non-CAXIX branch
    v_instr = width_sumer * C_LIGHT / L0
    center = L0 * (1.0 - velocity / C_LIGHT)
    sig = L0 * np.sqrt(K_BOLTZ * tbar / m + xi ** 2 + v_instr ** 2) / C_LIGHT
    basekernel = np.zeros((nvel, nflux))
    for j in range(nflux):
        for i in range(nvel):
            val = dvel / (_SQRT_2PI * sig) * \
                np.exp(-(wave[j] - center[i]) ** 2 / (2.0 * sig ** 2))
            basekernel[i, j] = 0.0 if val < 1.0e-15 else val
    chk_dump("01_kernel", basekernel=basekernel, velocity=velocity,
             sig=np.float64(sig), dvel=np.float64(dvel))

    # third-difference regularization matrix
    ke3 = np.zeros((nvel, nvel - 3))
    for n in range(nvel - 3):
        ke3[n, n] = -1.0; ke3[n + 1, n] = 3.0
        ke3[n + 2, n] = -3.0; ke3[n + 3, n] = 1.0
    baseh = ke3 @ ke3.T
    h = baseh * _f32(1.0e16) * _f32(1.0 / 6.0) / _f32(1.0e4)   # IDL: all float32 literals
    TrH = np.trace(h)
    chk_dump("02_h", h=h, trh=np.float64(TrH))

    invertvdem = np.zeros((nvel, specnum))
    evdem = np.zeros((nvel, specnum))
    retflux = np.zeros((nflux, specnum))
    for time in range(specnum):
        kernel = basekernel / eflux[:, time][None, :]          # (nvel,nflux)/eflux[f]
        A = kernel @ kernel.T                                   # (nvel,nvel)
        inten = flux[:, time] / eflux[:, time]                  # (nflux,)
        TrA = np.trace(A)
        initial = TrA / TrH
        data = kernel @ inten                                   # (nvel,)  = K^T g
        q = weight * initial                                    # WEIGHT mode
        matrix0 = A + q * h
        x = np.linalg.solve(matrix0, data)
        invertvdem[:, time] = x
        inverse = np.linalg.inv(matrix0)
        covariance = inverse @ A @ inverse
        evdem[:, time] = _SCALE_F32 * np.sqrt(np.diag(covariance))  # float32 4pi*ssd^2
        retflux[:, time] = (x @ basekernel) / _SCALE_F32           # float32 4pi*ssd^2
        if time == 0:
            chk_dump("03_setup", kernel=kernel, a=A, inten=inten,
                     data=data[None, :], initial=np.float64(initial),
                     q=np.float64(q), matrix0=matrix0)

    vd = np.squeeze(invertvdem)
    ev = np.squeeze(evdem)
    rf = np.squeeze(retflux)
    chk_dump("99_final", invertvdem=vd, evdem=ev, velocity=velocity, retflux=rf)
    return vd, ev, velocity, rf
