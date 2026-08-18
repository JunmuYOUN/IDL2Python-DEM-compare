"""eit_dem — Python port of the numerical core of eit_dem_tool.pro (Newmark EIT DEM).

Ports the deterministic DEM-computation kernels: e_interp (linear interp), e_fit2
(polynomial fit via matrix inversion), eit_kcorr /iter branch (log-ratio corrections),
and the DEM-map expansion (eit_dem_tool.pro L238-260). Verified against the IDL oracle.

Boundary (NOT ported here): eit_line_map (CHIANTI forward-model that predicts EIT DN
from a DEM) and eit_prep (image I/O) — these are the SSW/CHIANTI response scaffolding;
their outputs (dn_maps, images) are injected, exactly as aia_get_response is for the AIA
solvers. The non-/iter branch of eit_kcorr (which calls eit_line_map) is out of scope.

IDL idioms: element-wise >/< are np.maximum/np.minimum; IDL `invert(mat)##y` = (inv@y.T).T;
temperature-equality tests use float32 (temp is IDL fltarr).
"""
import numpy as np


def e_interp(img1, img2, t1, t2, t3):
    """Simple linear interpolation (eit_dem_tool.pro e_interp)."""
    slope = (img2 - img1) / (t2 - t1)
    return slope * (t3 - t1) + img1


def e_fit2(x, y):
    """Polynomial-coefficient fit via matrix inversion. x has 2/3/5 elements;
    y is logical (npix, nfit). Returns coefs (npix, nfit)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = x.size
    if n == 2:
        mat = np.array([[x[0], 1], [x[1], 1]], dtype=np.float64)
    elif n == 3:
        mat = np.array([[x[0] ** 2, x[0], 1], [x[1] ** 2, x[1], 1],
                        [x[2] ** 2, x[2], 1]], dtype=np.float64)
    elif n == 5:
        mat = np.array([[x[0], 1, 0, 0, 0], [x[1], 1, 0, 0, 0],
                        [0, 0, x[2] ** 2, x[2], 1], [0, 0, x[3] ** 2, x[3], 1],
                        [0, 0, x[4] ** 2, x[4], 1]], dtype=np.float64)
    else:
        raise ValueError("e_fit2 expects 2, 3 or 5 x-values")
    inv = np.linalg.inv(mat)
    return (inv @ y.T).T            # IDL invert(mat) ## y


def eit_kcorr_iter(eitimgs, dn_maps, x):
    """/iter branch of eit_kcorr: log-ratio correction maps from observed images and
    injected predicted-DN maps (dn_maps from eit_line_map). Returns k171,k195,k284,k304."""
    lowsig = [0.1, 0.1, 0.1, 0.1]
    dn_pred_171 = dn_maps["dn171"]
    dn_pred_195 = dn_maps["dn195"]
    cool_284 = dn_maps["cool_284"]
    hot_284 = dn_maps["hot_284"]
    cool_304 = dn_maps["cool_304"]
    warm_304 = dn_maps["warm_304"]
    hot_304 = dn_maps["hot_304"]
    dn_pred_304 = cool_304 + warm_304 + hot_304

    e171, e195, e284, e304 = (eitimgs["eit171"], eitimgs["eit195"],
                              eitimgs["eit284"], eitimgs["eit304"])
    kcorr171 = np.log10(np.maximum(np.maximum(e171, lowsig[0]) / np.maximum(dn_pred_171, lowsig[0]), 1e-5))
    kcorr195 = np.log10(np.maximum(np.maximum(e195, lowsig[1]) / np.maximum(dn_pred_195, lowsig[1]), 1e-5))
    kcorr304 = np.log10(np.maximum(np.maximum(e304, lowsig[3]) / np.maximum(dn_pred_304, lowsig[3]), 1e-5))

    k284_cool = e_interp(kcorr304, kcorr171, x[0], x[1], x[4])
    cool_284 = (10.0 ** k284_cool) * cool_284
    kcorr284 = np.log10(np.maximum((np.maximum(e284, lowsig[2]) - cool_284) / np.maximum(hot_284, 1e-5), 1e-5))

    k304_warm = e_interp(kcorr171, kcorr195, x[1], x[2], x[5])
    warm_304 = (10.0 ** k304_warm) * warm_304
    k304_hot = e_interp(kcorr195, kcorr284, x[2], x[3], x[6])
    hot_304 = (10.0 ** k304_hot) * hot_304
    hot_304 = hot_304 + warm_304
    kcorr304 = np.log10(np.maximum((np.maximum(e304, lowsig[3]) - hot_304) / cool_304, 1e-5))
    return kcorr171, kcorr195, kcorr284, kcorr304


def dem_map_expand(coefs, temp, dem, x):
    """DEM-map expansion (eit_dem_tool.pro L238-260): coefs (nx,ny,5) -> dem_map
    (nx,ny,26) via per-temperature linear/quadratic fit + clamping, then cool extension."""
    coefs = np.asarray(coefs, dtype=np.float64)
    temp = np.asarray(temp)
    dem = np.asarray(dem)
    nx, ny = coefs.shape[0], coefs.shape[1]
    nt = temp.size
    dem_m = np.zeros((nx, ny, nt), dtype=np.float32)
    k63 = None
    x1 = float(x[1])
    for i in range(nt):
        this_t = temp[i]
        if this_t <= x1:
            kcorr_t = coefs[:, :, 0] * this_t + coefs[:, :, 1]
        else:
            kcorr_t = coefs[:, :, 2] * this_t * this_t + coefs[:, :, 3] * this_t + coefs[:, :, 4]
        if this_t == np.float32(6.3):
            k63 = kcorr_t
        if this_t == np.float32(6.4):
            kcorr_t = np.minimum(np.maximum(kcorr_t, k63 - 3), k63 + 0.6)
        if this_t == np.float32(6.5):
            kcorr_t = np.minimum(np.maximum(kcorr_t, k63 - 6), k63 + 1.0)
        dem_m[:, :, i] = np.maximum(dem[i] + np.minimum(kcorr_t, 4), 0)

    dem_scale = dem_m[:, :, 0] / 21.897
    dem2 = np.array([26.391, 25.383, 24.469, 23.654, 22.946, 22.357])
    temp2 = (np.arange(26) * 0.1 + 4.0).astype(np.float32)
    dem_map = np.zeros((nx, ny, temp2.size), dtype=np.float32)
    dem_map[:, :, 6:] = dem_m
    for i in range(6):
        dem_map[:, :, i] = dem2[i] * dem_scale
    return dem_map, temp2
