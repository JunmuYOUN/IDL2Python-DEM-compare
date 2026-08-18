"""xrt_iter_demstat - Python port of the DEM forward-model + chi^2 core of
xrt_dem_iterative2.pro (Weber et al. iterative XRT DEM). Given per-line emissivities,
observed intensities/errors, a DEM curve and a logT grid, returns modeled intensities,
residuals, and total chi^2. i_mod = 10^dem ## (emis*10^t*ln(10^dt)); IDL uses FLTARR
(float32) for emis/p, matched here. The MPFIT spline-inversion loop that repeatedly
evaluates this function is optimizer-path-dependent and is NOT part of this parity core.
"""
import numpy as np
from chk_dump import chk_dump


def xrt_iter_demstat(line_emis, line_t, line_nt, i_obs, i_err, dem, t,
                     weights=None, abunds=None):
    n_line = int(i_obs.shape[0])
    nt = int(t.shape[0])
    if weights is None:
        weights = np.ones(n_line)
    if abunds is None:
        abunds = np.ones(n_line)
    dt = t[1] - t[0]
    emis = np.zeros((n_line, nt), dtype=np.float32)     # IDL FLTARR
    p = np.zeros((n_line, nt), dtype=np.float32)
    for i in range(n_line):
        lnt = int(line_nt[i])
        emis[i, :] = np.maximum(np.interp(t, line_t[i, :lnt], line_emis[i, :lnt]), 0.0)
        p[i, :] = emis[i, :] * 10.0 ** t * np.log(10.0 ** dt)
    ldem = 10.0 ** dem
    i_mod = (p @ ldem) * abunds
    di = i_mod - i_obs
    chisq_arr = di ** 2 * weights / i_err ** 2
    chisq_arr[weights == 0.0] = 0.0
    chisq = np.sum(chisq_arr)
    chk_dump('01_setup', weights=weights, abunds=abunds, emis=emis, p=p)
    chk_dump('99_final', i_mod=i_mod, di=di, chisq=np.float64(chisq))
    return i_mod, di, chisq
