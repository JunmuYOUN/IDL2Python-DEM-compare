"""trace_dem — Python port of trace_dem.pro (TRACE 3-channel linear DEM, Kankelborg 1998).

Ports the basis functions (sigmoid, trace_sbasis) and the core linear inversion
(emc = K_inv # p per pixel; Tavg; emtot). The setup constants K_inv/tau are computed by
trace_dem_setup via trace_t_resp (SSW response) + int_tabulated (IDL 5-point Newton-Cotes
on a spline-resampled grid); those are injected from the oracle — int_tabulated's exact
replication is the documented conversion boundary. Verified against the IDL oracle.
"""
import numpy as np

from chk_dump import chk_dump

PI = 3.14159265359


def sigmoid(theta):
    """IDL sigmoid: sin(theta)^2, zeroed where theta>pi or theta<0."""
    theta = np.asarray(theta, dtype=np.float64)
    result = np.sin(theta) ** 2
    result[(theta > PI) | (theta < 0.0)] = 0.0
    return result


def trace_sbasis(index, T):
    """IDL trace_sbasis: 3 sinusoidal bell basis functions (index 0,1,2)."""
    logT = np.log10(np.asarray(T, dtype=np.float64))
    logTa, logTb, logTc = 5.95, 6.31, 5.0
    logT0 = {0: logTa, 1: (logTa + logTb) / 2.0, 2: logTb}[index]
    omega = 0.5 * PI / (logTb - logTa)
    theta = omega * (logT - logT0) + 0.5 * PI
    result = sigmoid(theta)
    if index == 0:
        result[logT < logTa] = 1.0
        result[logT < logTc] = 0.0
    return result


def trace_dem(data, exptimes, K_inv, tau):
    """Return (emc, Tavg, emtot). data (nx,ny,3) DN; K_inv (3,3); tau (3,)."""
    data = np.asarray(data, dtype=np.float64)
    exptimes = np.asarray(exptimes, dtype=np.float64)
    K_inv = np.asarray(K_inv, dtype=np.float64)
    tau = np.asarray(tau, dtype=np.float64)
    chk_dump("00_input", data=data, exptimes=exptimes)

    nx, ny = data.shape[0], data.shape[1]
    p = data / exptimes[None, None, :]              # DN/s
    emc = np.zeros((nx, ny, 3))
    Tavg = np.zeros((nx, ny))
    for x in range(nx):
        for y in range(ny):
            emc[x, y, :] = K_inv @ p[x, y, :]        # IDL K_inv # p
            Tavg[x, y] = np.sum(emc[x, y, :] * tau) / np.sum(emc[x, y, :])
    emtot = np.sum(emc, axis=2)                      # IDL total(emc,3)
    chk_dump("99_final", emc=emc, Tavg=Tavg, emtot=emtot)
    return emc, Tavg, emtot
