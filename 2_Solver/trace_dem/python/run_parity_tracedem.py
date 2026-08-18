"""run_parity_tracedem.py — Python twin driver for trace_dem parity.

Injects data/exptimes (00_input), the T grid (01_basis), and K_inv/tau (02_setup) from
the IDL oracle. Verifies (a) the ported basis functions over T, (b) the linear inversion.
"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(oracle_dir, pid):
    rec = np.atleast_1d(readsav(os.path.join(oracle_dir, "probe_%s.sav" % pid),
                                verbose=False)["chk"])[0]
    return {n.lower(): np.asarray(rec[n]) for n in rec.dtype.names if n.lower() != "probe_id"}


def main():
    oracle_dir, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True)
    os.environ["CHK_DIR"] = py_dir
    from trace_dem import trace_dem, trace_sbasis

    inp = _load(oracle_dir, "00_input")
    basis = _load(oracle_dir, "01_basis")
    setup = _load(oracle_dir, "02_setup")
    data = np.transpose(inp["data"])            # (3,ny,nx)->(nx,ny,3)
    exptimes = inp["exptimes"]
    T = basis["t"]
    K_inv = np.transpose(setup["k_inv"])        # (3,3) logical
    tau = setup["tau"]

    from chk_dump import chk_dump
    # echo the injected setup constants (injection round-trip check)
    chk_dump("02_setup", a=setup["a"], tau=tau,
             k=np.transpose(setup["k"]), k_inv=K_inv)

    # (a) verify the ported basis functions over the injected T grid
    sb0 = trace_sbasis(0, T); sb1 = trace_sbasis(1, T); sb2 = trace_sbasis(2, T)
    chk_dump("01_basis", t=T, sb0=sb0, sb1=sb1, sb2=sb2)

    # (b) verify the linear inversion
    emc, Tavg, emtot = trace_dem(data, exptimes, K_inv, tau)
    print("injected data", data.shape, "K_inv", K_inv.shape, "emc range",
          float(np.min(emc)), float(np.max(emc)))


if __name__ == "__main__":
    main()
