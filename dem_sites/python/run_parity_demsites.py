"""run_parity_demsites.py — Python twin driver for dem_sites (SITES) parity.

Injects the 5 solver inputs from the IDL oracle 00_input probe (logical orientation),
runs the port with twin probes. Usage: python run_parity_demsites.py <oracle> <py_dir>
"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(oracle_dir, pid):
    sav = readsav(os.path.join(oracle_dir, "probe_%s.sav" % pid), verbose=False)
    rec = np.atleast_1d(sav["chk"])[0]
    return {n.lower(): np.asarray(rec[n]) for n in rec.dtype.names if n.lower() != "probe_id"}


def load_input(oracle_dir):
    d = _load(oracle_dir, "00_input")
    response = np.transpose(d["response"])   # (nwl,nt) -> (nt,nwl)
    # gaussian_function is an IDL/SSW library utility (not the SITES algorithm). Inject the
    # oracle's float32 kernel to test the algorithm at double precision; the port's own
    # gaussian_function reimplementation is separately verified to match (01_setup/ker).
    ker = _load(oracle_dir, "01_setup")["ker"]
    return d["obs_in"], d["err_in"], response, d["response_err"], d["delta_temp"], ker


def main():
    oracle_dir, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True)
    os.environ["CHK_DIR"] = py_dir
    from dem_sites import dem_sites

    obs_in, err_in, response, response_err, delta_temp, ker = load_input(oracle_dir)
    print("injected:", obs_in.shape, response.shape, "ker", ker.shape)
    demmain, demerr, obsmod, irep = dem_sites(obs_in, err_in, response, response_err,
                                              delta_temp, convergence=0.04, ker=ker)
    print("irep", irep, "dem peak", float(np.max(demmain)))


if __name__ == "__main__":
    main()
