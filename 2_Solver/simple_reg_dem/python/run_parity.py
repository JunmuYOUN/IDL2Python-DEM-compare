"""run_parity.py — Python twin driver for simple_reg_dem parity.

Injects the solver inputs from the IDL oracle's 00_input probe (parity-protocol
input-injection pattern), so the data-generation scaffolding is NOT part of the
compared surface — only simple_reg_dem itself is. readsav reverses IDL dims; under
orientation=logical we transpose back to IDL logical layout before feeding the port.

Usage:
  python run_parity.py <oracle_probe_dir> <py_probe_dir> <i0> <j0>
"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)  # simple_reg_dem.py, chk_dump.py live here (converted/)


def load_oracle_input(oracle_dir):
    sav = readsav(os.path.join(oracle_dir, "probe_00_input.sav"), verbose=False)
    rec = np.atleast_1d(sav["chk"])[0]
    d = {n.lower(): rec[n] for n in rec.dtype.names if n.lower() != "probe_id"}
    # logical orientation: reverse all axes to restore IDL index order
    data = np.transpose(np.asarray(d["data"]))      # (nchan,ny,nx) -> (nx,ny,nchan)
    errors = np.transpose(np.asarray(d["errors"]))
    tresps = np.transpose(np.asarray(d["tresps"]))  # (nchan,nt) -> (nt,nchan)
    exptimes = np.asarray(d["exptimes"])            # 1-D, unchanged
    logt = np.asarray(d["logt"])                    # 1-D, unchanged
    return data, errors, exptimes, logt, tresps


def main():
    oracle_dir, py_dir, i0, j0 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    os.makedirs(py_dir, exist_ok=True)
    os.environ["CHK_DIR"] = py_dir
    os.environ["PROBE_I0"] = str(i0)
    os.environ["PROBE_J0"] = str(j0)

    from simple_reg_dem import simple_reg_dem  # after CHK_DIR/PROBE_* set

    data, errors, exptimes, logt, tresps = load_oracle_input(oracle_dir)
    print("injected shapes:", data.shape, errors.shape, exptimes.shape,
          logt.shape, tresps.shape, "dtypes:", data.dtype, tresps.dtype)
    dems, chi2 = simple_reg_dem(data, errors, exptimes, logt, tresps)
    print("done dems", dems.shape, dems.dtype, "chi2", chi2.shape,
          "chi2 range", float(np.min(chi2)), float(np.max(chi2)))


if __name__ == "__main__":
    main()
