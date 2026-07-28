"""run_parity_demreg.py — Python twin driver for demreg re-verification.

Injects the wrapper inputs (dn,edn,tresp,tresp_logt,temps) from the IDL oracle's
00_input probe, then runs the existing python port (dn2dem_pos) with twin probes.
dem_norm0 defaults to ones inside the port (matches the oracle call).

Usage: python run_parity_demreg.py <oracle_probe_dir> <py_probe_dir>
"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # dn2dem_pos.py, demmap_pos.py, chk_dump.py live here


def load_oracle_input(oracle_dir):
    sav = readsav(os.path.join(oracle_dir, "probe_00_input.sav"), verbose=False)
    rec = np.atleast_1d(sav["chk"])[0]
    d = {n.lower(): np.asarray(rec[n]) for n in rec.dtype.names if n.lower() != "probe_id"}
    # logical orientation: reverse all axes of >=2D arrays (single pixel dn_in is 1D)
    dn_in = np.transpose(d["dn_in"]) if d["dn_in"].ndim >= 2 else d["dn_in"]
    edn_in = np.transpose(d["edn_in"]) if d["edn_in"].ndim >= 2 else d["edn_in"]
    tresp = np.transpose(d["tresp"])      # (nf,ngd) -> (ngd,nf)  [orientation: logical]
    tresp_logt = d["tresp_logt"]          # (ngd,)
    temps = d["temps"]                    # (nedge,)
    return dn_in, edn_in, tresp, tresp_logt, temps


def main():
    oracle_dir, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True)
    os.environ["CHK_DIR"] = py_dir
    from dn2dem_pos import dn2dem_pos

    dn_in, edn_in, tresp, tresp_logt, temps = load_oracle_input(oracle_dir)
    print("injected:", dn_in.shape, edn_in.shape, tresp.shape, tresp_logt.shape, temps.shape)
    dem, edem, elogt, chisq, dn_reg = dn2dem_pos(dn_in, edn_in, tresp, tresp_logt, temps)
    print("dem", np.asarray(dem).shape, "chisq", float(np.asarray(chisq)),
          "dem range", float(np.min(dem)), float(np.max(dem)))


if __name__ == "__main__":
    main()
