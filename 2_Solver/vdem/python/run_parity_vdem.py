"""run_parity_vdem.py — Python twin driver for vdem (velocity DEM) parity.
Injects 00_input, runs the port. Usage: python run_parity_vdem.py <oracle> <py_dir>"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load(oracle_dir):
    rec = np.atleast_1d(readsav(oracle_dir + "/probe_00_input.sav", verbose=False)["chk"])[0]
    d = {n.lower(): np.asarray(rec[n]) for n in rec.dtype.names if n.lower() != "probe_id"}
    return d["fin"], d["efin"], d["wavein"], float(d["mass"]), float(d["l0"]), float(d["tbar"])


def main():
    oracle, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True)
    os.environ["CHK_DIR"] = py_dir
    from vdem import vdem

    fin, efin, wavein, mass, l0, tbar = load(oracle)
    fin = np.atleast_2d(fin); efin = np.atleast_2d(efin)
    if fin.shape[0] != wavein.size:
        fin = fin.T; efin = efin.T
    vd, ev, velo, rf = vdem(fin, efin, wavein, mass, l0, tbar, weight=1.0)
    print("vdem", np.asarray(vd).shape, "min/max", float(np.min(vd)), float(np.max(vd)))


if __name__ == "__main__":
    main()
