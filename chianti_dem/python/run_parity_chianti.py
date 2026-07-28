import os, sys
import numpy as np
from scipy.io import readsav
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

def load(oracle_dir, pid):
    rec = np.atleast_1d(readsav(os.path.join(oracle_dir, "probe_%s.sav" % pid), verbose=False)["chk"])[0]
    return {n.lower(): np.asarray(rec[n]) for n in rec.dtype.names if n.lower() != "probe_id"}

def main():
    oracle_dir, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True); os.environ["CHK_DIR"] = py_dir
    from chianti_dem_fit import DemFit
    a = load(oracle_dir, "00_input"); b = load(oracle_dir, "00b_grids")
    ch = np.transpose(a["ch_tot_contr"])   # (n_obs,nt)->(nt,n_obs)
    df = DemFit(a["obs_int"], a["obs_sig"], ch, a["log_t_mesh"], a["log_dem_temp"],
                b["dem_temp"], b["d_dem_temp"], b["log_dem_mesh"])
    mesh, y, chisqr, it = df.fit(flambda=10.0, scale=1.0, niter=20)
    print("iters", it, "chisqr", float(chisqr), "fitted mesh", np.array2string(mesh, precision=5))

if __name__ == "__main__":
    main()
