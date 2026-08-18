import os, sys
import numpy as np
from scipy.io import readsav
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)

def load(oracle, pid):
    rec=np.atleast_1d(readsav(os.path.join(oracle,"probe_%s.sav"%pid),verbose=False)['chk'])[0]
    return {n.lower():np.asarray(rec[n]) for n in rec.dtype.names if n.lower()!='probe_id'}

def main():
    oracle, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True); os.environ["CHK_DIR"]=py_dir
    from aia_sparse_em import aia_sparse_em_init, aia_sparse_em_solve
    ini=load(oracle,"01_init"); inp=load(oracle,"00_input")
    tresp=np.transpose(ini["tresp"])           # (21,6)->(6,21)
    lgtaxis=ini["lgtaxis"]
    image=np.transpose(inp["image"])           # (6,1,1)->(1,1,6)
    Dict, basis=aia_sparse_em_init(tresp, lgtaxis)
    coeffs, oem, zmax, status=aia_sparse_em_solve(image, Dict, basis,
                                                  tolfac=float(inp["tolfac"]), eps=float(inp["eps"]))
    print("total(coeffs)=", float(coeffs.sum()), " peak oem=", float(oem.max()),
          " status=", float(status.ravel()[0]), " nonzero=", int((coeffs>0).sum()))
main()
