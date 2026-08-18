"""run_parity_firdem.py — Python twin driver for firdem (Plowman fast DEM) core parity.

Injects the a_struc / mapping-matrix scaffolding from the oracle probes (library
boundary), runs the verified core (triangle basis + regularize_data + first-pass DEM),
and dumps twin probes. Usage: python run_parity_firdem.py <oracle> <py_dir>
"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("CHK_DIR", "")


def load(oracle_dir, pid):
    rec = np.atleast_1d(readsav(os.path.join(oracle_dir, "probe_%s.sav" % pid),
                                verbose=False)["chk"])[0]
    out = {}
    for n in rec.dtype.names:
        if n.lower() == "probe_id":
            continue
        v = np.asarray(rec[n])
        out[n.lower()] = np.transpose(v) if v.ndim >= 2 else v   # logical orientation
    return out


def main():
    oracle_dir, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True)
    os.environ["CHK_DIR"] = py_dir
    from chk_dump import chk_dump
    from firdem import firdem_triangle_basis, firdem_regularize_data, first_pass_dem

    d00 = load(oracle_dir, "00_input")
    d01a = load(oracle_dir, "01a_astruc")
    d01b = load(oracle_dir, "01b_astruc")
    d02 = load(oracle_dir, "02_setup")
    d05 = load(oracle_dir, "05_a2array")

    # --- injected boundary (redump for trivial cross-check) ---
    chk_dump("00_input", dns0=d00["dns0"], errs=d00["errs"], exptimes=d00["exptimes"],
             tresps_arr=d00["tresps_arr"], logt=d00["logt"])
    chk_dump("01a_astruc", t=d01a["t"], tr=d01a["tr"], tr_norm=d01a["tr_norm"], basis=d01a["basis"])
    chk_dump("01b_astruc", a_array=d01b["a_array"], wvec=d01b["wvec"],
             a_inv=d01b["a_inv"], flat_coffs=d01b["flat_coffs"])
    chk_dump("05_a2array", a2_array=d05["a2_array"])

    # --- verified: triangle basis (float32-grid limited) ---
    nb2 = d02["basis2"].shape[1]
    logt2 = np.log10(np.asarray(d01a["t"], dtype=np.float64))
    basis2, t0a = firdem_triangle_basis(nb2, logt2)
    chk_dump("02_setup", basis2=basis2, t0a=t0a, basis22=d02["basis22"],
             normfac=d02["normfac"], chi2_reg_ends=d02["chi2_reg_ends"])

    # --- verified: regularized data (double precision inversion) ---
    datavec0 = d00["dns0"][0, 0, :]
    sigmas = d00["errs"][0, 0, :]
    tr_norm = d02["normfac"]                       # a_struc.tr_norm*exptimes
    a_inv = d01b["a_inv"]
    chi2_end = float(d02["chi2_reg_ends"][0])
    data_out = firdem_regularize_data(datavec0, sigmas, tr_norm, a_inv, chi2_end=chi2_end)
    chk_dump("03_regdata", data_out=data_out)

    # --- verified: first-pass DEM assembly ---
    dem_initial = first_pass_dem(data_out, a_inv, d02["normfac"], d02["basis22"])
    chk_dump("04_firstpass", dem_initial=dem_initial)

    print("data_out", np.array2string(data_out, precision=4),
          "dem_initial [%.4g, %.4g]" % (dem_initial.min(), dem_initial.max()))


if __name__ == "__main__":
    main()
