"""run_parity_eitdem.py — Python twin driver for the eit_dem numerical core.

Injects inputs from the IDL oracle probes (logical orientation), runs the port kernels
(e_fit2, eit_kcorr_iter, dem_map_expand), dumps output probes with matching names.
Usage: python run_parity_eitdem.py <oracle_dir> <py_dir>
"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from chk_dump import chk_dump          # noqa: E402
from eit_dem import e_fit2, eit_kcorr_iter, dem_map_expand   # noqa: E402


def load(oracle_dir, pid):
    rec = np.atleast_1d(readsav(os.path.join(oracle_dir, "probe_%s.sav" % pid), verbose=False)["chk"])[0]
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

    # ---- Part A: e_fit2 ----
    a = load(oracle_dir, "00_efit2_in")
    chk_dump("00_efit2_in", xt=a["xt"], y=a["y"])
    coefs_ef = e_fit2(a["xt"], a["y"])
    chk_dump("01_efit2_out", coefs_ef=coefs_ef)

    # ---- Part B: eit_kcorr /iter ----
    b = load(oracle_dir, "02_kcorr_in")
    b2 = load(oracle_dir, "02b_kcorr_in")
    chk_dump("02_kcorr_in", e171=b["e171"], e195=b["e195"], e284=b["e284"],
             e304=b["e304"], d171=b["d171"], d195=b["d195"])
    chk_dump("02b_kcorr_in", cool284=b2["cool284"], hot284=b2["hot284"],
             cool304=b2["cool304"], warm304=b2["warm304"], hot304=b2["hot304"], x=b2["x"])
    eitimgs = {"eit171": b["e171"], "eit195": b["e195"], "eit284": b["e284"], "eit304": b["e304"]}
    dn_maps = {"dn171": b["d171"], "dn195": b["d195"], "cool_284": b2["cool284"],
               "hot_284": b2["hot284"], "cool_304": b2["cool304"],
               "warm_304": b2["warm304"], "hot_304": b2["hot304"]}
    k171, k195, k284, k304 = eit_kcorr_iter(eitimgs, dn_maps, b2["x"])
    chk_dump("03_kcorr_out", k171=k171, k195=k195, k284=k284, k304=k304)

    # ---- Part C: DEM-map expansion ----
    c = load(oracle_dir, "04_expand_in")
    chk_dump("04_expand_in", coefs=c["coefs"], temp=c["temp"], dem=c["dem"], x=c["x"])
    dem_map, temp2 = dem_map_expand(c["coefs"], c["temp"], c["dem"], c["x"])
    chk_dump("99_dem_map", dem_map=dem_map, temp2=temp2)
    print("done: coefs_ef", coefs_ef.shape, "dem_map", dem_map.shape)


if __name__ == "__main__":
    main()
