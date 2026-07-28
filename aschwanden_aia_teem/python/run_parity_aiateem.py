"""run_parity_aiateem.py — Python twin driver for aschwanden aia_teem parity.

Injects oracle probes (logical orientation), builds the flux table (verify vs 01_flux),
runs the fit on the INJECTED oracle flux (so the discrete argmin is deterministic),
dumps twin probes. Usage: python run_parity_aiateem.py <oracle_dir> <py_dir>
"""
import os
import sys

import numpy as np
from scipy.io import readsav

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load(oracle_dir, pid):
    sav = readsav(os.path.join(oracle_dir, "probe_%s.sav" % pid), verbose=False)
    rec = np.atleast_1d(sav["chk"])[0]
    out = {}
    for n in rec.dtype.names:
        if n.lower() == "probe_id":
            continue
        v = np.asarray(rec[n])
        out[n.lower()] = np.transpose(v) if v.ndim >= 2 else v   # orientation: logical
    return out


def main():
    oracle_dir, py_dir = sys.argv[1], sys.argv[2]
    os.makedirs(py_dir, exist_ok=True)
    os.environ["CHK_DIR"] = py_dir
    from aia_teem import build_flux_table, fit
    from chk_dump import chk_dump

    r = load(oracle_dir, "00_resp")
    f = load(oracle_dir, "01_flux")
    inp = load(oracle_dir, "02_input")
    resp, telog, dte = r["resp"], r["telog"], r["dte"]
    flux_oracle, tsig = f["flux"], f["tsig_"]
    images, texp_ = inp["images"], inp["texp_"]

    chk_dump("00_resp", resp=resp, telog=telog, dte=dte)
    flux_py = build_flux_table(resp, telog, dte, tsig)
    chk_dump("01_flux", flux=flux_py, tsig_=tsig)          # verify table build vs oracle
    chk_dump("02_input", images=images, texp_=texp_)

    # fit on the injected oracle flux (isolates the discrete argmin from table float diffs)
    te_map, em_map, sig_map, chi_map = fit(images, texp_, flux_oracle, telog, tsig)
    chk_dump("99_final", te_map=te_map, em_map=em_map, sig_map=sig_map, chi_map=chi_map)
    print("shapes flux", flux_py.shape, "te_map", te_map.shape)
    print("te_map:\n", te_map)


if __name__ == "__main__":
    main()
