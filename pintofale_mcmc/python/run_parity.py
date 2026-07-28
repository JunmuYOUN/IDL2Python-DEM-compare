#!/usr/bin/env python
"""run_parity.py — CLI over parity_lib: compare converted mcmc_dem vs IDL oracle.

Usage:  python converted/run_parity.py [set1|set2|set3|all] [oracle_dir]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parity_lib import SETS, run_set        # noqa: E402

which = sys.argv[1] if len(sys.argv) > 1 else "all"
oracle = sys.argv[2] if len(sys.argv) > 2 else None
names = sorted(SETS) if which == "all" else [which]

fail = 0
for name in names:
    results, npass = run_set(name, oracle)
    for label, verd, det in results:
        if verd != "EXACT":
            print(f"  {verd:<7s} {label:<26s} {det}")
    bad = len(results) - npass
    fail += bad
    print(f"  [{name}] {npass}/{len(results)} bit-identical — {SETS[name]['desc']}")
    if bad:
        first = next(r for r in results if r[1] != "EXACT")
        print(f"  [{name}] FIRST DIVERGENCE: {first[0]} — {first[2]}")
sys.exit(1 if fail else 0)
