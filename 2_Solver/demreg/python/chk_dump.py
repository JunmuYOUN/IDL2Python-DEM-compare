"""Parity probe dumper (Python side). Twin of tools/idl/chk_dump.pro.

Usage — insert at points semantically identical to the IDL probes:

    from chk_dump import chk_dump
    chk_dump('05_ratio_masks', mas=mas, msk=msk, mak=mak)

No-op unless env CHK_DIR is set, so calls may remain in shipped code.
Never deletes anything (overwrites the probe file for the same id only).
"""
import os

import numpy as np


def chk_dump(probe_id, **variables):
    """Dump named variables to $CHK_DIR/probe_<probe_id>.npz (compressed)."""
    d = os.environ.get("CHK_DIR", "")
    if not d:
        return
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "probe_%s.npz" % probe_id)
    arrays = {}
    for name, value in variables.items():
        arrays[name.lower()] = np.asarray(value)
    np.savez_compressed(path, **arrays)
    print("chk_dump: wrote %s" % path)


def load_probe(path):
    """Load a probe .npz back as {name: array} (for injection use-cases)."""
    out = {}
    with np.load(path, allow_pickle=False) as z:
        for k in z.files:
            out[k.lower()] = z[k]
    return out
