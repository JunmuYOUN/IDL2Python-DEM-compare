"""Download and convert the tabulated SDO/AIA temperature response.

Run this with the repository's ``ssw`` Conda environment::

    conda run -n ssw python 1_Resp/fetch_aia_temperature_response.py

The upstream ASDF file is used by the SunPy/aiapy-based synthesizAR package and
was generated from SolarSoft ``aia_get_response(/temp, /dn)``.  A compact NPZ
copy is created so the DEM notebook only needs NumPy at runtime.
"""

from __future__ import annotations

import hashlib
from importlib.metadata import version
from pathlib import Path
from urllib.request import urlopen

import asdf
import numpy as np


SOURCE_COMMIT = "77aab1767e25bc10a200e1b26270da1de20922c2"
SOURCE_URL = (
    "https://raw.githubusercontent.com/wtbarnes/synthesizAR/"
    f"{SOURCE_COMMIT}/synthesizAR/instruments/data/aia_temperature_response.asdf"
)
EXPECTED_SHA256 = "aec75a29b2abe108a086932e5dd6a098cba628d69e8fc3ffa2e8ecfff4708740"
CHANNELS = np.array([94, 131, 171, 193, 211, 335], dtype=np.int64)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _download_verified(destination: Path) -> str:
    with urlopen(SOURCE_URL, timeout=60) as response:
        payload = response.read()
    digest = _sha256(payload)
    if digest != EXPECTED_SHA256:
        raise RuntimeError(
            f"Downloaded ASDF SHA-256 mismatch: expected {EXPECTED_SHA256}, got {digest}"
        )
    destination.write_bytes(payload)
    return digest


def _quantity_value(tree, key: str) -> np.ndarray:
    value = tree[key]
    if hasattr(value, "value"):
        value = value.value
    elif isinstance(value, dict) and "value" in value:
        value = value["value"]
    return np.asarray(value, dtype=np.float64)


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    asdf_path = output_dir / "aia_temperature_response.asdf"
    npz_path = output_dir / "aia_temperature_response.npz"

    digest = _download_verified(asdf_path)
    with asdf.open(asdf_path, mode="r", memmap=False, lazy_load=False) as af:
        temperature_k = _quantity_value(af.tree, "temperature")
        response = np.column_stack([_quantity_value(af.tree, str(ch)) for ch in CHANNELS])
        response_unit = str(getattr(af.tree["94"], "unit", "cm5 DN / (pix s)"))

    logt = np.log10(temperature_k)
    if temperature_k.shape != (101,) or response.shape != (101, CHANNELS.size):
        raise ValueError(
            f"Unexpected response dimensions: temperature={temperature_k.shape}, response={response.shape}"
        )
    if not np.all(np.isfinite(response)) or np.any(response < 0):
        raise ValueError("Temperature response must be finite and non-negative")
    if not np.all(np.diff(logt) > 0):
        raise ValueError("Temperature grid must be strictly increasing")

    np.savez_compressed(
        npz_path,
        temperature_k=temperature_k,
        logt=logt,
        wavelengths=CHANNELS,
        response=response,
        response_unit=np.array(response_unit),
        source_url=np.array(SOURCE_URL),
        source_commit=np.array(SOURCE_COMMIT),
        source_sha256=np.array(digest),
        upstream_generation=np.array("SolarSoft aia_get_response(/temp,/dn)"),
        time_dependent_correction=np.array(False),
        eve_normalization=np.array(False),
        sunpy_version=np.array(version("sunpy")),
        aiapy_version=np.array(version("aiapy")),
        asdf_version=np.array(version("asdf")),
    )

    print(f"Saved original: {asdf_path} ({asdf_path.stat().st_size:,} bytes)")
    print(f"Saved portable: {npz_path} ({npz_path.stat().st_size:,} bytes)")
    print(f"SHA-256: {digest}")
    print(f"Grid: log10(T/K)={logt[0]:.2f}..{logt[-1]:.2f}, {logt.size} points")
    print(f"Response shape: {response.shape}, unit={response_unit}")


if __name__ == "__main__":
    main()
