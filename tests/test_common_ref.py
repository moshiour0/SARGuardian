from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import common_ref


def _write_tif(path: Path, *, transform, values=None, crs="EPSG:32645", shape=(300, 300)):
    import rasterio
    from rasterio.transform import Affine

    data = values if values is not None else np.ones(shape, dtype=np.float32)
    if isinstance(transform, tuple):
        transform = Affine(*transform)
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path, "w", driver="GTiff", height=shape[0], width=shape[1],
        count=1, dtype="float32", crs=crs, transform=transform, nodata=-9999,
    ) as ds:
        ds.write(np.asarray(data, dtype=np.float32), 1)


def _target_csv(path: Path):
    path.write_text(
        "id,lat,lon,reading\n"
        "4,28.27799,85.52983,DETACHMENT-like\n"
        "9,28.25460,85.47591,DETACHMENT-like\n",
        encoding="utf-8",
    )


def test_common_reference_refuses_duplicate_acquisition_pairs_from_multiple_polarizations(
    tmp_path, monkeypatch, capsys
):
    from rasterio.transform import from_origin

    out = tmp_path / "export_goff_src"
    transform = from_origin(535000, 3155000, 80, 80)
    _write_tif(out / "GOFF_20260702_20260714_PR_HH-layer2.tif", transform=transform)
    _write_tif(out / "GOFF_20260702_20260714_PR_VV-layer2.tif", transform=transform)
    targets = tmp_path / "targets.csv"
    _target_csv(targets)
    monkeypatch.setattr(sys, "argv", [
        "common_ref.py", "--dir", str(out), "--match", "layer2",
        "--targets", str(targets),
    ])
    assert common_ref.main() == 1
    assert "multiple exports for the same acquisition pair" in capsys.readouterr().err.lower()


def test_common_reference_refuses_same_shape_rasters_with_different_transforms(
    tmp_path, monkeypatch, capsys
):
    from rasterio.transform import from_origin

    out = tmp_path / "export_goff_src"
    t1 = from_origin(535000, 3155000, 80, 80)
    t2 = from_origin(535080, 3155000, 80, 80)
    _write_tif(out / "GOFF_20260702_20260714_PR_HH-layer2.tif", transform=t1)
    _write_tif(out / "GOFF_20260714_20260726_PR_HH-layer2.tif", transform=t2)
    targets = tmp_path / "targets.csv"
    _target_csv(targets)
    monkeypatch.setattr(sys, "argv", [
        "common_ref.py", "--dir", str(out), "--match", "HH-layer2",
        "--targets", str(targets),
    ])
    assert common_ref.main() == 1
    assert "grid mismatch" in capsys.readouterr().err.lower()


def test_common_reference_excludes_out_of_aoi_candidate_and_uses_exact_polygon(
    tmp_path, monkeypatch
):
    from rasterio.transform import from_origin

    out = tmp_path / "export_goff_src"
    transform = from_origin(535000, 3155000, 80, 80)
    rng = np.random.default_rng(4)
    # Full rectangular windows intentionally contain finite cells outside the
    # AOI. The common-reference algorithm must explicitly clip to SOURCE_RING.
    for name, offset in [
        ("GOFF_20260702_20260714_PR_HH-layer2.tif", 0.0),
        ("GOFF_20260714_20260726_PR_HH-layer2.tif", 2.0),
    ]:
        values = rng.normal(0, 1, (300, 300)).astype(np.float32) + offset
        _write_tif(out / name, transform=transform, values=values)
    targets = tmp_path / "targets.csv"
    _target_csv(targets)
    result = tmp_path / "common_ref.csv"
    monkeypatch.setattr(sys, "argv", [
        "common_ref.py", "--dir", str(out), "--match", "HH-layer2",
        "--targets", str(targets), "--reading", "DETACHMENT-like",
        "--buffer-km", "2", "--csv", str(result),
    ])
    assert common_ref.main() == 0
    with result.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    assert {row["target_id"] for row in rows} == {"4"}
