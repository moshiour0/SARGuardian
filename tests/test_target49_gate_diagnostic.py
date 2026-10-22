from __future__ import annotations
import sys
from pathlib import Path
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from gunw_reader import decode_gunw_mask
from tools.target49_gate_diagnostic import audit_gunw, audit_goff, decode_wrs_mask

def test_documented_decimal_wrs_mask_digits_are_audited_separately():
    raw=np.array([[11,12,99,111,112,0,10,1,255]],dtype=np.uint8)
    before=raw.copy()
    gates=decode_wrs_mask(raw)
    assert np.array_equal(raw,before)
    assert gates["both_subswaths_valid_ignoring_water"].tolist()==[[True,True,True,True,True,False,False,False,False]]
    assert gates["current_reader_mask_gate"].tolist()==[[True,True,True,False,False,False,False,False,False]]
    assert gates["water_flag_nonzero"].tolist()==[[False,False,False,True,True,False,False,False,False]]
    assert gates["fill_255"].tolist()==[[False,False,False,False,False,False,False,False,True]]

def test_source_aoi_reports_target9_outside_without_relabeling_it_as_raw_no_data(tmp_path):
    from pyproj import Transformer
    file=tmp_path/"synthetic_NISAR_L2_PR_GUNW.h5"
    tf=Transformer.from_crs(4326,32645,always_xy=True)
    x0,y0=tf.transform(85.47591,28.25460)
    # A wider synthetic grid makes the square target window intersect
    # SOURCE_RING even though Target 9's centre is outside it. This reproduces
    # the real GOFF case where one nearby in-AOI pixel must not qualify T9.
    xs=x0+(np.arange(13)-6)*100.0
    ys=y0+(6-np.arange(13))*100.0
    prefix="science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram"
    group=f"{prefix}/HH"
    with h5py.File(file,"w") as f:
        f.create_dataset(f"{group}/unwrappedPhase",data=np.ones((13,13),dtype=np.float32))
        f.create_dataset(f"{group}/coherenceMagnitude",data=np.full((13,13),0.9,dtype=np.float32))
        f.create_dataset(f"{group}/connectedComponents",data=np.ones((13,13),dtype=np.uint8))
        mask=np.full((13,13),11,dtype=np.uint8); mask[6,6]=111
        f.create_dataset(f"{group}/mask",data=mask)
        f.create_dataset(f"{prefix}/xCoordinates",data=xs)
        f.create_dataset(f"{prefix}/yCoordinates",data=ys)
        p=f.create_dataset(f"{prefix}/projection",data=np.array(32645,dtype=np.int32))
        p.attrs["epsg_code"]=32645
    result=audit_gunw(file,aoi="source",radius=6,coh_min=0.3)
    t9=result["targets"]["9"]; t4=result["targets"]["4"]
    assert t9["inside_grid_extent"] is True
    assert t9["inside_requested_aoi_polygon"] is False
    assert t9["point_eligible_for_requested_aoi"] is False
    assert t9["window_cell_count"]==169
    assert t9["has_finite_nonzero_phase_support"] is True
    assert t9["has_connected_component_support"] is True
    assert t9["has_current_reader_mask_support"] is True
    assert t9["has_nearby_aoi_gated_support"] is True
    assert t9["nearest_grid_cell_support"]["nearest_grid_cell_after_aoi"] is False
    assert t9["has_aoi_gated_support"] is False
    assert t9["has_target_point_support"] is False
    assert t9["cumulative_stage_counts"]["finite_nonzero_phase"]==169
    assert t9["cumulative_stage_counts"]["after_current_reader_mask_gate"]==168
    assert 0 < t9["cumulative_stage_counts"]["after_aoi"] < 168
    assert "mask_numeric_layout_sampled_rows" in result
    assert "mask_numeric_layout_exact" in result
    assert not any(key.startswith(" ") for key in result)
    assert t4["inside_grid_extent"] is False
    assert any("TARGET_9_OUTSIDE_REQUESTED_AOI" in warning for warning in result["warnings"])

def test_target4_is_explicitly_inside_source_polygon():
    from tools.target49_gate_diagnostic import TARGETS
    from gunw_reader import point_in_ring, AOIS
    t=TARGETS["4"]
    assert point_in_ring(t["lon"],t["lat"],AOIS["source"])
    assert not point_in_ring(TARGETS["9"]["lon"],TARGETS["9"]["lat"],AOIS["source"])


def test_uint32_product_mask_does_not_treat_bit24_as_decimal_water_code():
    """A 0x01000000 ionosphere-screen flag must be separate from low-byte digits."""
    raw = np.array([[
        11, 12, 111, 16777228, 16777216, 16777340, 255,
        0x0000010C, 0x0001000C, 0x0200000C,
    ]], dtype=np.uint32)
    before = raw.copy()
    corrected = decode_wrs_mask(raw, apply_iono=True)
    uncorrected = decode_wrs_mask(raw, apply_iono=False)
    assert np.array_equal(raw, before)
    assert corrected["land_valid_subswaths_low_byte"].tolist() == [[
        True, True, False, True, False, False, False, True, True, True
    ]]
    assert corrected["ionosphere_interpolated_mask"].tolist() == [[
        False, False, False, True, True, True, False, False, False, False
    ]]
    assert corrected["anomaly_flags_nonzero"].tolist() == [[
        False, False, False, False, False, False, False, True, True, False
    ]]
    assert corrected["reserved_bits_nonzero"].tolist() == [[
        False, False, False, False, False, False, False, False, False, True
    ]]
    assert corrected["current_reader_mask_gate"].tolist() == [[
        True, True, False, False, False, False, False, False, False, False
    ]]
    assert uncorrected["current_reader_mask_gate"].tolist() == [[
        True, True, False, True, False, False, False, False, False, False
    ]]


def test_exact_mask_bit_summary_comes_from_the_complete_global_histogram(tmp_path):
    file = tmp_path / "uint32-mask_NISAR_L2_PR_GUNW.h5"
    from pyproj import Transformer
    from tools.target49_gate_diagnostic import TARGETS
    # Put Target 4 at the centre of the grid so we can check both local support
    # and the full-file histogram. The raster itself is deliberately tiny.
    tf = Transformer.from_crs(4326, 32645, always_xy=True)
    x0, y0 = tf.transform(TARGETS["4"]["lon"], TARGETS["4"]["lat"])
    xs = x0 + (np.arange(13) - 6) * 80.0
    ys = y0 + (6 - np.arange(13)) * 80.0
    prefix = "science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram"
    group = f"{prefix}/HH"
    mask = np.full((13, 13), 11, dtype=np.uint32)
    mask[6, 6] = 0x0100000C
    mask[6, 7] = 0x0000010C
    with h5py.File(file, "w") as h5:
        h5.create_dataset(f"{group}/unwrappedPhase", data=np.ones((13, 13), np.float32))
        h5.create_dataset(f"{group}/coherenceMagnitude", data=np.full((13, 13), 0.9, np.float32))
        h5.create_dataset(f"{group}/connectedComponents", data=np.ones((13, 13), np.uint8))
        h5.create_dataset(f"{prefix}/mask", data=mask)
        h5.create_dataset(f"{prefix}/xCoordinates", data=xs)
        h5.create_dataset(f"{prefix}/yCoordinates", data=ys)
        p = h5.create_dataset(f"{prefix}/projection", data=np.array(32645, dtype=np.int32))
        p.attrs["epsg_code"] = 32645
    result = audit_gunw(file, aoi="source", radius=6, coh_min=0.3)
    exact = result["mask_numeric_layout_exact"]
    assert exact["coverage"].startswith("exact;")
    assert exact["dtype"] == "uint32"
    assert exact["bit_set_counts_lsb0"]["24"] == 1
    assert exact["bit_set_counts_lsb0"]["8"] == 1
    assert exact["low_byte_histogram_decimal"]["12"] == 2
    assert result["targets"]["4"]["has_current_reader_mask_support"] is True


def test_raw_diagnostic_mask_gate_matches_the_production_reader_for_fill_and_packed_bits():
    """The diagnostic must never drift from the exact production quality gate."""
    cases = [
        np.array([[0, 11, 12, 111, 255]], dtype=np.uint8),
        np.array([[
            11, 12, 111, 255, 0x0100000C, 0x0000010C,
            0x0001000C, 0x0200000C,
        ]], dtype=np.uint32),
    ]
    pairs = [
        ("water_flag_nonzero", "water_flag_nonzero"),
        ("ref_subswath_zero", "reference_subswath_zero"),
        ("secondary_subswath_zero", "secondary_subswath_zero"),
        ("both_subswaths_valid_ignoring_water", "both_subswaths_valid_ignoring_water"),
        ("land_valid_subswaths_low_byte", "land_valid_subswaths_low_byte"),
        ("anomaly_flags_nonzero", "anomaly_flags_nonzero"),
        ("ionosphere_interpolated_mask", "ionosphere_interpolated_mask"),
        ("reserved_bits_nonzero", "reserved_bits_nonzero"),
        ("current_reader_mask_gate", "current_reader_mask_gate"),
    ]
    for raw in cases:
        for apply_iono in (True, False):
            production = decode_gunw_mask(raw, apply_iono=apply_iono)
            diagnostic = decode_wrs_mask(raw, apply_iono=apply_iono)
            for diagnostic_key, production_key in pairs:
                assert np.array_equal(
                    diagnostic[diagnostic_key], production[production_key]
                ), (
                    f"{raw.dtype} mask decoder drift in {diagnostic_key}; "
                    f"diagnostic={diagnostic[diagnostic_key]}, "
                    f"production={production[production_key]}"
                )


def test_raw_gunw_diagnostic_rejects_unvalidated_mask_dtypes():
    import pytest
    with pytest.raises(TypeError, match="Unsupported GUNW mask dtype"):
        decode_wrs_mask(np.array([[12]], dtype=np.uint16))


def test_goff_target9_nearby_polygon_pixel_is_not_target_point_support(tmp_path):
    from pyproj import Transformer

    file = tmp_path / "synthetic_NISAR_L2_PR_GOFF.h5"
    tf = Transformer.from_crs(4326, 32645, always_xy=True)
    x0, y0 = tf.transform(85.47591, 28.25460)
    xs = x0 + (np.arange(13) - 6) * 100.0
    ys = y0 + (6 - np.arange(13)) * 100.0
    parent = "science/LSAR/GOFF/grids/frequencyA/pixelOffsets/HH"
    group = f"{parent}/layer1"
    with h5py.File(file, "w") as f:
        f.create_dataset(f"{group}/slantRangeOffset", data=np.ones((13, 13), dtype=np.float32))
        f.create_dataset(f"{group}/alongTrackOffset", data=np.ones((13, 13), dtype=np.float32))
        f.create_dataset(f"{group}/correlationSurfacePeak", data=np.full((13, 13), 0.9, dtype=np.float32))
        f.create_dataset(f"{group}/snr", data=np.full((13, 13), 5.0, dtype=np.float32))
        f.create_dataset(f"{group}/xCoordinates", data=xs)
        f.create_dataset(f"{group}/yCoordinates", data=ys)
        p = f.create_dataset(f"{parent}/projection", data=np.array(32645, dtype=np.int32))
        p.attrs["epsg_code"] = 32645

    result = audit_goff(file, aoi="source", radius=6, corr_min=0.3, snr_min=3.0)
    t9 = result["layers"][0]["targets"]["9"]
    assert t9["inside_grid_extent"] is True
    assert t9["inside_requested_aoi_polygon"] is False
    assert t9["cumulative_stage_counts"]["after_aoi"] > 0
    assert t9["has_nearby_aoi_gated_support"] is True
    assert t9["nearest_grid_cell_support"]["nearest_grid_cell_after_aoi"] is False
    assert t9["has_aoi_gated_support"] is False
    assert t9["has_target_point_support"] is False
