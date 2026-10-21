from __future__ import annotations
import sys
from pathlib import Path
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from tools.target49_gate_diagnostic import audit_gunw, decode_wrs_mask

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
    xs=x0+(np.arange(13)-6)*80.0
    ys=y0+(6-np.arange(13))*80.0
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
    assert t9["window_cell_count"]==169
    assert t9["cumulative_stage_counts"]["finite_nonzero_phase"]==169
    assert t9["cumulative_stage_counts"]["after_current_reader_mask_gate"]==168
    assert t9["cumulative_stage_counts"]["after_aoi"] < 168
    assert t4["inside_grid_extent"] is False
    assert any("TARGET_9_OUTSIDE_REQUESTED_AOI" in warning for warning in result["warnings"])

def test_target4_is_explicitly_inside_source_polygon():
    from tools.target49_gate_diagnostic import TARGETS
    from gunw_reader import point_in_ring, AOIS
    t=TARGETS["4"]
    assert point_in_ring(t["lon"],t["lat"],AOIS["source"])
    assert not point_in_ring(TARGETS["9"]["lon"],TARGETS["9"]["lat"],AOIS["source"])
