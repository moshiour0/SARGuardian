#!/usr/bin/env python3
"""Read-only raw NISAR HDF5 gate diagnostic for SARGuardian targets 4 and 9.

This tool does not alter science outputs or relax any threshold. It separately
reports native-grid support, quality-gate counts and requested-AOI eligibility.
"""
from __future__ import annotations
import argparse, hashlib, json, math, re, sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gunw_reader import AOIS, build_aoi_mask, point_in_ring, walk
from goff_reader import layer_groups

TARGETS = {
    "4": {"lon": 85.52983, "lat": 28.27799},
    "9": {"lon": 85.47591, "lat": 28.25460},
}
CHUNK_ROWS = 128

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def decode_wrs_mask(raw, *, apply_iono=True):
    """Decode packed GUNW mask fields without converting the whole word to decimal."""
    raw = np.asarray(raw)
    if not np.issubdtype(raw.dtype, np.integer):
        raise TypeError(f"GUNW mask must have integer dtype; got {raw.dtype}")
    m = raw.astype(np.uint64, copy=False)
    low = m & np.uint64(0xFF)
    fill = m == 255
    water = (low // 100) % 10
    ref_sub = (low // 10) % 10
    sec_sub = low % 10
    anomaly = (m & np.uint64(0x00FFFF00)) != 0
    iono_interpolated = (m & np.uint64(0x01000000)) != 0
    reserved = (m & np.uint64(0xFE000000)) != 0
    both = (ref_sub > 0) & (sec_sub > 0) & ~fill
    land_subswaths = (water == 0) & both
    usable = land_subswaths & ~anomaly & ~reserved
    if apply_iono:
        usable &= ~iono_interpolated
    return {
        "fill_255": fill,
        "water_flag_nonzero": (water != 0) & ~fill,
        "ref_subswath_zero": (ref_sub == 0) & ~fill,
        "secondary_subswath_zero": (sec_sub == 0) & ~fill,
        "both_subswaths_valid_ignoring_water": both,
        "land_valid_subswaths_low_byte": land_subswaths,
        "anomaly_flags_nonzero": anomaly,
        "ionosphere_interpolated_mask": iono_interpolated,
        "reserved_bits_nonzero": reserved,
        "current_reader_mask_gate": usable,
    }

def matching(datasets, group, parent, names, shape=None, ndim=2):
    lower = {p.lower(): p for p in datasets}
    for name in names:
        for base in (group, parent):
            p = lower.get(f"{base}/{name}".lower())
            if p and (shape is None or datasets[p][0] == shape) and (ndim is None or len(datasets[p][0]) == ndim):
                return p
    return None

def epsg_from(h5, datasets, group, parent):
    p = matching(datasets, group, parent, ("projection",), shape=None, ndim=None)
    if not p: return None
    try:
        obj = h5[p]
        v = obj.attrs.get("epsg_code")
        return int(np.ravel(v)[0]) if v is not None else int(np.ravel(obj[()])[0])
    except (ValueError, TypeError, OSError, IndexError):
        return None

def projected_targets(xs, ys, epsg, aoi, radius):
    if epsg is None:
        raise ValueError("Product EPSG metadata missing; refusing to guess grid CRS")
    from pyproj import Transformer
    tf = Transformer.from_crs(4326, epsg, always_xy=True)
    out = {}
    for tid, t in TARGETS.items():
        x, y = tf.transform(t["lon"], t["lat"])
        col, row = int(np.argmin(np.abs(xs - x))), int(np.argmin(np.abs(ys - y)))
        out[tid] = {
            "lat": t["lat"], "lon": t["lon"], "row": row, "col": col,
            "radius_px": radius, "window_shape": [2 * radius + 1, 2 * radius + 1],
            "nearest_grid_cell_distance_m": round(math.hypot(float(xs[col]-x), float(ys[row]-y)), 3),
            "inside_grid_extent": bool(xs.min() <= x <= xs.max() and ys.min() <= y <= ys.max()),
            "inside_requested_aoi_polygon": bool(point_in_ring(t["lon"], t["lat"], AOIS[aoi])),
        }
    return out

def window_slice(center, radius, length):
    return slice(max(0, center-radius), min(length, center+radius+1))

def hist_update(counter, values):
    vals, counts = np.unique(values, return_counts=True)
    for v, c in zip(vals.tolist(), counts.tolist()): counter[str(v)] += int(c)

def make_aoi_mask(xs, ys, aoi, epsg):
    m = build_aoi_mask(xs, ys, AOIS[aoi], epsg)
    if m is None: raise ValueError("Cannot build AOI mask; check pyproj and the product EPSG code")
    return m

def _json_safe(value):
    """Convert HDF5 attribute values to JSON while preserving scalar/array data."""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, np.ndarray):
        return [_json_safe(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _dataset_metadata(h5, path):
    if not path or path not in h5:
        return None
    ds = h5[path]
    return {
        "path": path,
        "dtype": str(ds.dtype),
        "shape": list(ds.shape),
        "attrs": {str(k): _json_safe(v) for k, v in ds.attrs.items()},
        "parent_group_attrs": {str(k): _json_safe(v) for k, v in ds.parent.attrs.items()},
    }


def _mask_numeric_layout(mask):
    """Report raw number decomposition only; does not assign scientific meaning."""
    flat = np.asarray(mask).astype(np.uint64, copy=False).ravel()
    low_byte = flat & np.uint64(0xFF)
    high_word = flat & np.uint64(0xFFFFFF00)
    return {
        "dtype": str(np.asarray(mask).dtype),
        "high_word_histogram_hex": {
            f"0x{int(k):08x}": int(v) for k, v in sorted(Counter(high_word.tolist()).items())
        },
        "low_byte_histogram_decimal": {
            str(k): int(v) for k, v in sorted(Counter(low_byte.tolist()).items())
        },
        "bit_set_counts_lsb0": {
            str(bit): int(((flat >> np.uint64(bit)) & np.uint64(1)).sum())
            for bit in range(32)
            if np.any((flat >> np.uint64(bit)) & np.uint64(1))
        },
        "warning": "Numeric decomposition only; do not infer mask semantics without product-format verification.",
    }


def _mask_numeric_layout_from_hist(histogram, dtype):
    """Exact bit/byte summary reconstructed from the complete raw-value histogram."""
    high, low, bits = Counter(), Counter(), Counter()
    for raw, count in histogram.items():
        v, n = int(raw), int(count)
        high[f"0x{(v & 0xFFFFFF00):08x}"] += n
        low[str(v & 0xFF)] += n
        for bit in range(32):
            if v & (1 << bit):
                bits[str(bit)] += n
    return {
        "dtype": str(dtype),
        "high_word_histogram_hex": dict(sorted(high.items())),
        "low_byte_histogram_decimal": dict(sorted(low.items(), key=lambda kv: int(kv[0]))),
        "bit_set_counts_lsb0": dict(sorted(bits.items(), key=lambda kv: int(kv[0]))),
        "coverage": "exact; derived from full-file raw mask histogram",
        "warning": "Numeric decomposition only; no bit/byte semantics inferred.",
    }


def audit_gunw(path, aoi, radius, coh_min, apply_iono=True):
    result_type = "GUNW"
    with h5py.File(path, "r") as h5:
        d = walk(h5)
        pphase = next((p for p,(s,_) in sorted(d.items(), key=lambda kv: (0 if "frequencya" in kv[0].lower() else 1, 0 if re.search(r"/(hh|vv)/",kv[0].lower()) else 1, len(kv[0]))) if "unwrappedphase" in p.lower() and len(s)==2), None)
        if not pphase: raise ValueError("2-D unwrappedPhase dataset not found")
        group, parent = pphase.rsplit("/",1)[0], pphase.rsplit("/",1)[0].rsplit("/",1)[0]
        shape = d[pphase][0]
        coh = matching(d, group, parent, ("coherenceMagnitude","coherence"), shape)
        comp = matching(d, group, parent, ("connectedComponents",), shape)
        mask = matching(d, group, parent, ("mask",), shape)
        xp = matching(d, group, parent, ("xCoordinates",), ndim=1)
        yp = matching(d, group, parent, ("yCoordinates",), ndim=1)
        epsg = epsg_from(h5, d, group, parent)
        if not xp or not yp: raise ValueError(f"Coordinate vectors missing for {group}")
        xs, ys = np.asarray(h5[xp][()], dtype=float).ravel(), np.asarray(h5[yp][()], dtype=float).ravel()
        if (len(ys),len(xs)) != shape: raise ValueError(f"Phase shape {shape} != coordinate-vector shape {(len(ys),len(xs))}")
        targets = projected_targets(xs,ys,epsg,aoi,radius)
        counts, mhist = Counter(), Counter()
        tcounts = {t:Counter() for t in targets}
        tmaskhist = {t:Counter() for t in targets}
        aoi_cells = 0
        phase_ds, coh_ds = h5[pphase], (h5[coh] if coh else None)
        comp_ds, mask_ds = (h5[comp] if comp else None), (h5[mask] if mask else None)
        for y0 in range(0,shape[0],CHUNK_ROWS):
            y1 = min(shape[0],y0+CHUNK_ROWS); rs = slice(y0,y1)
            phase = np.asarray(phase_ds[rs,:],dtype=float)
            stage = np.isfinite(phase) & (phase != 0)
            valid = stage.copy()
            counts["finite_nonzero_phase"] += int(stage.sum())
            if coh_ds is not None:
                cv = np.asarray(coh_ds[rs,:],dtype=float) >= coh_min
                valid &= cv; counts[f"after_coherence_ge_{coh_min:g}"] += int(valid.sum())
            else: counts["coherence_dataset_missing"] += 1
            if comp_ds is not None:
                valid &= np.asarray(comp_ds[rs,:]) > 0
                counts["after_connected_components_gt_0"] += int(valid.sum())
            else: counts["connected_components_dataset_missing"] += 1
            raw_mask = np.asarray(mask_ds[rs,:]) if mask_ds is not None else None
            gates = decode_wrs_mask(raw_mask, apply_iono=apply_iono) if raw_mask is not None and raw_mask.shape == phase.shape else None
            if gates is not None:
                hist_update(mhist,raw_mask)
                before = valid.copy()
                counts["after_mask_ignoring_water"] += int((before & gates["both_subswaths_valid_ignoring_water"]).sum())
                valid &= gates["land_valid_subswaths_low_byte"]
                counts["after_land_valid_subswaths_low_byte"] += int(valid.sum())
                valid &= ~gates["anomaly_flags_nonzero"]
                counts["after_anomaly_flags_clear"] += int(valid.sum())
                if apply_iono:
                    valid &= ~gates["ionosphere_interpolated_mask"]
                    counts["after_iono_interpolated_mask_excluded"] += int(valid.sum())
                else:
                    counts["iono_interpolated_pixels_retained_no_iono"] += int((valid & gates["ionosphere_interpolated_mask"]).sum())
                valid &= ~gates["reserved_bits_nonzero"]
                counts["after_reserved_bits_clear"] += int(valid.sum())
                counts["after_current_reader_mask_gate"] += int(valid.sum())
            else: counts["mask_missing_or_shape_mismatch"] += 1
            am = make_aoi_mask(xs,ys[rs],aoi,epsg)
            counts["after_aoi"] += int((valid & am).sum())
            aoi_cells += int(am.sum())
            for tid,t in targets.items():
                if not t["inside_grid_extent"]: continue
                rr = window_slice(t["row"],radius,shape[0]); cc=window_slice(t["col"],radius,shape[1])
                start,stop=max(rr.start,y0),min(rr.stop,y1)
                if start>=stop: continue
                local_sl = slice(start-y0,stop-y0); aoi_sl=am[local_sl,cc]
                pv=stage[local_sl,cc].copy()
                tc=tcounts[tid]; tc["window_cells"] += int(pv.size); tc["finite_nonzero_phase"] += int(pv.sum())
                if coh_ds is not None:
                    pv &= np.asarray(coh_ds[rs,:][local_sl,cc],dtype=float) >= coh_min
                    tc[f"after_coherence_ge_{coh_min:g}"] += int(pv.sum())
                if comp_ds is not None:
                    pv &= np.asarray(comp_ds[rs,:][local_sl,cc]) > 0
                    tc["after_connected_components_gt_0"] += int(pv.sum())
                if gates is not None:
                    gm={k:v[local_sl,cc] for k,v in gates.items()}
                    tc["after_mask_ignoring_water"] += int((pv & gm["both_subswaths_valid_ignoring_water"]).sum())
                    pv &= gm["land_valid_subswaths_low_byte"]
                    tc["after_land_valid_subswaths_low_byte"] += int(pv.sum())
                    pv &= ~gm["anomaly_flags_nonzero"]
                    tc["after_anomaly_flags_clear"] += int(pv.sum())
                    if apply_iono:
                        pv &= ~gm["ionosphere_interpolated_mask"]
                        tc["after_iono_interpolated_mask_excluded"] += int(pv.sum())
                    else:
                        tc["iono_interpolated_pixels_retained_no_iono"] += int((pv & gm["ionosphere_interpolated_mask"]).sum())
                    pv &= ~gm["reserved_bits_nonzero"]
                    tc["after_reserved_bits_clear"] += int(pv.sum())
                    tc["after_current_reader_mask_gate"] += int(pv.sum())
                    hist_update(tmaskhist[tid],raw_mask[local_sl,cc])
                pv &= aoi_sl; tc["after_aoi"] += int(pv.sum())
        target_obj={}
        for tid,t in targets.items():
            target_obj[tid]={**t,"window_cell_count":int(tcounts[tid]["window_cells"]),
                "cumulative_stage_counts":{k:int(v) for k,v in tcounts[tid].items() if k!="window_cells"},
                "raw_mask_value_histogram":dict(tmaskhist[tid])}
        attrs={}
        for key in ("identification/granuleName","identification/productVersion","identification/productLevel"):
            if key in h5:
                try:
                    raw=h5[key][()]; attrs[key]=raw.decode("utf-8","replace") if isinstance(raw,bytes) else str(np.asarray(raw).item())
                except Exception: pass
        return {"file":path.name,"type":result_type,"sha256":sha256(path),"bytes":path.stat().st_size,
            "crid":next((x for x in re.findall(r"_([A-Z]\d{5})_",path.name)),None),
            "phase_dataset":pphase,"coherence_dataset":coh,"components_dataset":comp,"mask_dataset":mask,
            "mask_dtype":str(mask_ds.dtype) if mask_ds is not None else None,
            "mask_dataset_metadata":_dataset_metadata(h5, mask),
            "           mask_numeric_layout_sampled_rows":_mask_numeric_layout(mask_ds[::CHUNK_ROWS, :] if mask_ds is not None else np.array([], dtype=np.uint8)),
            "phase_dataset_metadata":_dataset_metadata(h5, pphase),
            "coherence_dataset_metadata":_dataset_metadata(h5, coh),
            "components_dataset_metadata":_dataset_metadata(h5, comp),
            "x_dataset_metadata":_dataset_metadata(h5, xp),
            "y_dataset_metadata":_dataset_metadata(h5, yp),
            "x_dataset":xp,"y_dataset":yp,"epsg":epsg,"shape":list(shape),"attributes":attrs,
            "thresholds":{"coherence_min":coh_min,"radius_px":radius,"requested_aoi":aoi,"apply_iono":apply_iono},
            "global_stage_counts":dict(counts),"global_source_aoi_grid_cells":aoi_cells,
            "global_raw_mask_histogram":dict(mhist),"targets":target_obj,
            "warnings":[
                f"TARGET_{tid}_OUTSIDE_REQUESTED_AOI: finite cells after AOI clipping are not evidence of target-point support."
                for tid,t in targets.items() if not t["inside_requested_aoi_polygon"]],
            "interpretation_limits":["Gate counts are diagnostic, not physical validation.",
              "Mask hypotheses are reported side-by-side; do not change masks based on counts alone.",
              "Products from different frame IDs are not a matched footprint control."]}

def audit_goff(path, aoi, radius, corr_min, snr_min):
    with h5py.File(path,"r") as h5:
        d=walk(h5); groups=layer_groups(d)
        if not groups: raise ValueError("No pixelOffsets layers found")
        layers=[]
        for group in groups:
            names={p.rsplit("/",1)[-1]:p for p in d if p.startswith(group+"/")}
            rngp=names.get("slantRangeOffset"); azip=names.get("alongTrackOffset")
            corrp=names.get("correlationSurfacePeak"); snrp=names.get("snr")
            xp,yp=names.get("xCoordinates"),names.get("yCoordinates")
            if not rngp or not xp or not yp: continue
            shape=d[rngp][0]; xs=np.asarray(h5[xp][()],float).ravel(); ys=np.asarray(h5[yp][()],float).ravel()
            epsg=epsg_from(h5,d,group,group.rsplit("/",1)[0])
            targets=projected_targets(xs,ys,epsg,aoi,radius); c=Counter(); tc={k:Counter() for k in targets}; aoicells=0
            for y0 in range(0,shape[0],CHUNK_ROWS):
                y1=min(shape[0],y0+CHUNK_ROWS); rs=slice(y0,y1)
                rng=np.asarray(h5[rngp][rs,:],float); valid=np.isfinite(rng)
                c["finite_range_offset"]+=int(valid.sum())
                if azip:
                    azi=np.asarray(h5[azip][rs,:],float); valid &= np.isfinite(azi)
                    c["after_finite_azimuth"]+=int(valid.sum())
                else: azi=None
                corr=np.asarray(h5[corrp][rs,:],float) if corrp else None
                if corr is not None:
                    valid &= corr>=corr_min; c[f"after_correlation_ge_{corr_min:g}"]+=int(valid.sum())
                else: c["correlation_dataset_missing"]+=1
                snr=np.asarray(h5[snrp][rs,:],float) if snrp else None
                if snr is not None:
                    valid &= snr>=snr_min; c[f"after_snr_ge_{snr_min:g}"]+=int(valid.sum())
                else: c["snr_dataset_missing"]+=1
                am=make_aoi_mask(xs,ys[rs],aoi,epsg); c["after_aoi"]+=int((valid&am).sum()); aoicells+=int(am.sum())
                for tid,t in targets.items():
                    if not t["inside_grid_extent"]: continue
                    rr=window_slice(t["row"],radius,shape[0]); cc=window_slice(t["col"],radius,shape[1])
                    start,stop=max(rr.start,y0),min(rr.stop,y1)
                    if start>=stop: continue
                    sl=slice(start-y0,stop-y0); win=(sl,cc)
                    vv=np.isfinite(rng[win]); tc[tid]["window_cells"]+=int(vv.size); tc[tid]["finite_range_offset"]+=int(vv.sum())
                    if azi is not None: vv &= np.isfinite(azi[win]); tc[tid]["after_finite_azimuth"]+=int(vv.sum())
                    if corr is not None: vv &= corr[win]>=corr_min; tc[tid][f"after_correlation_ge_{corr_min:g}"]+=int(vv.sum())
                    if snr is not None: vv &= snr[win]>=snr_min; tc[tid][f"after_snr_ge_{snr_min:g}"]+=int(vv.sum())
                    vv &= am[sl,cc]; tc[tid]["after_aoi"]+=int(vv.sum())
            layers.append({"group":group,"range_dataset":rngp,"azimuth_dataset":azip,"correlation_dataset":corrp,"snr_dataset":snrp,
                "shape":list(shape),"epsg":epsg,"global_stage_counts":dict(c),"global_aoi_grid_cells":aoicells,
                "targets":{k:{**v,"window_cell_count":int(tc[k]["window_cells"]),
                    "cumulative_stage_counts":{s:int(n) for s,n in tc[k].items() if s!="window_cells"}} for k,v in targets.items()},
                "warning":"GOFF is distributed unculled; finite/correlation/SNR gates do not prove physical correctness."})
        if not layers: raise ValueError("No GOFF layer with offsets and coordinate vectors was audited")
        return {"file":path.name,"type":"GOFF","sha256":sha256(path),"bytes":path.stat().st_size,"layers":layers,
            "thresholds":{"correlation_min":corr_min,"snr_min":snr_min,"radius_px":radius,"requested_aoi":aoi}}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("files",nargs="+",type=Path)
    ap.add_argument("--aoi",choices=sorted(AOIS),default="source")
    ap.add_argument("--radius",type=int,default=6)
    ap.add_argument("--coherence",type=float,default=0.3)
    ap.add_argument("--no-iono",action="store_true",help="diagnose without ionosphere correction; retain the screen-only mask bit")
    ap.add_argument("--correlation",type=float,default=0.3)
    ap.add_argument("--snr",type=float,default=3.0)
    ap.add_argument("--out",type=Path,default=Path("outputs/target49_gate_diagnostic.json"))
    args=ap.parse_args()
    if args.radius<0: ap.error("--radius must be non-negative")
    report={"schema_version":1,"created_utc":datetime.now(timezone.utc).isoformat(),"aoi":args.aoi,
        "target_coordinates":TARGETS,"limits":["Diagnostic stage counts are not validation of geophysical truth."],"results":[],"errors":[]}
    for path in args.files:
        try:
            item=audit_gunw(path,args.aoi,args.radius,args.coherence,apply_iono=not args.no_iono) if "GUNW" in path.name.upper() else audit_goff(path,args.aoi,args.radius,args.correlation,args.snr)
            report["results"].append(item); print("PASS diagnostic",path.name,item["type"],item["sha256"])
        except Exception as exc:
            report["errors"].append({"file":path.name,"type":type(exc).__name__,"message":str(exc)})
            print(f"ERROR {path}: {type(exc).__name__}: {exc}",file=sys.stderr)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(f"Wrote {args.out}; success={len(report['results'])}; errors={len(report['errors'])}")
    return 1 if report["errors"] or not report["results"] else 0

if __name__=="__main__":
    raise SystemExit(main())
