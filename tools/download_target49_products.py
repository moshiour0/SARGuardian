#!/usr/bin/env python3
"""Download exact NISAR HDF5 granules logged by SARGuardian-backend Run #55.

Uses ASF's supported Earthdata Login/.netrc session handling. Never put a
password or token on the command line, in this script, or in the output manifest.
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys
from datetime import datetime, timezone
from pathlib import Path

RUN55_GRANULES = {
    "GUNW-006-BETA-control": "NISAR_L2_PR_GUNW_009_098_A_016_010_4000_SH_20260103T233921_20260103T233956_20260115T233922_20260115T233956_X05010_N_F_J_001.h5",
    "GUNW-007-PROVISIONAL": "NISAR_L2_PR_GUNW_023_098_A_016_024_4000_SH_20260620T233921_20260620T233956_20260702T233920_20260702T233955_P05023_N_F_J_001.h5",
    "GUNW-009-PROVISIONAL": "NISAR_L2_PR_GUNW_024_098_A_016_025_4000_SH_20260702T233920_20260702T233955_20260714T233920_20260714T233954_P05023_F_F_J_001.h5",
    "GUNW-010-PROVISIONAL": "NISAR_L2_PR_GUNW_025_048_D_074_026_2000_SH_20260711T125814_20260711T125848_20260723T125813_20260723T125848_P05023_N_F_J_001.h5",
    "GUNW-011-TARGET4-LOW-SUPPORT": "NISAR_L2_PR_GUNW_025_098_A_016_026_4000_SH_20260714T233920_20260714T233954_20260726T233919_20260726T233954_P05023_N_F_J_001.h5",
    "GUNW-012-TARGET4-LOW-SUPPORT": "NISAR_L2_PR_GUNW_026_048_D_074_028_4000_SH_20260723T125813_20260723T125848_20260816T125812_20260816T125846_P05023_N_F_J_001.h5",
    "GUNW-013-PROVISIONAL-FOLLOWUP": "NISAR_L2_PR_GUNW_026_098_A_016_028_4000_SH_20260726T233919_20260726T233954_20260819T233918_20260819T233953_P05023_N_F_J_001.h5",
    "GOFF-023-BETA-control": "NISAR_L2_PR_GOFF_009_098_A_016_010_4000_SH_20260103T233921_20260103T233956_20260115T233922_20260115T233956_X05010_N_F_J_001.h5",
    "GOFF-025-PROVISIONAL": "NISAR_L2_PR_GOFF_024_098_A_016_025_4000_SH_20260702T233920_20260702T233955_20260714T233920_20260714T233954_P05023_F_F_J_001.h5",
    "GOFF-026-PROVISIONAL": "NISAR_L2_PR_GOFF_025_048_D_074_026_2000_SH_20260711T125814_20260711T125848_20260723T125813_20260723T125848_P05023_N_F_J_001.h5",
    "GOFF-027-PROVISIONAL": "NISAR_L2_PR_GOFF_025_098_A_016_026_4000_SH_20260714T233920_20260714T233954_20260726T233919_20260726T233954_P05023_N_F_J_001.h5",
    "GOFF-028-PROVISIONAL-MATCHED-DESCENDING": "NISAR_L2_PR_GOFF_026_048_D_074_028_4000_SH_20260723T125813_20260723T125848_20260816T125812_20260816T125846_P05023_N_F_J_001.h5",
    "GOFF-029-PROVISIONAL": "NISAR_L2_PR_GOFF_026_098_A_016_028_4000_SH_20260726T233919_20260726T233954_20260819T233918_20260819T233953_P05023_N_F_J_001.h5",
}
PROFILES = {
    "target4-9": list(RUN55_GRANULES),
    "provisional-gunw": ["GUNW-007-PROVISIONAL", "GUNW-009-PROVISIONAL"],
    "gunw-maturity-control": ["GUNW-006-BETA-control", "GUNW-007-PROVISIONAL", "GUNW-009-PROVISIONAL"],
    "target4-gunw-low-support": ["GUNW-009-PROVISIONAL", "GUNW-010-PROVISIONAL", "GUNW-011-TARGET4-LOW-SUPPORT", "GUNW-012-TARGET4-LOW-SUPPORT", "GUNW-013-PROVISIONAL-FOLLOWUP"],
    "target4-9-focused": ["GUNW-009-PROVISIONAL", "GUNW-010-PROVISIONAL", "GUNW-011-TARGET4-LOW-SUPPORT", "GUNW-012-TARGET4-LOW-SUPPORT", "GUNW-013-PROVISIONAL-FOLLOWUP", "GOFF-025-PROVISIONAL", "GOFF-026-PROVISIONAL", "GOFF-027-PROVISIONAL", "GOFF-028-PROVISIONAL-MATCHED-DESCENDING", "GOFF-029-PROVISIONAL"],
    "target4-9-matched": ["GUNW-006-BETA-control", "GUNW-009-PROVISIONAL", "GUNW-010-PROVISIONAL", "GUNW-011-TARGET4-LOW-SUPPORT", "GUNW-012-TARGET4-LOW-SUPPORT", "GUNW-013-PROVISIONAL-FOLLOWUP", "GOFF-023-BETA-control", "GOFF-025-PROVISIONAL", "GOFF-026-PROVISIONAL", "GOFF-027-PROVISIONAL", "GOFF-028-PROVISIONAL-MATCHED-DESCENDING", "GOFF-029-PROVISIONAL"],
    "goff-target4-9": ["GOFF-023-BETA-control", "GOFF-025-PROVISIONAL", "GOFF-026-PROVISIONAL", "GOFF-027-PROVISIONAL", "GOFF-028-PROVISIONAL-MATCHED-DESCENDING", "GOFF-029-PROVISIONAL"],
}
def file_sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""): h.update(block)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--profile",choices=sorted(PROFILES),default="target4-9-matched")
    ap.add_argument("--outdir",type=Path,default=Path("data/target49_raw"))
    ap.add_argument("--dry-run",action="store_true")
    ap.add_argument("--manifest",type=Path,default=None)
    args=ap.parse_args()
    selected=[(k,RUN55_GRANULES[k]) for k in PROFILES[args.profile]]
    print(f"profile={args.profile}; products={len(selected)}")
    for label,name in selected: print(f"{label}\t{name}")
    if args.dry_run:
        print("DRY_RUN: no network calls or downloads performed")
        print("Requires network + Earthdata Login; ASF docs recommend ~/.netrc or ASFSession.")
        return 0
    try:
        import asf_search as asf
    except ImportError:
        print("ERROR: install pinned requirements.txt (asf_search is required)",file=sys.stderr)
        return 2
    args.outdir.mkdir(parents=True,exist_ok=True)
    report={"schema_version":1,"created_utc":datetime.now(timezone.utc).isoformat(),
      "source":"Exact basenames in SARGuardian-backend Run #55 worker logs","profile":args.profile,
      "authentication":"Delegated to ASF/EDL; credentials are never written to this manifest.",
      "products":[],"errors":[]}
    for label,granule in selected:
        try:
            dest=args.outdir/granule
            if dest.is_file() and dest.stat().st_size:
                status="already_present"; actual=dest
            else:
                results=asf.granule_search([granule])
                if len(results)!=1:
                    raise RuntimeError(f"Exact basename lookup returned {len(results)} results; expected 1")
                print(f"DOWNLOADING {granule}",flush=True)
                token = os.environ.get("EARTHDATA_TOKEN", "").strip()
                if token:
                    session = asf.ASFSession().auth_with_token(token)
                    results.download(path=str(args.outdir), processes=1, session=session)
                else:
                    # ASF's preferred fallback is ~/.netrc. Auth is handled by
                    # the library; no credential is written to the manifest.
                    results.download(path=str(args.outdir), processes=1)
                matches=list(args.outdir.rglob(granule))
                if len(matches)!=1 or matches[0].stat().st_size==0:
                    raise RuntimeError(f"Exact non-empty basename not uniquely present; matches={len(matches)}")
                actual=matches[0]; status="downloaded"
            item={"label":label,"granule":granule,"path":str(actual),"status":status,
                  "bytes":actual.stat().st_size,"sha256":file_sha256(actual)}
            report["products"].append(item)
            print(f"{status.upper()} {granule} bytes={item['bytes']} sha256={item['sha256']}")
        except Exception as exc:
            report["errors"].append({"label":label,"granule":granule,"type":type(exc).__name__,"message":str(exc)})
            print(f"ERROR {label}: {type(exc).__name__}: {exc}",file=sys.stderr)
    target=args.manifest or (args.outdir/"download_manifest.json")
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(f"WROTE_MANIFEST {target}; present_or_downloaded={len(report['products'])}; errors={len(report['errors'])}")
    return int(bool(report["errors"]) or len(report["products"])!=len(selected))
if __name__=="__main__":
    raise SystemExit(main())
