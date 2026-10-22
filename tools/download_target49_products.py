#!/usr/bin/env python3
"""Download the exact NISAR granules used by SARGuardian-backend Run #55.

Uses ASF's documented direct HTTPS product layout for a known exact granule
instead of catalogue exact-name search (which can omit older BETA revisions).
Authentication is delegated to ASFSession/.netrc. Credentials and signed URLs
are never logged or written to the manifest.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def direct_product_url(granule_name: str, product_type: str, maturity: str) -> str:
    """Build the official direct product URL for a fixed Run #55 granule name."""
    name = Path(granule_name).name
    if name != granule_name or not name.lower().endswith(".h5"):
        raise ValueError("Granule must be an exact HDF5 basename.")
    if product_type not in {"GUNW", "GOFF"} or maturity not in {"BETA", "PROVISIONAL"}:
        raise ValueError("Unsupported product type or maturity.")
    stem = name[:-3]
    collection = f"NISAR_L2_{product_type}_{maturity}_V1"
    return f"https://nisar.asf.earthdatacloud.nasa.gov/NISAR/{collection}/{stem}/{name}"


def authenticated_session(asf):
    token = os.environ.get("EARTHDATA_TOKEN", "").strip()
    username = os.environ.get("EARTHDATA_USERNAME", "").strip()
    password = os.environ.get("EARTHDATA_PASSWORD", "")
    if token:
        return asf.ASFSession().auth_with_token(token)
    if username and password:
        return asf.ASFSession().auth_with_creds(username, password)
    if username or password:
        raise RuntimeError("Set both Earthdata username and password, or provide EARTHDATA_TOKEN.")
    # ASFSession will use the host's supported .netrc credentials when present.
    return asf.ASFSession()


def download_exact(session, url: str, dest: Path, *, retries: int = 4) -> tuple[int, str]:
    """Stream one known product with timeouts, length checking and safe failure cleanup."""
    parsed = urlsplit(url)
    if parsed.scheme != "https" or parsed.hostname != "nisar.asf.earthdatacloud.nasa.gov":
        raise ValueError("Refusing a non-canonical NISAR product URL.")
    import requests
    last_status = None
    for attempt in range(1, retries + 1):
        response = None
        try:
            response = session.get(url, allow_redirects=True, stream=True, timeout=(30, 120))
            response.raise_for_status()
            header = response.headers.get("Content-Length")
            try:
                expected = int(header) if header else None
            except (TypeError, ValueError, OverflowError):
                expected = None
            dest.parent.mkdir(parents=True, exist_ok=True)
            downloaded = 0
            try:
                with dest.open("wb") as target:
                    for chunk in response.iter_content(chunk_size=8 * 1024 * 1024):
                        if chunk:
                            target.write(chunk)
                            downloaded += len(chunk)
            except Exception:
                dest.unlink(missing_ok=True)
                raise
            if downloaded <= 0:
                dest.unlink(missing_ok=True)
                raise RuntimeError("Empty product response.")
            if expected and expected != downloaded:
                dest.unlink(missing_ok=True)
                raise RuntimeError(f"Content-Length mismatch: expected {expected}, received {downloaded}.")
            return downloaded, sha256_file(dest)
        except Exception as exc:
            dest.unlink(missing_ok=True)
            response_for_status = getattr(exc, "response", None)
            last_status = getattr(response_for_status, "status_code", None)
            if attempt >= retries:
                break
            time.sleep(min(attempt * 5, 15))
        finally:
            if response is not None:
                response.close()
    raise RuntimeError(f"Exact direct download failed after {retries} attempts (HTTP status={last_status or 'unknown'}).")


def product_info(label: str) -> tuple[str, str]:
    product_type = label.split("-", 1)[0]
    maturity = "BETA" if "BETA-control" in label else "PROVISIONAL"
    return product_type, maturity


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--profile", choices=sorted(PROFILES), default="target4-9-matched")
    ap.add_argument("--label", action="append", choices=sorted(RUN55_GRANULES),
                    help="Download only this exact Run #55 product label; may be repeated. Overrides --profile.")
    ap.add_argument("--outdir", type=Path, default=Path("data/target49_raw"))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--manifest", type=Path, default=None)
    args = ap.parse_args()
    labels = args.label if args.label else PROFILES[args.profile]
    selected = [(label, RUN55_GRANULES[label]) for label in labels]
    print(f"profile={'explicit-labels' if args.label else args.profile}; products={len(selected)}")
    for label, name in selected:
        print(f"{label}\t{name}")
    if args.dry_run:
        print("DRY_RUN: no network calls or downloads performed")
        print("Live downloads use the canonical ASF direct-product URL and Earthdata credentials/.netrc.")
        return 0
    try:
        import asf_search as asf
    except ImportError:
        print("ERROR: install pinned requirements.txt (asf_search is required)", file=sys.stderr)
        return 2
    args.outdir.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": "Exact basenames from SARGuardian-backend Run #55 worker logs",
        "profile": "explicit-labels" if args.label else args.profile,
        "authentication": "Delegated to ASFSession/.netrc; credential values and signed URLs are never stored.",
        "products": [],
        "errors": [],
    }
    for label, granule in selected:
        product_type, maturity = product_info(label)
        dest = args.outdir / granule
        try:
            if dest.is_file() and dest.stat().st_size > 0:
                status = "already_present"
                nbytes = dest.stat().st_size
                digest = sha256_file(dest)
            else:
                url = direct_product_url(granule, product_type, maturity)
                session = authenticated_session(asf)
                try:
                    print(f"DOWNLOADING label={label}", flush=True)
                    nbytes, digest = download_exact(session, url, dest)
                finally:
                    session.close()
                status = "downloaded"
            item = {
                "label": label, "granule": granule, "product_type": product_type,
                "maturity": maturity, "path": str(dest), "status": status,
                "bytes": nbytes, "sha256": digest, "signed_url_stored": False,
            }
            report["products"].append(item)
            print(f"{status.upper()} label={label} bytes={nbytes} sha256={digest}")
        except Exception as exc:
            # Do not echo raw exceptions: HTTP clients can contain signed URLs.
            report["errors"].append({"label": label, "granule": granule,
                                     "type": type(exc).__name__, "message": "Direct exact-product download failed."})
            print(f"ERROR label={label} type={type(exc).__name__}", file=sys.stderr)
    manifest = args.manifest or (args.outdir / "download_manifest.json")
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"WROTE_MANIFEST {manifest}; products={len(report['products'])}; errors={len(report['errors'])}")
    return int(bool(report["errors"]) or len(report["products"]) != len(selected))


if __name__ == "__main__":
    raise SystemExit(main())
