"""
timeseries.py
-------------
Assemble NISAR GUNW pairs into a cumulative LOS displacement time series.

Each interferogram measures a *difference*:

    phi(i,j)  =  d(t_j) - d(t_i)

Given a network of pairs, invert for d(t) at every acquisition epoch. This is
the small-baseline (SBAS) idea, reduced to what NISAR L2 actually needs -
the products are already unwrapped and geocoded, so all that is left is the
temporal inversion.

Three things this tool refuses to fudge
=======================================
1. It never mixes orbit geometries. Ascending path 98 and descending path 48
   measure different components of the same 3-D motion. They are inverted
   separately and reported separately.

2. It never bridges a disconnected network. If no chain of interferograms links
   two blocks of epochs, their displacements are not comparable and no amount
   of least squares will make them so. Each connected component gets its own
   series with its own zero.

3. It reports redundancy honestly. A component with as many pairs as unknowns
   has an exact solution and therefore no residual and no error estimate. That
   is not a good fit - it is no fit at all.

Usage
-----
    # Network structure straight from the catalogue, before downloading
    python src/timeseries.py --from-catalogue

    # Network structure of what you already have
    python src/timeseries.py --dir data/nisar_l2 --network

    # Region-mean inversion (fast, robust - start here)
    python src/timeseries.py --dir data/nisar_l2 --invert \
        --ref-lat 28.21 --ref-lon 85.47 --csv outputs/ts.csv --plot outputs/ts.png

Requires: numpy.  --dir modes also need h5py (via gunw_reader).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import re
import sys
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("timeseries")

from paths import NISAR, resolve  # noqa: E402

EVENT_DATE = date(2026, 8, 26)

GRANULE_HEAD_RE = re.compile(
    r"^NISAR_L2_(?P<processing>PR|UR)_(?P<product>GUNW|GOFF)_"
    r"(?P<cycle>\d+)_(?P<path>\d+)_(?P<direction>[AD])_"
    r"(?P<frame>\d+)_(?P<subswath>\d+)_(?P<bandwidth>\d+)_"
    r"(?P<polarisation>[A-Z]{2})_"
)
CRID_TAIL_RE = re.compile(
    r"_([A-Z]\d{5})_[A-Z]_[A-Z]_[A-Z]_\d{3}(?:\.h5)?$", re.IGNORECASE
)


def granule_metadata(name: str) -> dict | None:
    """Parse exact NISAR L2 footprint and processing/release identity."""
    basename = Path(str(name)).name
    match = GRANULE_HEAD_RE.match(basename)
    if not match:
        return None
    fields = match.groupdict()
    for key in ("cycle", "path", "frame", "subswath", "bandwidth"):
        fields[key] = int(fields[key])
    crid_match = CRID_TAIL_RE.search(basename)
    fields["crid"] = crid_match.group(1).upper() if crid_match else "UNKNOWN"
    return fields


# ---------------------------------------------------------------------------
# Pair inventory
# ---------------------------------------------------------------------------
class Pair:
    """One interferogram with explicit spatial footprint and processing identity."""

    __slots__ = (
        "path", "direction", "ref", "sec", "value", "n_px", "coherence", "source",
        "product", "cycle", "frame", "subswath", "bandwidth", "polarisation",
        "crid", "processing",
    )

    def __init__(self, path, direction, ref, sec, value=None, n_px=None,
                 coherence=None, source="", metadata=None):
        self.path = path
        self.direction = direction
        self.ref = ref
        self.sec = sec
        self.value = value
        self.n_px = n_px
        self.coherence = coherence
        self.source = source
        meta = metadata or granule_metadata(source) or {}
        for key in ("product", "cycle", "frame", "subswath", "bandwidth",
                    "polarisation", "crid", "processing"):
            setattr(self, key, meta.get(key))

    @property
    def span_days(self) -> int:
        return (self.sec - self.ref).days

    def network_key(self) -> tuple:
        """Separate frame/subswath/release series; never silently cross a CRID."""
        if self.frame is None:
            return (self.direction, self.path, None, None, None, None, None, self.product)
        return (self.direction, self.path, self.frame, self.subswath,
                self.bandwidth, self.polarisation, self.crid, self.product)

    def network_label(self, *, short: bool = False) -> str:
        label = f"{'ASCENDING' if self.direction == 'A' else 'DESCENDING'} path {self.path}"
        if self.frame is not None:
            label += (f" frame {self.frame} subswath {self.subswath} "
                      f"bandwidth {self.bandwidth} pol {self.polarisation}")
        if self.crid and self.crid != "UNKNOWN":
            label += f" CRID {self.crid}"
        if self.product:
            label += f" {self.product}"
        return label.replace("ASCENDING", "ASC").replace("DESCENDING", "DESC") if short else label

    def __repr__(self):
        return f"<{self.network_label()} {self.ref}->{self.sec}>"


def pairs_from_granule_names(names: list[str], paths: list, directions: list) -> list[Pair]:
    out = []
    for name, path, direction in zip(names, paths, directions):
        stamps = re.findall(r"_(\d{8})T\d{6}", name)
        if len(stamps) < 4:
            continue
        metadata = granule_metadata(name)
        if metadata is None:
            logger.warning("Skipping catalogue item without complete NISAR footprint identity: %s", name)
            continue
        out.append(Pair(
            path=metadata["path"],
            direction=metadata["direction"],
            ref=datetime.strptime(stamps[0], "%Y%m%d").date(),
            sec=datetime.strptime(stamps[2], "%Y%m%d").date(),
            source=name, metadata=metadata,
        ))
    return out

def pairs_from_catalogue() -> list[Pair]:
    """Query ASF for GUNW pairs over the AOI without downloading anything."""
    import nisar_acquisition as acq

    rows = acq.query_by_year("2025-01-01", date.today().isoformat(),
                             platform="NISAR", processingLevel="GUNW")
    return pairs_from_granule_names(
        [r["granuleName"] for r in rows],
        [r.get("path") for r in rows],
        [r.get("flightDirection") for r in rows],
    )


def pairs_from_stats(
    path: Path,
    fringe_mm: float = 121.95,
    layer: str | None = None,
    *,
    product: str | None = None,
    processing: str | None = None,
    crid: str | None = None,
    min_valid_pixels: int = 6,
) -> list[Pair]:
    """Build an epoch network from reader CSVs without silently mixing products.

    Duplicate PR/UR GUNW products may differ by whole phase cycles; those are
    resolved only within the same frame, subswath, bandwidth, polarisation,
    CRID and acquisition pair, using other intervals from that exact geometry.
    Different CRIDs are treated as separate processing releases. Conflicting
    GOFF PR/UR offsets are never resolved using a phase-fringe heuristic.
    """
    if min_valid_pixels < 1:
        raise ValueError("min_valid_pixels must be positive")
    if processing is not None:
        processing = processing.upper()
        if processing not in {"PR", "UR"}:
            raise ValueError("processing must be PR, UR or None")
    requested_product = product.upper() if product else None
    requested_crid = crid.upper() if crid else None

    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        logger.error("No rows in %s", path)
        return []

    if layer and "layer" in rows[0]:
        before = len(rows)
        rows = [r for r in rows if layer in (r.get("layer") or "")]
        logger.info("Layer filter '%s': %d of %d rows", layer, len(rows), before)
        if not rows:
            logger.error("No rows match layer '%s'", layer)
            return []

    filtered = []
    for row in rows:
        name = row.get("file", "")
        meta = granule_metadata(name)
        if meta is None:
            logger.warning("Cannot parse full NISAR product identity from %s; skipped",
                           name or "<blank>")
            continue
        row_product = meta["product"]
        row_processing = meta["processing"]
        row_crid = meta["crid"]

        # Reader metadata and the product basename must describe the same file.
        # A stale or hand-edited metadata column must not override the filename.
        supplied_processing = (row.get("processing") or "").strip().upper()
        supplied_crid = (row.get("crid") or "").strip().upper()
        if supplied_processing and supplied_processing not in {"NAN", row_processing}:
            raise ValueError(
                f"Stats metadata processing={supplied_processing} disagrees with filename "
                f"processing={row_processing}: {name}"
            )
        if supplied_crid and supplied_crid not in {"NAN", row_crid}:
            raise ValueError(
                f"Stats metadata CRID={supplied_crid} disagrees with filename CRID={row_crid}: {name}"
            )
        if requested_product and row_product != requested_product:
            continue
        if processing and row_processing != processing:
            continue
        if requested_crid and row_crid != requested_crid:
            continue
        row["_metadata"] = meta
        filtered.append(row)
    rows = filtered
    if not rows:
        logger.error("No rows remain after product/processing/CRID filters.")
        return []
    if processing:
        logger.info("Processing filter %s retained %d rows", processing, len(rows))
    if requested_crid:
        logger.info("CRID filter %s retained %d rows", requested_crid, len(rows))

    head = rows[0]
    if "median" in head:
        value_col, coh_col, wrapped = "median", "mean_coherence", True
        actual_product = "GUNW"
    elif "range_median_mm" in head:
        value_col, coh_col, wrapped = "range_median_mm", "mean_correlation", False
        actual_product = "GOFF"
    else:
        logger.error("%s has neither 'median' (GUNW) nor 'range_median_mm' (GOFF). "
                     "Columns present: %s", path, ", ".join(sorted(head)))
        return []
    if any(row["_metadata"]["product"] != actual_product for row in rows):
        raise ValueError("A stats CSV mixes GUNW and GOFF product schemas; create separate product CSVs.")
    if requested_product and requested_product != actual_product:
        logger.error("CSV values are %s, but --product %s was requested.",
                     actual_product, requested_product)
        return []
    logger.info("Reading %s values from column '%s'", actual_product, value_col)

    parsed = []
    invalid_rows = 0
    for row in rows:
        meta = row["_metadata"]
        name = row.get("file", "")
        try:
            stamps = re.findall(r"_(\d{8})T\d{6}", name)
            if len(stamps) < 4:
                raise ValueError("fewer than four acquisition timestamps")
            raw_value = (row.get(value_col) or "").strip()
            value = float(raw_value)
            raw_valid_px = (row.get("valid_px") or "").strip()
            if not raw_valid_px:
                raise ValueError("valid_px is missing; re-export using the current reader")
            valid_px = int(float(raw_valid_px))
            if not math.isfinite(value):
                raise ValueError(f"{value_col} is not finite")
            if valid_px < min_valid_pixels:
                raise ValueError(
                    f"only {valid_px} valid pixels; minimum is {min_valid_pixels}"
                )
            if "usable" in row and str(row["usable"]).strip().lower() not in {"true", "1", "yes"}:
                raise ValueError("row is marked unusable")
            a = datetime.strptime(stamps[0], "%Y%m%d").date()
            b = datetime.strptime(stamps[2], "%Y%m%d").date()
            if b <= a:
                raise ValueError(f"invalid acquisition interval {a} -> {b}")
            coherence_value = None
            raw_coh = (row.get(coh_col) or "").strip()
            if raw_coh:
                parsed_coherence = float(raw_coh)
                if math.isfinite(parsed_coherence):
                    coherence_value = parsed_coherence
            pair = Pair(
                path=meta["path"], direction=meta["direction"],
                ref=a, sec=b, value=value, n_px=valid_px,
                coherence=coherence_value, source=name, metadata=meta,
            )
            geometry = (
                meta["product"], meta["cycle"], meta["path"], meta["direction"],
                meta["frame"], meta["subswath"], meta["bandwidth"],
                meta["polarisation"],
            )
            parsed.append((meta["processing"], pair, geometry, meta["crid"]))
        except (TypeError, ValueError, OverflowError) as exc:
            invalid_rows += 1
            logger.warning("Stats row skipped (%s): %s", name or "<blank>", exc)
    if invalid_rows:
        logger.warning(
            "Skipped %d invalid/unusable stats row(s); no empty/NaN measurement was admitted.",
            invalid_rows,
        )
    if not parsed:
        logger.error("No finite stats rows with at least %d valid pixels.", min_valid_pixels)
        return []

    # One exact footprint/interval cannot silently mix CRIDs. The appropriate
    # action for a reprocessed product is an explicit --crid choice, not a
    # magnitude heuristic between different scientific releases.
    crids_by_pair: dict[tuple, set[str]] = defaultdict(set)
    for _proc, pair, geometry, pair_crid in parsed:
        crids_by_pair[(geometry, pair.ref, pair.sec)].add(pair_crid)
    conflicts = [
        (key, sorted(values)) for key, values in crids_by_pair.items()
        if len(values) > 1
    ]
    if conflicts:
        first_key, versions = conflicts[0]
        geometry, ref_date, sec_date = first_key
        raise ValueError(
            "Same product footprint and acquisition pair appears with multiple CRIDs "
            f"({', '.join(versions)}): {ref_date} -> {sec_date}, "
            f"frame={geometry[4]} subswath={geometry[5]} pol={geometry[7]}. "
            "Do not resolve a reprocessed-release conflict by magnitude; rerun with --crid."
        )

    groups: dict[tuple, list] = defaultdict(list)
    for proc, pair, geometry, pair_crid in parsed:
        groups[(geometry, pair.ref, pair.sec, pair_crid)].append((proc, pair))

    out, dropped = [], []
    for key, candidates in groups.items():
        if len(candidates) == 1:
            out.append(candidates[0][1])
            continue
        geometry, ref_date, sec_date, pair_crid = key
        processing_values = [proc for proc, _ in candidates]
        if len(set(processing_values)) != len(processing_values):
            raise ValueError(
                f"Duplicate {actual_product} rows for {ref_date}->{sec_date} have the same "
                f"processing code ({processing_values}); verify CSV uniqueness before inversion."
            )
        values = [pair.value for _, pair in candidates]
        same_value = all(
            math.isclose(values[0], value, rel_tol=1e-9, abs_tol=1e-6)
            for value in values[1:]
        )
        if same_value:
            # No scientific branch choice is required when values agree.
            best = next((candidate for candidate in candidates if candidate[0] == "PR"),
                        candidates[0])
        elif not wrapped:
            raise ValueError(
                f"Conflicting GOFF PR/UR offsets for {ref_date}->{sec_date}, "
                f"frame={geometry[4]} subswath={geometry[5]} CRID={pair_crid}. "
                "Offsets are not phase-wrapped; choose --processing PR or --processing UR explicitly."
            )
        else:
            raise ValueError(
                f"Conflicting GUNW PR/UR phase values for {ref_date}->{sec_date}, "
                f"frame={geometry[4]} subswath={geometry[5]} CRID={pair_crid}. "
                "A magnitude-only guess is not an adequate phase-ambiguity resolution; "
                "choose --processing PR or --processing UR explicitly, or resolve the ambiguity "
                "with a validated redundant-network closure analysis."
            )

        for proc, pair in candidates:
            if pair is best[1]:
                continue
            gap = abs(pair.value - best[1].value)
            if wrapped:
                logger.warning(
                    "%s -> %s frame %s sub %s CRID %s: dropped %s (%+.1f mm), "
                    "%.1f mm from kept %s = %.2f fringes",
                    ref_date, sec_date, geometry[4], geometry[5], pair_crid,
                    proc, pair.value, gap, best[0], gap / fringe_mm,
                )
            else:
                logger.warning(
                    "%s -> %s frame %s sub %s CRID %s: duplicate GOFF values agree; "
                    "kept %s, dropped %s",
                    ref_date, sec_date, geometry[4], geometry[5], pair_crid,
                    best[0], proc,
                )
            dropped.append((key, proc, pair.value, gap))
        out.append(best[1])

    by_interval = defaultdict(set)
    for pair in out:
        by_interval[(pair.direction, pair.path, pair.ref, pair.sec)].add(pair.network_key())
    if any(len(geometries) > 1 for geometries in by_interval.values()):
        logger.warning(
            "One or more acquisition intervals contain multiple frame/subswath/CRID "
            "footprints. They were retained as separate network identities rather than silently dropping all but one."
        )

    logger.info("Loaded %d %s pairs from %s (%d duplicate product(s) dropped)",
                len(out), actual_product, path.name, len(dropped))
    return out

def pairs_from_dir(directory: Path, product: str = "GUNW") -> list[Pair]:
    """
    Collect pairs from a folder tree. product is GUNW or GOFF.

    Both are interferometric pairs with the same naming, so the network logic
    is identical - only the reader differs. Searching recursively so the
    organised layout (GUNW/2025-11_winter/...) works without extra flags.
    """
    product = product.upper()
    files = sorted(p for p in directory.rglob("*.h5") if product in p.name.upper())
    if not files:
        logger.error("No *%s*.h5 under %s", product, directory)
        return []
    out = []
    for fp in files:
        stamps = re.findall(r"_(\d{8})T\d{6}", fp.name)
        metadata = granule_metadata(fp.name)
        if len(stamps) < 4 or metadata is None:
            logger.warning("Cannot parse full NISAR footprint/release identity from %s", fp.name)
            continue
        if metadata["product"] != product:
            logger.warning("Product filename does not match requested %s: %s", product, fp.name)
            continue
        out.append(
            Pair(
                path=metadata["path"],
                direction=metadata["direction"],
                ref=datetime.strptime(stamps[0], "%Y%m%d").date(),
                sec=datetime.strptime(stamps[2], "%Y%m%d").date(),
                source=str(fp),
                metadata=metadata,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Network structure
# ---------------------------------------------------------------------------
def connected_components(pairs: list[Pair]) -> list[list[date]]:
    """Group epochs into blocks linked by interferograms (union-find)."""
    parent: dict[date, date] = {}

    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for p in pairs:
        union(p.ref, p.sec)

    groups: dict[date, list[date]] = defaultdict(list)
    for epoch in parent:
        groups[find(epoch)].append(epoch)
    return sorted((sorted(v) for v in groups.values()), key=lambda g: g[0])


def describe_network(pairs: list[Pair]) -> dict:
    by_geom: dict[tuple, list[Pair]] = defaultdict(list)
    representative: dict[tuple, Pair] = {}
    for p in pairs:
        key = p.network_key()
        by_geom[key].append(p)
        representative.setdefault(key, p)

    print(f"\n{len(pairs)} interferograms across {len(by_geom)} spatial/release networks\n")
    summary = {}

    for key, group in sorted(by_geom.items(), key=lambda item: str(item[0])):
        comps = connected_components(group)
        epochs = sorted({e for p in group for e in (p.ref, p.sec)})
        label = representative[key].network_label()
        print("=" * 68)
        print(f"{label}   {len(group)} pairs, {len(epochs)} epochs, "
              f"{len(comps)} connected component(s)")
        print("=" * 68)

        if len(comps) > 1:
            print("  !! NETWORK IS DISCONNECTED - these blocks cannot be joined.")
            print("     Displacements in one block are NOT comparable to another.")

        comp_info = []
        for k, comp in enumerate(comps, 1):
            inside = [p for p in group if p.ref in comp and p.sec in comp]
            unknowns = len(comp) - 1
            redundancy = len(inside) - unknowns
            print(f"\n  Component {k}: {comp[0]} .. {comp[-1]}  "
                  f"({(comp[-1]-comp[0]).days} d, {len(comp)} epochs, {len(inside)} pairs)")
            print(f"    unknowns {unknowns}, redundancy {redundancy}", end="")
            if redundancy == 0:
                print("  <- exactly determined: no residual, no error estimate")
            elif redundancy < 0:
                print("  <- UNDER-DETERMINED")
            else:
                print()
            for p in sorted(inside, key=lambda x: x.sec):
                flag = "  <- last before event" if p.sec <= EVENT_DATE and \
                    p.sec == max(q.sec for q in inside if q.sec <= EVENT_DATE) else ""
                print(f"      {p.ref} -> {p.sec}  ({p.span_days:>2} d){flag}")
            comp_info.append({"epochs": [str(e) for e in comp], "n_pairs": len(inside),
                              "redundancy": redundancy})
        summary[label] = comp_info
        print()
    return summary


# ---------------------------------------------------------------------------
# Inversion
# ---------------------------------------------------------------------------
def invert_component(pairs: list[Pair], epochs: list[date]) -> dict:
    """
    Least-squares solve for cumulative displacement at each epoch.

    d(epochs[0]) is fixed at 0 - every series is relative to its own first
    acquisition, which is why disconnected components can never be compared.
    """
    epochs = sorted(set(epochs))
    index = {e: i for i, e in enumerate(epochs)}
    n_unknown = len(epochs) - 1
    if n_unknown < 1:
        return {"error": "single epoch"}

    # Only finite, measured edges can constrain a numeric inversion. Catalog
    # edges or pairs whose measurement failed must not artificially connect
    # components in the design matrix.
    usable = []
    for pair in pairs:
        if pair.value is None:
            continue
        try:
            value = float(pair.value)
        except (TypeError, ValueError, OverflowError):
            continue
        if not np.isfinite(value):
            continue
        if pair.ref not in index or pair.sec not in index or pair.sec <= pair.ref:
            return {"error": "pair dates are inconsistent with the supplied epoch set"}
        usable.append((pair, value))
    if not usable:
        return {"error": "no finite pairs carry a measured value"}

    G = np.zeros((len(usable), n_unknown))
    obs = np.zeros(len(usable))
    for r, (p, value) in enumerate(usable):
        i, j = index[p.ref], index[p.sec]
        if i > 0:
            G[r, i - 1] = -1.0
        if j > 0:
            G[r, j - 1] = +1.0
        obs[r] = value

    # numpy.linalg.lstsq returns a minimum-norm answer for a rank-deficient
    # graph instead of raising. Such an answer bridges unobserved epoch gaps
    # and must never be published as a cumulative displacement series.
    design_rank = int(np.linalg.matrix_rank(G))
    if design_rank < n_unknown:
        return {
            "error": f"rank-deficient measured network ({design_rank}/{n_unknown}); inversion withheld",
            "design_matrix_rank": design_rank,
            "unknown_epoch_count": n_unknown,
            "n_pairs": len(usable),
        }
    design_condition = float(np.linalg.cond(G))
    if not np.isfinite(design_condition) or design_condition > 1e10:
        return {
            "error": f"ill-conditioned measured network (condition={design_condition:.3g}); inversion withheld",
            "design_matrix_rank": design_rank,
            "design_matrix_condition": design_condition,
            "n_pairs": len(usable),
        }

    solution, *_ = np.linalg.lstsq(G, obs, rcond=None)
    cumulative = np.concatenate([[0.0], solution])
    residuals = G @ solution - obs

    dof = len(usable) - n_unknown
    if dof > 0:
        sigma = float(np.sqrt(np.sum(residuals ** 2) / dof))
        try:
            cov = sigma ** 2 * np.linalg.inv(G.T @ G)
            errors = np.concatenate([[0.0], np.sqrt(np.diag(cov))])
        except np.linalg.LinAlgError:
            errors = np.full(len(epochs), np.nan)
    else:
        sigma = float("nan")
        errors = np.full(len(epochs), np.nan)

    # linear velocity over the component
    days = np.array([(e - epochs[0]).days for e in epochs], dtype=float)
    velocity = velocity_err = float("nan")
    if len(epochs) >= 3:
        A = np.vstack([days, np.ones_like(days)]).T
        coef, *_ = np.linalg.lstsq(A, cumulative, rcond=None)
        velocity = float(coef[0])
        fit_res = A @ coef - cumulative
        if len(epochs) > 2:
            s = np.sqrt(np.sum(fit_res ** 2) / (len(epochs) - 2))
            try:
                velocity_err = float(s * np.sqrt(np.linalg.inv(A.T @ A)[0, 0]))
            except np.linalg.LinAlgError:
                pass

    return {
        "epochs": epochs,
        "cumulative_mm": cumulative,
        "error_mm": errors,
        "residual_rms_mm": float(np.sqrt(np.mean(residuals ** 2))) if len(residuals) else float("nan"),
        "sigma_mm": sigma,
        "dof": dof,
        "velocity_mm_per_day": velocity,
        "velocity_err_mm_per_day": velocity_err,
        "n_pairs": len(usable),
        "design_matrix_rank": design_rank,
        "design_matrix_condition": design_condition,
        "network_status": "REDUNDANT" if dof > 0 else "EXACTLY_DETERMINED_NO_RESIDUAL_DOF",
    }


def jackknife(label: str, comp_index: int, pairs: list[Pair],
              epochs: list[date]) -> dict | None:
    """
    Refit the linear velocity with each pair left out in turn.

    A network with no redundancy fits every observation exactly, so the formal
    error bars are undefined and a trend can look significant while resting
    entirely on one interferogram. Leave-one-out is the cheapest honest test
    available: if dropping a single pair changes the sign of the velocity, the
    trend is that pair, not the ground.

    This is not a refinement. On the current stack it is the step that
    separates a real signal from three artefacts.
    """
    usable = [p for p in pairs if p.value is not None]
    if len(usable) < 4:
        return None

    full = invert_component(usable, epochs)
    v0 = full.get("velocity_mm_per_day")
    if v0 is None or not np.isfinite(v0):
        return None

    print(f"\n  jackknife ({len(usable)} pairs, full fit {v0:+.3f} mm/day)")
    vs, flips, load_bearing = [], [], []
    for k, p in enumerate(usable):
        sub = usable[:k] + usable[k + 1:]
        # Removing a link from a chain splits it. The two halves then have
        # independent zeros and lstsq, handed a rank-deficient system, returns
        # a minimum-norm solution that silently bridges the break. Refitting
        # that number and calling it a jackknife would be worse than not
        # testing at all - it invents agreement or disagreement at random.
        if len(connected_components(sub)) > len(connected_components(usable)):
            load_bearing.append(p)
            print(f"    without {p.ref} -> {p.sec}   network splits - "
                  f"no trend can be fitted")
            continue
        eps = sorted({e for q in sub for e in (q.ref, q.sec)})
        r = invert_component(sub, eps)
        v = r.get("velocity_mm_per_day")
        if v is None or not np.isfinite(v):
            continue
        vs.append(v)
        flip = v * v0 < 0
        if flip:
            flips.append(p)
        print(f"    without {p.ref} -> {p.sec}   {v:+.3f} mm/day"
              + ("   <- SIGN FLIPS" if flip else ""))

    if load_bearing:
        print(f"\n    {len(load_bearing)} of {len(usable)} pairs are load-bearing: "
              f"removing any one\n    disconnects the network. This block is a "
              f"chain with redundancy 0, so\n    leave-one-out cannot test it - "
              f"there is nothing to leave out.")
        print(f"    The trend {v0:+.3f} mm/day is UNTESTED, not confirmed.")
        print(f"    To make it testable, add pairs that close loops: a 24-day\n"
              f"    interferogram spanning two existing 12-day steps gives\n"
              f"    redundancy 1 and a real residual. Those products already\n"
              f"    exist in the archive.")

    if not vs:
        return {"geometry": label, "component": comp_index,
                "velocity_mm_per_day": v0, "testable": False}
    robust = min(vs) * max(vs) > 0
    print(f"    range {min(vs):+.3f} to {max(vs):+.3f} mm/day   ", end="")
    if robust:
        print("ROBUST - every subset agrees in sign")
    else:
        who = ", ".join(f"{p.ref}->{p.sec}" for p in flips)
        print(f"NOT ROBUST\n    the sign depends on a single pair ({who}). "
              f"Treat the trend as\n    undetermined until another pair or the "
              f"other geometry supports it.")
    return {"geometry": label, "component": comp_index, "velocity_mm_per_day": v0,
            "jackknife_min": min(vs), "jackknife_max": max(vs),
            "robust": robust, "testable": True}


def measure_pairs(pairs: list[Pair], coh_threshold: float,
                  ref_lat: float | None, ref_lon: float | None,
                  clip_aoi: bool, flip_sign: bool,
                  target_lat: float | None = None, target_lon: float | None = None,
                  target_radius_px: int = 5, auto_ref: bool = False,
                  goff_layer: str = "layer2", common_ref: bool = False,
                  buffer_km: float = 2.0) -> None:
    """
    Fill Pair.value with a displacement measurement from each GUNW.

    Without a target, the median over the whole valid AOI is used. That is only
    meaningful if the AOI is dominated by the signal you care about - a small
    landslide inside a large stable AOI will be averaged into nothing.

    With --target-lat/--target-lon, the median over a small window at the target
    is used instead. Combined with a reference point on stable ground, that is
    the displacement OF the target RELATIVE TO the reference, which is the only
    thing InSAR can actually measure.
    """
    from gunw_reader import read_gunw

    lattices: dict[int, np.ndarray] = {}
    for p in pairs:
        if not p.source or not Path(p.source).exists():
            continue
        is_goff = "GOFF" in Path(p.source).name.upper()
        try:
            if is_goff:
                # GOFF has no phase. Slant-range offset IS the LOS component,
                # so it drops straight into the same inversion - the network
                # maths does not care which product measured the difference.
                from goff_reader import read_goff
                g = read_goff(Path(p.source), layer=goff_layer, clip_aoi=clip_aoi)
                key = next((k for k in g["layers"] if k.endswith(goff_layer)),
                           sorted(g["layers"])[0])
                L = g["layers"][key]
                r = {"valid": L["valid"],
                     "displacement_mm": L["range_m"] * 1000.0,
                     "coherence": L["correlation"],
                     "xs": L["xs"], "ys": L["ys"], "epsg": L["epsg"]}
            else:
                r = read_gunw(Path(p.source), coherence_threshold=coh_threshold,
                              ref_lat=ref_lat, ref_lon=ref_lon, auto_ref=auto_ref,
                              clip_aoi=clip_aoi, flip_sign=flip_sign)
        except Exception as exc:
            logger.error("%s: %s", Path(p.source).name, exc)
            continue
        valid = r["valid"]
        disp = r["displacement_mm"]

        if common_ref:
            # Onto the fixed AOI lattice, so every pair of a geometry shares
            # pixel centres and one datum can be set across all of them.
            import gunw_reader as _g
            grid = _g.aoi_grid(r["xs"], r["ys"], _g.AOI_RING, r["epsg"])
            if grid is None:
                logger.error("%s: off the AOI lattice, cannot join a common datum; "
                             "excluded", Path(p.source).name)
                continue
            lattices[id(p)] = (_g.place_on_grid(np.where(valid, disp, np.nan), grid),
                               grid, r["epsg"])
            p.n_px = int(valid.sum())
            p.coherence = (float(r["coherence"][valid].mean())
                           if r["coherence"] is not None and valid.any() else None)
            continue

        if target_lat is not None and target_lon is not None and r["xs"] is not None:
            xs, ys = r["xs"], r["ys"]
            tx, ty = target_lon, target_lat
            if r["epsg"] and r["epsg"] != 4326:
                try:
                    from pyproj import Transformer
                    tx, ty = Transformer.from_crs(4326, r["epsg"], always_xy=True).transform(
                        target_lon, target_lat)
                except ImportError:
                    logger.warning("pyproj missing - cannot place target on a projected grid.")
            j = int(np.argmin(np.abs(xs - tx)))
            i = int(np.argmin(np.abs(ys - ty)))
            t = target_radius_px
            sel = np.zeros_like(valid)
            sel[max(0, i - t): i + t + 1, max(0, j - t): j + t + 1] = True
            valid = valid & sel

        n = int(valid.sum())
        if n == 0:
            logger.warning("%s: no valid pixels in the measurement window, excluded",
                           Path(p.source).name)
            continue
        p.value = float(np.median(disp[valid]))
        p.n_px = n
        p.coherence = (float(r["coherence"][valid].mean())
                       if r["coherence"] is not None else None)
        # Include the exact granule basename. Candidate diagnostics must join
        # the target-window measurement to the AOI median from the same product;
        # dates alone collide across ascending/descending tracks and CRIDs.
        granule_identity = Path(p.source).name if p.source else "UNKNOWN"
        logger.info("%s -> %s [%s] : %+8.2f mm  (%d px, coh %.2f)",
                    p.ref, p.sec, granule_identity, p.value, n,
                    p.coherence if p.coherence is not None else float("nan"))

    if common_ref:
        apply_common_datum(pairs, lattices, target_lat, target_lon,
                           target_radius_px, buffer_km)


def apply_common_datum(pairs: list[Pair], lattices: dict, target_lat, target_lon,
                       radius_px: int, buffer_km: float) -> None:
    """
    Fill Pair.value on ONE datum per geometry, not one per pair.

    Every pair used to carry its own reference - an auto-reference block for
    GUNW, a deramp plane for GOFF - and the inversion summed those differences
    as though they shared a zero. Here the datum is the median over pixels
    valid in every pair of the geometry and more than buffer_km from the
    target, so it is never set on the ground being measured. The constant
    removed from each pair is logged: it is exactly the step that pair's
    private reference would have put into the series.
    """
    from rasterio.transform import rowcol
    from common_ref import buffer_mask, common_mask, rereference, window_median

    groups: dict[tuple, list[Pair]] = defaultdict(list)
    for p in pairs:
        if id(p) in lattices:
            groups[p.network_key()].append(p)

    for key, group in sorted(groups.items(), key=lambda item: str(item[0])):
        label = group[0].network_label()
        stack = [lattices[id(p)][0] for p in group]
        _, grid, epsg = lattices[id(group[0])]
        if any(lattices[id(p)][0].shape != stack[0].shape for p in group):
            logger.error("%s: pairs landed on different lattices - no common datum", label)
            continue
        centre = None
        if target_lat is not None and target_lon is not None:
            tx, ty = target_lon, target_lat
            if epsg and epsg != 4326:
                from pyproj import Transformer
                tx, ty = Transformer.from_crs(4326, epsg, always_xy=True).transform(
                    target_lon, target_lat)
            centre = rowcol(grid["transform"], tx, ty)
        stable = common_mask(stack)
        if centre is not None:
            px = grid["res_x"]
            stable &= ~buffer_mask(stable.shape, [centre], int(round(buffer_km * 1000 / px)))
        try:
            fixed, consts = rereference(stack, stable)
        except ValueError as exc:
            logger.error("%s: %s - no common datum, pairs left unmeasured",
                         label, exc)
            continue
        logger.info("%s: common datum on %d px valid in all %d pairs",
                    label, int(stable.sum()), len(group))
        for p, f, k in zip(group, fixed, consts):
            if centre is not None:
                p.value, n = window_median(f, centre, radius_px)
            else:
                v = f[np.isfinite(f)]
                p.value, n = (float(np.median(v)) if v.size else float("nan")), int(v.size)
            if not np.isfinite(p.value):
                p.value = None
                logger.warning("%s -> %s: no valid pixel at the target, excluded",
                               p.ref, p.sec)
                continue
            logger.info("%s -> %s : %+8.2f mm on the common datum (%d px), "
                        "datum error %+.1f mm = %+.2f mm/day",
                        p.ref, p.sec, p.value, n, k, k / max(p.span_days, 1))


def report_series(label: str, comp_index: int, result: dict) -> list[dict]:
    if "error" in result:
        print(f"\n{label} component {comp_index}: {result['error']}")
        return []

    print(f"\n{label}  -  component {comp_index}")
    print(f"  {result['n_pairs']} pairs, {len(result['epochs'])} epochs, dof {result['dof']}")
    if result["dof"] > 0:
        print(f"  residual RMS {result['residual_rms_mm']:.2f} mm  (sigma {result['sigma_mm']:.2f} mm)")
    else:
        print("  exactly determined - residual and error are undefined")

    print(f"\n  {'EPOCH':<12}{'CUMULATIVE mm':>15}{'+/- mm':>10}{'DAYS':>7}")
    print("  " + "-" * 44)
    t0 = result["epochs"][0]
    rows = []
    for e, d, err in zip(result["epochs"], result["cumulative_mm"], result["error_mm"]):
        errs = "   n/a" if not np.isfinite(err) else f"{err:>10.2f}"
        print(f"  {str(e):<12}{d:>15.2f}{errs}{(e-t0).days:>7}")
        rows.append({"geometry": label, "component": comp_index, "epoch": str(e),
                     "days_from_start": (e - t0).days,
                     "cumulative_mm": round(float(d), 3),
                     "error_mm": None if not np.isfinite(err) else round(float(err), 3)})

    v, ve = result["velocity_mm_per_day"], result["velocity_err_mm_per_day"]
    if np.isfinite(v):
        span = (result["epochs"][-1] - t0).days
        pm = f" +/- {ve:.3f}" if np.isfinite(ve) else ""
        print(f"\n  linear velocity {v:+.3f}{pm} mm/day   "
              f"({v*365:+.0f} mm/yr, over {span} d)")
        if np.isfinite(ve) and ve > 0 and abs(v) < 2 * ve:
            print("  NOT significant at 2 sigma - consistent with no motion.")
    return rows


# ---------------------------------------------------------------------------
def plot_series(all_rows: list[dict], out: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib required for --plot:  pip install matplotlib")
        return

    groups = defaultdict(list)
    for r in all_rows:
        groups[(r["geometry"], r["component"])].append(r)

    fig, ax = plt.subplots(figsize=(10, 6))
    for (geom, comp), rows in sorted(groups.items()):
        xs = [datetime.strptime(r["epoch"], "%Y-%m-%d") for r in rows]
        ys = [r["cumulative_mm"] for r in rows]
        es = [r["error_mm"] or 0 for r in rows]
        ax.errorbar(xs, ys, yerr=es, marker="o", capsize=3,
                    label=f"{geom} (block {comp})")
    ax.axvline(datetime(2026, 8, 26), color="crimson", ls="--", lw=1.2)
    ax.text(datetime(2026, 8, 26), ax.get_ylim()[1], " 26 Aug event",
            color="crimson", va="top", fontsize=9)
    ax.set_ylabel("Cumulative LOS displacement (mm)\n+ve away from satellite")
    ax.set_title("NISAR L-band cumulative displacement\n"
                 "each block has its own zero - blocks are not comparable",
                 fontsize=11)
    ax.grid(alpha=.3)
    ax.legend(fontsize=9)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    logger.info("Wrote %s", out)


def main() -> int:
    ap = argparse.ArgumentParser(description="GUNW pairs -> cumulative LOS time series")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--dir", metavar="DIR", nargs="?", const=str(NISAR),
                     help="folder of products (default: data/nisar_l2)")
    src.add_argument("--from-catalogue", action="store_true",
                     help="analyse network from ASF without downloading")
    src.add_argument("--from-stats", metavar="STATS.csv",
                     help="invert from a gunw_reader --csv summary, no products "
                          "needed (implies --invert)")
    ap.add_argument("--network", action="store_true", help="network structure only")
    ap.add_argument("--jackknife", action="store_true",
                    help="refit leaving each pair out; flags trends that rest "
                         "on a single interferogram")
    ap.add_argument("--invert", action="store_true", help="run the inversion (needs --dir)")
    ap.add_argument("--product", choices=("GUNW", "GOFF"), default="GUNW",
                    help="which product to build the series from")
    ap.add_argument("--processing", choices=("PR", "UR"),
                    help="with --from-stats, explicitly select routine (PR) or urgent (UR) products")
    ap.add_argument("--crid", metavar="CRID",
                    help="with --from-stats, select one exact release (e.g. X05010 or P05023); prevents mixed-release duplicate selection")
    ap.add_argument("--min-valid-pixels", type=int, default=6,
                    help="with --from-stats, skip rows with fewer usable pixels (default: 6)")
    ap.add_argument("--goff-layer", default="layer2",
                    help="GOFF correlation-window layer (layer3 is usually quietest)")
    ap.add_argument("--aoi", choices=("source", "langtang", "lhende"), default="source",
                    help="analysis AOI. Defaults to the confirmed detachment zone; it used to default to the wider langtang box")
    ap.add_argument("--auto-ref", action="store_true",
                    help="pick the reference automatically (recommended)")
    ap.add_argument("--coh-threshold", type=float, default=0.3)
    ap.add_argument("--ref-lat", type=float)
    ap.add_argument("--ref-lon", type=float)
    ap.add_argument("--target-lat", type=float, help="measure here instead of the AOI median")
    ap.add_argument("--target-lon", type=float)
    ap.add_argument("--target-radius", type=int, default=5, help="target window half-width in px")
    ap.add_argument("--common-ref", action="store_true",
                    help="reference every pair of a geometry to ONE datum: pixels "
                         "valid in all of them, --buffer-km clear of the target. "
                         "Without it each pair keeps its own reference and the "
                         "series inherits a step at every epoch")
    ap.add_argument("--buffer-km", type=float, default=2.0,
                    help="with --common-ref, no datum pixel within this of the target")
    ap.add_argument("--no-clip", action="store_true")
    ap.add_argument("--flip-sign", action="store_true")
    ap.add_argument("--csv", metavar="OUT.csv")
    ap.add_argument("--plot", metavar="OUT.png")
    ap.add_argument("--json", metavar="OUT.json", help="dump the network summary")
    args = ap.parse_args()

    from gunw_reader import AOIS, point_in_ring, set_aoi
    set_aoi(args.aoi)

    if (args.target_lat is None) != (args.target_lon is None):
        ap.error("--target-lat and --target-lon must be supplied together")
    if args.target_lat is not None:
        if not (np.isfinite(args.target_lat) and np.isfinite(args.target_lon)
                and -90 <= args.target_lat <= 90 and -180 <= args.target_lon <= 180):
            ap.error("target coordinates must be finite WGS84 latitude/longitude")
        if not point_in_ring(args.target_lon, args.target_lat, AOIS[args.aoi]):
            ap.error(f"target point ({args.target_lat:.6f}, {args.target_lon:.6f}) is outside the exact {args.aoi} AOI polygon; select the AOI that actually contains the target")

    if args.target_radius < 0:
        ap.error("--target-radius must be non-negative")
    if args.buffer_km < 0:
        ap.error("--buffer-km must be non-negative")

    if (args.processing or args.crid or args.min_valid_pixels != 6) and not args.from_stats:
        ap.error("--processing, --crid and --min-valid-pixels are supported only with --from-stats")

    if args.from_stats:
        if args.min_valid_pixels < 1:
            ap.error("--min-valid-pixels must be positive")
        try:
            pairs = pairs_from_stats(
                resolve(args.from_stats),
                layer=args.goff_layer if args.product == "GOFF" else None,
                product=args.product,
                processing=args.processing,
                crid=args.crid,
                min_valid_pixels=args.min_valid_pixels,
            )
        except ValueError as exc:
            ap.error(str(exc))
        args.invert = True          # the values are already in the CSV
    elif args.from_catalogue:
        pairs = pairs_from_catalogue()
    else:
        pairs = pairs_from_dir(resolve(args.dir), args.product)
    if not pairs:
        logger.error("No pairs found.")
        return 1

    summary = describe_network(pairs)
    if args.json:
        Path(args.json).write_text(json.dumps(summary, indent=1))
        logger.info("Wrote %s", args.json)

    if args.network or not args.invert:
        if args.from_catalogue:
            print("Nothing downloaded yet - run nisar_acquisition.py --download GUNW "
                  "then re-run with --dir to invert.")
        return 0

    if args.from_catalogue:
        logger.error("--invert needs local files. Use --dir.")
        return 1

    if not args.from_stats:
        if args.ref_lat is None or args.ref_lon is None:
            logger.warning(
                "No --ref-lat/--ref-lon. Each interferogram keeps its own arbitrary "
                "offset, so the assembled series will be meaningless. Strongly advised."
            )

        logger.info("Measuring displacement for %d pairs...", len(pairs))
        if args.target_lat is None:
            logger.warning(
                "No --target-lat/--target-lon: using the median over the whole AOI. "
                "A localised signal will be averaged away. Give a target for hazard work."
            )
        measure_pairs(pairs, args.coh_threshold, args.ref_lat, args.ref_lon,
                      not args.no_clip, args.flip_sign,
                      args.target_lat, args.target_lon, args.target_radius,
                      auto_ref=args.auto_ref, goff_layer=args.goff_layer,
                      common_ref=args.common_ref, buffer_km=args.buffer_km)

    by_geom = defaultdict(list)
    for p in pairs:
        by_geom[p.network_key()].append(p)

    all_rows: list[dict] = []
    representatives = {key: group[0] for key, group in by_geom.items()}
    for key, group in sorted(by_geom.items(), key=lambda item: str(item[0])):
        label = representatives[key].network_label(short=True)
        for k, comp in enumerate(connected_components(group), 1):
            inside = [p for p in group if p.ref in comp and p.sec in comp]
            all_rows += report_series(label, k, invert_component(inside, comp))
            if args.jackknife:
                jackknife(label, k, inside, comp)

    if args.csv and all_rows:
        keys = ["geometry", "component", "epoch", "days_from_start", "cumulative_mm", "error_mm"]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(all_rows)
        logger.info("Wrote %s (%d rows)", args.csv, len(all_rows))

    if args.plot and all_rows:
        plot_series(all_rows, Path(args.plot))

    return 0


if __name__ == "__main__":
    sys.exit(main())
