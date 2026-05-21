# BreakTrace — https://github.com/Farre-lab/BreakTrace
"""
breakpoints.py
==============
Define evolutionary breakpoint regions (EBRs) from pairwise synteny files
and build a cross-species overlap matrix.

Pipeline
--------
1. find_breakpoints()   – for each consecutive block pair, compute the gap
                          interval on the reference.  Intervals are labelled:
                            break  – raw gap <= resolution (or adjacent)
                            gap    – raw gap >  resolution (unresolved region)
                          increase is applied only to break intervals.
2. merge_breakpoints()  – greedy interval merge using break intervals only.
                          gap intervals are excluded from grouping but kept
                          in the output.
3. build_matrix()       – one row per merged group (breaks) plus one row per
                          gap interval.  Each species column contains:
                            "break:start--end"  for break intervals
                            "gap:start--end"    for gap intervals
                            0                   if absent

Input file format (tab-separated, no header)
--------------------------------------------
col 0  ref_species
col 1  ref_chr
col 2  ref_start      (int)
col 3  ref_end        (int)
col 4  target_id      (not used)
col 5  target_start   (not used)
col 6  target_end     (not used)
col 7  strand         (not used)
col 8  target_species
col 9  genome_type    (not used)

Usage
-----
    from breakpoints import run_pipeline

    matrix_df = run_pipeline(
        synteny_dir = "synteny/",
        resolution  = 300000,   # bp — gaps wider than this are excluded from merge
        increase    = 0,        # bp padding added to break intervals only
        output_tsv  = "matrix.tsv",
    )

Dependencies
------------
    pip install pandas
"""

from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional

import pandas as pd


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SyntenicBlock:
    ref_species   : str
    ref_chr       : str
    ref_start     : int
    ref_end       : int
    target_species: str


@dataclass
class Breakpoint:
    """
    An interval on the reference genome between two consecutive syntenic blocks.

    kind      : "break" if raw gap <= resolution (participates in merge)
                "gap"   if raw gap >  resolution (recorded but excluded from merge)
    ref_start : start of the interval (with increase applied for breaks)
    ref_end   : end   of the interval (with increase applied for breaks)
    raw_start : start of the raw gap (no increase)
    raw_end   : end   of the raw gap (no increase)
    """
    ref_species   : str
    ref_chr       : str
    target_species: str
    ref_start     : int
    ref_end       : int
    raw_start     : int
    raw_end       : int
    kind          : str   # "break" or "gap"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_synteny_file(path: str) -> List[SyntenicBlock]:
    """Parse one pairwise synteny file into a list of SyntenicBlock."""
    blocks: List[SyntenicBlock] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.lower().split("\t")
            cols = [c.strip() for c in cols]
            if len(cols) < 10:
                continue
            blocks.append(SyntenicBlock(
                ref_species    = cols[0],
                ref_chr        = cols[1],
                ref_start      = int(cols[2]),
                ref_end        = int(cols[3]),
                target_species = cols[8],
            ))
    return blocks


# ---------------------------------------------------------------------------
# Step 1 – find_breakpoints
# ---------------------------------------------------------------------------

def find_breakpoints(
    synteny_file: str,
    resolution  : int,
    increase    : int = 0,
) -> List[Breakpoint]:
    """
    Given one pairwise synteny file, return a list of Breakpoint objects
    representing the intervals between consecutive syntenic blocks.

    Gap size is measured as: nxt.ref_start - cur.ref_end  (raw, before increase).
    Adjacent blocks (gap size <= 0) are always treated as breaks.

    Labelling:
        gap_size <= resolution  →  "break": increase applied, used in merge
        gap_size >  resolution  →  "gap":   raw coordinates stored, excluded from merge

    Parameters
    ----------
    synteny_file : path to the input synteny TSV
    resolution   : DESCHRAMBLER resolution in bp; gaps wider than this → "gap"
    increase     : bp expansion applied to each side of break intervals only
    """
    blocks = _read_synteny_file(synteny_file)
    if not blocks:
        return []

    blocks.sort(key=lambda b: (b.ref_chr, b.ref_start, b.ref_end))

    breakpoints: List[Breakpoint] = []

    chrom_groups: Dict[str, List[SyntenicBlock]] = defaultdict(list)
    for b in blocks:
        chrom_groups[b.ref_chr].append(b)

    for chrom, chrom_blocks in chrom_groups.items():
        ref_sp = chrom_blocks[0].ref_species
        tar_sp = chrom_blocks[0].target_species

        for i in range(len(chrom_blocks) - 1):
            cur = chrom_blocks[i]
            nxt = chrom_blocks[i + 1]

            # Raw gap size (before any expansion)
            raw_gap  = nxt.ref_start - cur.ref_end   # <= 0 means adjacent
            raw_start = cur.ref_end
            raw_end   = max(nxt.ref_start, cur.ref_end) + 1

            if raw_gap > resolution:
                # Gap — store raw coordinates, no increase
                breakpoints.append(Breakpoint(
                    ref_species    = ref_sp,
                    ref_chr        = chrom,
                    target_species = tar_sp,
                    ref_start      = raw_start,
                    ref_end        = raw_end,
                    raw_start      = raw_start,
                    raw_end        = raw_end,
                    kind           = "gap",
                ))
            else:
                # Break — apply increase
                breakpoints.append(Breakpoint(
                    ref_species    = ref_sp,
                    ref_chr        = chrom,
                    target_species = tar_sp,
                    ref_start      = raw_start - increase,
                    ref_end        = raw_end   + increase,
                    raw_start      = raw_start,
                    raw_end        = raw_end,
                    kind           = "break",
                ))

    return breakpoints


# ---------------------------------------------------------------------------
# Step 2 – merge_breakpoints
# ---------------------------------------------------------------------------

def merge_breakpoints(
    all_breakpoints: List[Breakpoint],
) -> List[dict]:
    """
    Greedy interval merge using break intervals only.

    gap intervals are collected separately and returned as single-member
    groups (so they appear as rows in the matrix).

    Each group dict contains:
        ref_chr      – chromosome
        union_start  – min(ref_start) of break members
        union_end    – max(ref_end)   of break members
        inter_start  – max(ref_start) — start of shared overlap (pd.NA if none)
        inter_end    – min(ref_end)   — end   of shared overlap (pd.NA if none)
        kind         – "break" or "gap"
        members      – list of Breakpoint objects
    """
    if not all_breakpoints:
        return []

    breaks = [b for b in all_breakpoints if b.kind == "break"]
    gaps   = [b for b in all_breakpoints if b.kind == "gap"]

    # ── greedy merge of breaks ────────────────────────────────────────────
    sorted_breaks = sorted(breaks,
                           key=lambda b: (b.ref_chr, b.ref_start, b.ref_end))

    break_groups   : List[dict]       = []
    current_members: List[Breakpoint] = []
    current_union_s: int = 0
    current_union_e: int = 0
    current_chr    : str = ""

    def _flush():
        if not current_members:
            return
        inter_s = max(b.ref_start for b in current_members)
        inter_e = min(b.ref_end   for b in current_members)
        break_groups.append({
            "ref_chr"     : current_chr,
            "union_start" : current_union_s,
            "union_end"   : current_union_e,
            "inter_start" : inter_s if inter_s < inter_e else pd.NA,
            "inter_end"   : inter_e if inter_s < inter_e else pd.NA,
            "kind"        : "break",
            "members"     : list(current_members),
        })

    for bp in sorted_breaks:
        if bp.ref_chr != current_chr:
            _flush()
            current_members = [bp]
            current_union_s = bp.ref_start
            current_union_e = bp.ref_end
            current_chr     = bp.ref_chr
            continue

        if bp.ref_start <= current_union_e:
            current_members.append(bp)
            current_union_e = max(current_union_e, bp.ref_end)
        else:
            _flush()
            current_members = [bp]
            current_union_s = bp.ref_start
            current_union_e = bp.ref_end

    _flush()

    # ── gaps: collapse identical coordinates into one group ─────────────
    # Multiple species with the same gap coordinates → one row, not N rows.
    gap_index: Dict[tuple, dict] = {}
    for bp in gaps:
        key = (bp.ref_chr, bp.raw_start, bp.raw_end)
        if key not in gap_index:
            gap_index[key] = {
                "ref_chr"     : bp.ref_chr,
                "union_start" : bp.raw_start,
                "union_end"   : bp.raw_end,
                "inter_start" : pd.NA,
                "inter_end"   : pd.NA,
                "kind"        : "gap",
                "members"     : [],
            }
        gap_index[key]["members"].append(bp)
    gap_groups = list(gap_index.values())

    # ── inject gaps into overlapping break groups ────────────────────────
    # For each gap, find every break group whose union interval overlaps the
    # gap's raw coordinates.  If any overlap is found, add the gap as a
    # member of each overlapping break group and drop the gap's own row.
    # If no break group overlaps, keep the gap as its own row.

    # Index break groups by chromosome for fast lookup
    breaks_by_chr: Dict[str, List[dict]] = defaultdict(list)
    for g in break_groups:
        breaks_by_chr[g["ref_chr"]].append(g)

    orphan_gaps: List[dict] = []   # gaps that don't overlap any break group

    for gap_group in gap_groups:
        chrom     = gap_group["ref_chr"]
        gap_s     = gap_group["union_start"]
        gap_e     = gap_group["union_end"]
        gap_bp    = gap_group["members"][0]   # all members share same coords
        injected  = False

        for brk_group in breaks_by_chr.get(chrom, []):
            # Check overlap: gap and break group share any coordinate
            if gap_s <= brk_group["union_end"] and brk_group["union_start"] <= gap_e:
                # Inject all gap members into this break group
                for bp in gap_group["members"]:
                    brk_group["members"].append(bp)
                injected = True

        if not injected:
            orphan_gaps.append(gap_group)

    # ── merge and sort all groups ─────────────────────────────────────────
    all_groups = break_groups + orphan_gaps
    all_groups.sort(key=lambda g: (g["ref_chr"], g["union_start"]))
    return all_groups


# ---------------------------------------------------------------------------
# Step 3 – build_matrix
# ---------------------------------------------------------------------------

def build_matrix(
    merged_groups: List[dict],
    all_species  : List[str],
    ref_species  : str,
) -> pd.DataFrame:
    """
    Build the cross-species breakpoint matrix from merged groups.

    One row per merged group. Species columns contain:
        "break:start--end"   if the species has a break interval in this group
        "gap:start--end"     if the species has a gap interval in this group
        0                    if absent

    Columns:
        ref_species, ref_chr, union_start, union_end, inter_start, inter_end,
        <target_species_1>, ...
    """
    rows: List[dict] = []

    for group in merged_groups:
        row: dict = {
            "ref_species" : ref_species,
            "ref_chr"     : group["ref_chr"],
            "union_start" : group["union_start"],
            "union_end"   : group["union_end"],
            "inter_start" : group["inter_start"],
            "inter_end"   : group["inter_end"],
            **{sp: 0 for sp in all_species},
        }

        for bp in group["members"]:
            sp = bp.target_species
            if sp not in row:
                continue
            coord_str = f"{bp.kind}:{bp.raw_start}--{bp.raw_end}"
            if row[sp] == 0:
                row[sp] = coord_str
            else:
                existing = set(row[sp].split(","))
                existing.add(coord_str)
                row[sp] = ",".join(sorted(existing))

        rows.append(row)

    if not rows:
        cols = (["ref_species", "ref_chr",
                 "union_start", "union_end",
                 "inter_start", "inter_end"]
                + all_species)
        return pd.DataFrame(columns=cols)

    df = pd.DataFrame(rows)

    fixed_cols = ["ref_species", "ref_chr",
                  "union_start", "union_end",
                  "inter_start", "inter_end"]
    sp_cols = [s for s in all_species if s in df.columns]
    df = df[fixed_cols + sp_cols]

    return df.sort_values(["ref_chr", "union_start"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def _collect_input_files(
    synteny_dir: Optional[str],
    input_files: Optional[List[str]],
) -> List[str]:
    if synteny_dir is not None and input_files is not None:
        raise ValueError("Supply either synteny_dir or input_files, not both.")
    if synteny_dir is None and input_files is None:
        raise ValueError("Supply either synteny_dir (folder path) or input_files (list of paths).")
    if synteny_dir is not None:
        import glob as _glob
        pattern = os.path.join(synteny_dir, "*.txt")
        files   = sorted(_glob.glob(pattern))
        if not files:
            raise FileNotFoundError(
                f"No .txt files found in '{synteny_dir}'. "
                "Check the folder path and that files have a .txt extension."
            )
        print(f"Found {len(files)} .txt file(s) in '{synteny_dir}'")
        return files
    return list(input_files)


def run_pipeline(
    synteny_dir: Optional[str]       = None,
    input_files: Optional[List[str]] = None,
    resolution : int                 = 300_000,
    increase   : int                 = 0,
    output_tsv : Optional[str]       = None,
) -> pd.DataFrame:
    """
    Run the full breakpoint pipeline.

    Parameters
    ----------
    synteny_dir : folder containing pairwise synteny .txt files.
                  Supply this OR input_files, not both.
    input_files : explicit list of synteny file paths.
                  Supply this OR synteny_dir, not both.
    resolution  : DESCHRAMBLER resolution in bp (default 300000).
                  Gaps wider than this are labelled "gap" and excluded
                  from the greedy merge.
    increase    : bp padding added to each side of break intervals only.
    output_tsv  : if given, write the matrix DataFrame to this path.

    Returns
    -------
    pd.DataFrame – the breakpoint matrix
    """
    files = _collect_input_files(synteny_dir, input_files)

    all_breakpoints: List[Breakpoint] = []
    all_species    : List[str]        = []
    seen_sp        : set              = set()
    ref_species    : str              = ""

    for fpath in files:
        basename = os.path.basename(fpath)
        print(f"Finding breakpoints: {basename}")
        bps = find_breakpoints(fpath, resolution, increase)
        n_breaks = sum(1 for b in bps if b.kind == "break")
        n_gaps   = sum(1 for b in bps if b.kind == "gap")
        print(f"  → {n_breaks} breaks, {n_gaps} gaps")
        all_breakpoints.extend(bps)
        for bp in bps:
            if not ref_species:
                ref_species = bp.ref_species
            if bp.target_species not in seen_sp:
                all_species.append(bp.target_species)
                seen_sp.add(bp.target_species)

    n_breaks_total = sum(1 for b in all_breakpoints if b.kind == "break")
    n_gaps_total   = sum(1 for b in all_breakpoints if b.kind == "gap")
    print(f"Total: {n_breaks_total} breaks, {n_gaps_total} gaps "
          f"(resolution={resolution} bp)")

    print("Merging overlapping breaks …")
    merged = merge_breakpoints(all_breakpoints)
    n_break_groups = sum(1 for g in merged if g["kind"] == "break")
    n_gap_groups   = sum(1 for g in merged if g["kind"] == "gap")
    print(f"  → {n_break_groups} merged break groups, {n_gap_groups} gap rows")

    print("Building matrix …")
    matrix_df = build_matrix(merged, all_species, ref_species)
    print(f"Matrix shape: {matrix_df.shape}")

    if output_tsv:
        matrix_df.to_csv(output_tsv, sep="\t", index=False)
        print(f"Matrix written to {output_tsv}")

    return matrix_df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="BreakTrace: find and merge evolutionary breakpoint regions from pairwise synteny files.",
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--synteny-dir", dest="synteny_dir", metavar="DIR",
        help="Folder containing pairwise synteny .txt files.",
    )
    input_group.add_argument(
        "--input-files", dest="input_files", nargs="+", metavar="FILE",
        help="Explicit list of pairwise synteny TSV files.",
    )
    parser.add_argument(
        "--resolution", type=int, default=300_000,
        help="DESCHRAMBLER resolution in bp (default: 300000). "
             "Gaps wider than this are labelled 'gap' and excluded from merge.",
    )
    parser.add_argument(
        "--increase", type=int, default=0,
        help="bp padding added to break intervals only (default: 0).",
    )
    parser.add_argument(
        "--output", default="matrix.tsv",
        help="Output TSV file (default: matrix.tsv).",
    )
    args = parser.parse_args()

    run_pipeline(
        synteny_dir = args.synteny_dir,
        input_files = args.input_files,
        resolution  = args.resolution,
        increase    = args.increase,
        output_tsv  = args.output,
    )
