# BreakTrace — https://github.com/Farre-lab/BreakTrace
"""
classify_ebrs.py
================
Classify Evolutionary Breakpoint Regions (EBRs) using a phylogenetic tree
and the parsimony principle, relative to a named reconstructed ancestor.

All species names (tree leaf names, matrix column names, ancestor argument)
are normalised to lowercase internally so capitalisation in the tree or
matrix never causes mismatches.  Internal node names used as classification
labels are preserved as written in the tree file.

Background
----------
The matrix produced by breakpoints.py uses a reconstructed ancestral genome
(e.g. from DESCHRAMBLER) as the reference coordinate system.  A non-zero
cell means that species' genome differs from the ancestor at that position.

Species are divided into two groups derived from the ancestor node:

    ingroup  – leaves descending FROM the named ancestor (descendants)
    outgroup – all other leaves (diverged BEFORE the ancestor)

Classification scheme
---------------------
Let S     = all species with a non-zero value
    S_in  = S ∩ ingroup
    S_out = S ∩ outgroup

    no_signal
        S is empty.

    lineage_specific  →  species name
        |S_in| == 1 and S_out is empty.  Single ingroup lineage event.
        |S_in| == 0 and |S_out| == 1.    Single outgroup lineage event.

    ancestral_<ancestor>  →  ancestor name
        |S_in| == 0 and S_out == all outgroups.
        Every outgroup differs from the ancestor → the rearrangement
        occurred in the ancestor itself before the outgroups diverged.

    ancestral  →  MRCA node name
    partial_mrca  →  MRCA node name + fraction
    convergent
        All other cases with |S| >= 2.
        MRCA is computed on S_in ∪ S_out (all positive species).
        fraction = |S| / |all leaves under MRCA(S)|
        fraction == 1.0              → ancestral
        fraction >= threshold        → partial_mrca
        fraction <  threshold        → convergent

Parameters
----------
matrix_tsv   : TSV from breakpoints.run_pipeline()
newick_file  : Newick file with named internal nodes
ancestor     : name of the reconstructed ancestor node (case-insensitive)
threshold    : float in (0,1].  Default 0.5.
output_tsv   : optional output path

Usage
-----
    from classify_ebrs import classify_ebrs

    df = classify_ebrs(
        matrix_tsv  = "matrix.tsv",
        newick_file = "primates.nwk",
        ancestor    = "MAMMAL",
        threshold   = 0.5,
        output_tsv  = "matrix_classified.tsv",
    )
"""

from __future__ import annotations

from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import pandas as pd
from Bio import Phylo


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------

def load_tree(newick_file: str) -> Phylo.BaseTree.Tree:
    """Parse a Newick file and return a Bio.Phylo tree."""
    tree = Phylo.read(newick_file, "newick")
    _normalise_leaf_names(tree)
    _validate_internal_names(tree)
    return tree


def _normalise_leaf_names(tree: Phylo.BaseTree.Tree) -> None:
    """Lowercase all leaf names in-place so comparisons are case-insensitive."""
    for terminal in tree.get_terminals():
        terminal.name = terminal.name.lower()


def _validate_internal_names(tree: Phylo.BaseTree.Tree) -> None:
    unnamed = [
        c for c in tree.find_clades()
        if not c.is_terminal() and not c.name
    ]
    if unnamed:
        print(
            f"  Warning: {len(unnamed)} internal node(s) have no name. "
            "These will be labelled 'unnamed_ancestor' in the output."
        )


def get_ingroup_outgroup(
    tree     : Phylo.BaseTree.Tree,
    ancestor : str,
) -> Tuple[FrozenSet[str], FrozenSet[str]]:
    """
    Derive ingroup and outgroup leaf sets from the named ancestor node.
    The ancestor name lookup is case-insensitive.

    ingroup  = leaves descending from the ancestor node (already lowercased)
    outgroup = all other leaves
    """
    # Search case-insensitively among internal nodes
    anc_clade = None
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.name and \
                clade.name.lower() == ancestor.lower():
            anc_clade = clade
            break

    if anc_clade is None:
        raise ValueError(
            f"Ancestor '{ancestor}' not found in tree. "
            "Check that the name matches an internal node label "
            "(comparison is case-insensitive)."
        )

    ingroup    = frozenset(t.name for t in anc_clade.get_terminals())
    all_leaves = frozenset(t.name for t in tree.get_terminals())
    outgroup   = all_leaves - ingroup
    return ingroup, outgroup


def get_leaves_under(clade: Phylo.BaseTree.Clade) -> FrozenSet[str]:
    """All leaf names under a clade (already lowercased)."""
    return frozenset(t.name for t in clade.get_terminals())


def find_mrca(
    tree    : Phylo.BaseTree.Tree,
    species : Set[str],
) -> Phylo.BaseTree.Clade:
    """Return the MRCA clade of the given species set."""
    return tree.common_ancestor(*species)


def node_label(clade: Phylo.BaseTree.Clade) -> str:
    """Return the node label as written in the tree (case preserved)."""
    return clade.name if clade.name else "unnamed_ancestor"


# ---------------------------------------------------------------------------
# MRCA-based classification
# ---------------------------------------------------------------------------

def _mrca_classify(
    positive  : Set[str],
    tree      : Phylo.BaseTree.Tree,
    threshold : float,
) -> Tuple[str, str]:
    """
    Run MRCA classification on an arbitrary set of positive species.
    fraction = |positive ∩ MRCA_leaves| / |MRCA_leaves|
    """
    mrca        = find_mrca(tree, positive)
    mrca_leaves = get_leaves_under(mrca)
    mrca_name   = node_label(mrca)

    fraction = len(positive & mrca_leaves) / len(mrca_leaves)

    if fraction == 1.0:
        return ("ancestral", mrca_name)
    if fraction >= threshold:
        return ("partial_mrca", f"{mrca_name} ({fraction:.2f})")
    # For convergent, list the positive species so the user can see
    # which lineages independently acquired the break.
    species_list = ",".join(sorted(positive))
    return ("convergent", f"convergent ({species_list})")


# ---------------------------------------------------------------------------
# Classification of a single EBR row
# ---------------------------------------------------------------------------

def _classify_one(
    positive_species : Set[str],
    tree             : Phylo.BaseTree.Tree,
    ingroup          : FrozenSet[str],
    outgroup         : FrozenSet[str],
    ancestor         : str,
    threshold        : float,
) -> Tuple[str, str]:
    """
    Classify one EBR row.  Returns (ebr_class, ebr_label).

    Decision table
    --------------
    S_in  S_out   action
    ----  -----   ------
    0     0       no_signal
    0     1       lineage_specific → outgroup species
    0     all     ancestral_ancestor
    0     2+      MRCA on S_out
    1     0       lineage_specific → ingroup species
    1+    any     MRCA on S_in ∪ S_out
    """
    s_in  = positive_species & ingroup
    s_out = positive_species & outgroup

    if not s_in and not s_out:
        return ("no_signal", "no_signal")

    if not s_in:
        if len(s_out) == 1:
            return ("lineage_specific", next(iter(s_out)))
        if s_out == outgroup:
            return ("ancestral_ancestor", f"ancestral_{ancestor}")
        return _mrca_classify(s_out, tree, threshold)

    if len(s_in) == 1 and not s_out:
        return ("lineage_specific", next(iter(s_in)))

    return _mrca_classify(s_in | s_out, tree, threshold)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def classify_ebrs(
    matrix_tsv  : str,
    newick_file : str,
    ancestor    : str,
    threshold   : float = 0.5,
    output_tsv  : Optional[str] = None,
) -> pd.DataFrame:
    """
    Classify EBRs in a breakpoint matrix using phylogenetic parsimony.

    Parameters
    ----------
    matrix_tsv  : path to the TSV produced by breakpoints.run_pipeline()
    newick_file : path to Newick file with named internal nodes
    ancestor    : name of the reconstructed ancestor node (case-insensitive)
    threshold   : partial_mrca threshold (default 0.5)
    output_tsv  : optional output path

    Returns
    -------
    pd.DataFrame with all original columns plus ebr_class and ebr_label
    """
    # ── load inputs ──────────────────────────────────────────────────────────
    print(f"Loading matrix:  {matrix_tsv}")
    df = pd.read_csv(matrix_tsv, sep="\t", dtype=str).fillna("0")

    # Lowercase species column names in the matrix
    META_COLS = {"ref_species", "ref_chr", "union_start", "union_end",
                 "inter_start", "inter_end", "ref_start", "ref_end", "decision"}
    df.columns = [
        c.lower() if c not in META_COLS else c
        for c in df.columns
    ]

    print(f"Loading tree:    {newick_file}")
    tree = load_tree(newick_file)   # leaf names lowercased inside load_tree

    ingroup, outgroup = get_ingroup_outgroup(tree, ancestor)
    print(f"Ancestor:        {ancestor}")
    print(f"Ingroup  ({len(ingroup):2d}):  {sorted(ingroup)}")
    print(f"Outgroup ({len(outgroup):2d}):  {sorted(outgroup)}")
    print(f"Threshold:       {threshold}")

    # ── identify usable species columns ──────────────────────────────────────
    species_cols = [c for c in df.columns if c.lower() not in
                    {m.lower() for m in META_COLS}]
    matrix_sp   = set(species_cols)
    tree_leaves = ingroup | outgroup

    not_in_tree   = matrix_sp - tree_leaves
    not_in_matrix = tree_leaves - matrix_sp
    if not_in_tree:
        print(f"  Warning: species in matrix but not in tree (ignored): "
              f"{sorted(not_in_tree)}")
    if not_in_matrix:
        print(f"  Warning: tree leaves absent from matrix: "
              f"{sorted(not_in_matrix)}")

    usable = matrix_sp & tree_leaves

    # ── classify each row ─────────────────────────────────────────────────────
    classes      : List[str]  = []
    labels       : List[str]  = []
    classes_gaps : List[str]  = []
    labels_gaps  : List[str]  = []
    uncertain    : List[bool] = []

    for _, row in df.iterrows():
        # Conservative: only confirmed break species
        positive_breaks = {
            col for col in usable
            if str(row[col]).strip() not in ("0", "", "nan")
            and str(row[col]).strip().startswith("break:")
        }

        # Optimistic: breaks + gaps (possible rearrangements)
        positive_with_gaps = {
            col for col in usable
            if str(row[col]).strip() not in ("0", "", "nan")
        }

        cls_brk, lbl_brk = _classify_one(
            positive_species=positive_breaks,
            tree=tree, ingroup=ingroup, outgroup=outgroup,
            ancestor=ancestor, threshold=threshold,
        )
        cls_gap, lbl_gap = _classify_one(
            positive_species=positive_with_gaps,
            tree=tree, ingroup=ingroup, outgroup=outgroup,
            ancestor=ancestor, threshold=threshold,
        )

        classes.append(cls_brk)
        labels.append(lbl_brk)
        classes_gaps.append(cls_gap)
        labels_gaps.append(lbl_gap)
        uncertain.append(cls_brk != cls_gap)

    df["ebr_class"]           = classes
    df["ebr_label"]           = labels
    df["ebr_class_with_gaps"] = classes_gaps
    df["ebr_label_with_gaps"] = labels_gaps
    df["ebr_uncertain"]       = uncertain

    # ── reorder: metadata + classification columns before species ─────────────
    meta_cols  = [c for c in ["ref_species", "ref_chr",
                               "union_start", "union_end",
                               "inter_start", "inter_end"]
                  if c in df.columns]
    class_cols = ["ebr_class", "ebr_label",
                  "ebr_class_with_gaps", "ebr_label_with_gaps",
                  "ebr_uncertain"]
    sp_cols    = [c for c in df.columns
                  if c not in set(meta_cols + class_cols)]
    df = df[meta_cols + class_cols + sp_cols]

    # ── summary ───────────────────────────────────────────────────────────────
    n_total     = len(df)
    n_uncertain = sum(uncertain)
    print("\nClassification summary (conservative — breaks only):")
    for cls, cnt in df["ebr_class"].value_counts().items():
        pct = 100 * cnt / max(n_total, 1)
        print(f"  {cls:<25} {cnt:>6}  ({pct:.1f}%)")
    print(f"\n  Rows where gaps change classification: "
          f"{n_uncertain} / {n_total} ({100*n_uncertain/max(n_total,1):.1f}%)")

    if output_tsv:
        df.to_csv(output_tsv, sep="\t", index=False)
        print(f"\nClassified matrix written to {output_tsv}")

        # ── summary tables written alongside the main output ──────────────────
        # Two files: one for class, one for label.
        # Each merges conservative and with-gaps counts side by side.
        stem = output_tsv.rsplit(".", 1)[0] if "." in output_tsv else output_tsv

        def _merged_summary(col_cons: str, col_gaps: str, path: str) -> None:
            cons = (
                df[col_cons].value_counts()
                .rename_axis("label").reset_index(name="count_conservative")
            )
            cons["pct_conservative"] = (
                cons["count_conservative"] / max(n_total, 1) * 100
            ).round(1)
            gaps = (
                df[col_gaps].value_counts()
                .rename_axis("label").reset_index(name="count_with_gaps")
            )
            gaps["pct_with_gaps"] = (
                gaps["count_with_gaps"] / max(n_total, 1) * 100
            ).round(1)
            merged = (
                cons.merge(gaps, on="label", how="outer")
                .fillna(0)
                .sort_values("count_conservative", ascending=False)
                .reset_index(drop=True)
            )
            merged["count_conservative"] = merged["count_conservative"].astype(int)
            merged["count_with_gaps"]    = merged["count_with_gaps"].astype(int)
            merged.to_csv(path, sep="\t", index=False)
            print(f"  {path}")

        print("\nSummary tables:")
        _merged_summary("ebr_class", "ebr_class_with_gaps", f"{stem}_summary_class.tsv")
        _merged_summary("ebr_label", "ebr_label_with_gaps", f"{stem}_summary_label.tsv")
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="BreakTrace: classify EBRs using parsimony on a phylogenetic tree."
    )
    parser.add_argument("matrix_tsv",  help="TSV matrix from breakpoints.py")
    parser.add_argument("newick_file", help="Newick file with named internal nodes")
    parser.add_argument(
        "--ancestor", required=True,
        help="Name of the reconstructed ancestor node in the tree (case-insensitive)"
    )
    parser.add_argument(
        "--threshold", type=float, default=0.5,
        help="Fraction threshold for partial_mrca vs convergent (default: 0.5)"
    )
    parser.add_argument(
        "--output", default="matrix_classified.tsv",
        help="Output TSV path (default: matrix_classified.tsv)"
    )
    args = parser.parse_args()

    classify_ebrs(
        matrix_tsv  = args.matrix_tsv,
        newick_file = args.newick_file,
        ancestor    = args.ancestor,
        threshold   = args.threshold,
        output_tsv  = args.output,
    )
