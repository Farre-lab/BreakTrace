# BreakTrace

A Python pipeline for tracing Evolutionary Breakpoint Regions (EBRs) from pairwise synteny data back to their origin on a phylogenetic tree. Uses a reconstructed ancestral genome as the reference coordinate system.

Designed to work with output from [DESCHRAMBLER](https://github.com/deschrambler), but compatible with any tool that produces pairwise synteny blocks in the expected format.

[![GitHub](https://img.shields.io/badge/GitHub-Farre--lab%2FBreakTrace-blue)](https://github.com/Farre-lab/BreakTrace)

---

## Contents

- [Background](#background)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Format](#input-format)
- [breakpoints.py](#breakpointspy)
- [classify_ebrs.py](#classify_ebrspy)
- [Output](#output)
- [Test Data](#test-data)
- [Citation](#citation)

---

## Background

An EBR is the reference-genome interval between two consecutive syntenic blocks — the position where chromosomal order has changed between lineages. This tool:

1. **Finds breakpoints** between consecutive syntenic blocks for each species, classifying them as `break` (gap ≤ resolution threshold, confirmed rearrangement boundary) or `gap` (gap > resolution threshold, unresolved region).
2. **Merges overlapping breaks** from different species into consensus groups using a greedy interval merge, reporting both the union (widest span) and intersection (shared overlap) of each group.
3. **Classifies EBRs** using parsimony on a phylogenetic tree, relative to the named reconstructed ancestor.

Because the reference is a reconstructed ancestor rather than an extant species, a non-zero entry in the matrix means that species' genome **differs from the ancestor** at that position — i.e. a rearrangement occurred somewhere on the lineage from ancestor to that species.

### Classification classes

| Class | Condition | Interpretation |
|---|---|---|
| `ancestral_ancestor` | No ingroup signal; all outgroups positive | Break occurred in the ancestor itself |
| `ancestral` | All MRCA leaves positive | Single event in a named ancestral lineage |
| `partial_mrca` | ≥ threshold fraction of MRCA leaves positive | Ancestral event with some secondary losses |
| `convergent` | < threshold fraction of MRCA leaves positive | Independent events in unrelated lineages |
| `lineage_specific` | Single species positive | Event unique to one lineage |
| `no_signal` | No species positive | No rearrangement detected |

Each EBR is classified **twice**: conservatively (breaks only) and optimistically (breaks + gaps). An `ebr_uncertain` flag marks rows where the two differ.

---

## Installation

```bash
pip3 install pandas intervaltree biopython
```

Python 3.8 or later required.

---

## Quick Start

```bash
# Step 1 — find breakpoints and build the matrix
python3 breakpoints.py \
    --synteny-dir synteny/ \
    --resolution 300000 \
    --output results/matrix.tsv

# Step 2 — classify EBRs using the phylogenetic tree
python3 classify_ebrs.py \
    results/matrix.tsv \
    tree.nwk \
    --ancestor MAMMAL \
    --threshold 0.5 \
    --output results/matrix_classified.tsv
```

---

## Input Format

### Synteny files

One tab-separated file per species pair, no header row. Only columns 0–3 and 8 are used; the rest are accepted but ignored.

| Col | Index | Content | Example |
|---|---|---|---|
| ref_species | 0 | Reference (ancestor) name | mammalia |
| ref_chr | 1 | Reference chromosome | 1 |
| ref_start | 2 | Block start on reference | 4024542 |
| ref_end | 3 | Block end on reference | 7809330 |
| target_id | 4 | Target chr/scaffold (unused) | chr1 |
| target_start | 5 | Block start in target (unused) | 81825 |
| target_end | 6 | Block end in target (unused) | 245058708 |
| strand | 7 | Orientation (unused) | + |
| target_species | 8 | Target species name | hlallmiss2 |
| genome_type | 9 | Assembly type (unused) | chromosomes |

**Naming convention:** the species key is derived from the filename prefix before the first underscore. For example `hlallmiss2_mammalia.txt` → key `hlallmiss2`. This must match the target_species name in the file and the leaf name in the Newick tree. Species names are lowercased automatically.

### Newick tree

Standard Newick format with **named internal nodes**. Internal node names are used as classification labels and for defining the ingroup/outgroup split. The `--ancestor` argument must match an internal node name (case-insensitive).

```
((((human,chimp)ANC1,gorilla)ANC2,macaque)ANC3,mouse)root;
```

Branch lengths are accepted but ignored — only topology matters.

---

## breakpoints.py

Finds EBRs from pairwise synteny files and builds a cross-species matrix.

```
python3 breakpoints.py [-h] (--synteny-dir DIR | --input-files FILE [FILE ...])
                      --resolution INT [--increase INT] [--output FILE]
```

| Option | Default | Description |
|---|---|---|
| `--synteny-dir DIR` | — | Folder of `.txt` synteny files (mutually exclusive with `--input-files`) |
| `--input-files FILE…` | — | Explicit list of synteny files (mutually exclusive with `--synteny-dir`) |
| `--resolution INT` | 300000 | DESCHRAMBLER resolution in bp. Gaps wider than this are labelled `gap` and excluded from the greedy merge |
| `--increase INT` | 0 | bp padding added to each side of `break` intervals only |
| `--output FILE` | matrix.tsv | Output TSV path |

### Break vs gap

The raw gap between consecutive blocks is `nxt.ref_start - nxt.ref_end`:

- **`break`** — raw gap ≤ `--resolution` (or 0 for adjacent blocks). Interval is `[cur.ref_end, max(nxt.ref_start, cur.ref_end) + 1)`, expanded by `--increase` on each side.
- **`gap`** — raw gap > `--resolution`. Stored with raw coordinates, no expansion. Excluded from greedy merge. If a gap overlaps one or more break groups it is injected into those groups; otherwise it keeps its own row.

### Greedy merge

Break intervals from all species are sorted by `(ref_chr, ref_start)` and merged: a break joins the current group if it overlaps the running union interval. This is transitive — A overlaps B, B overlaps C → all three merge even if A and C don't directly overlap.

Each merged group reports:
- `union_start / union_end` — widest span (min start to max end)
- `inter_start / inter_end` — shared overlap (max start to min end); `NA` if no common region

---

## classify_ebrs.py

Classifies EBRs in the matrix using parsimony on a phylogenetic tree.

```
python3 classify_ebrs.py [-h] matrix_tsv newick_file
                        --ancestor NAME [--threshold FLOAT] [--output FILE]
```

| Option | Default | Description |
|---|---|---|
| `matrix_tsv` | — | TSV matrix from `breakpoints.py` |
| `newick_file` | — | Newick file with named internal nodes |
| `--ancestor NAME` | — | Reconstructed ancestor node name (case-insensitive) |
| `--threshold FLOAT` | 0.5 | Fraction of MRCA descendants required for `partial_mrca` vs `convergent` |
| `--output FILE` | matrix_classified.tsv | Output TSV path |

### Ingroup and outgroup

Derived automatically from `--ancestor`:

- **Ingroup** — leaves descending from the ancestor node. A break here means a rearrangement after the ancestor.
- **Outgroup** — all other leaves. A break here means the species differs from the ancestor's already-rearranged genome.

### Conservative vs optimistic classification

Each row is classified twice:

- **Conservative** (`ebr_class` / `ebr_label`) — only `break:` species count as positive.
- **Optimistic** (`ebr_class_with_gaps` / `ebr_label_with_gaps`) — `break:` and `gap:` species both count as positive.
- **`ebr_uncertain`** — `True` when the two classifications differ.

### Summary output files

When `--output` is provided, two summary TSVs are written alongside the main output:

- `*_summary_class.tsv` — EBR class counts, conservative vs with-gaps side by side
- `*_summary_label.tsv` — EBR label counts, conservative vs with-gaps side by side

Each has columns: `label`, `count_conservative`, `pct_conservative`, `count_with_gaps`, `pct_with_gaps`.

---

## Output

### matrix.tsv (from breakpoints.py)

| Column | Content |
|---|---|
| ref_species | Reference (ancestor) species name |
| ref_chr | Reference chromosome |
| union_start | Start of merged group (min across all members) |
| union_end | End of merged group (max across all members) |
| inter_start | Start of shared overlap region (NA if none) |
| inter_end | End of shared overlap region (NA if none) |
| \<species\> | `break:start--end` / `gap:start--end` / `0` (absent) |

### matrix_classified.tsv (from classify_ebrs.py)

All original matrix columns plus, inserted before species columns:

| Column | Content |
|---|---|
| ebr_class | Conservative class (breaks only) |
| ebr_label | Conservative label (node name / species name) |
| ebr_class_with_gaps | Optimistic class (breaks + gaps) |
| ebr_label_with_gaps | Optimistic label |
| ebr_uncertain | `True` if conservative and optimistic classifications differ |

---

## Test Data

A small synthetic test dataset is included in `test_data/` to verify the installation:

```bash
python3 breakpoints.py \
    --synteny-dir test_data/synteny/ \
    --resolution 100000 \
    --output test_data/matrix.tsv

python3 classify_ebrs.py \
    test_data/matrix.tsv \
    test_data/tree.nwk \
    --ancestor ANC2 \
    --output test_data/matrix_classified.tsv
```

The full test suite (requires `pytest`) covers all pipeline steps:

```bash
pip3 install pytest
pytest test_breakpoints.py test_classify_ebrs.py -v
```

---

## Citation

If you use this tool, please cite:

> [citation to be added on publication]

Code available at: https://github.com/Farre-lab/BreakTrace

---

## Licence

MIT Licence. See `LICENCE` for details.
