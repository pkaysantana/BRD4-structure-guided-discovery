# BRD4-structure-guided-discovery

Reproducible experimental baseline for BRD4 bromodomain-1 (BD1) ligand
complexes, inspired by Vidler et al. (2013), doi:10.1021/jm4011302.

**Phase 1** analyses PDB entries **4MEN** (ligand 25K) and **4MEO**
(ligand 25V): ligand identity verification, physicochemical descriptors,
and heavy-atom protein–ligand contact geometry at 3.5 / 4.0 / 5.0 Å.
Scientific report: `docs/REPORT.md`.

## Quickstart

```bash
python -m venv .venv
.venv/Scripts/activate          # Windows; use .venv/bin/activate on Unix
pip install -r requirements.txt
pip install -e .
python -m brd4_analysis.cli all # download raw data, run full analysis
pytest                          # 14 tests
```

## Commands

| Command | Effect |
|---|---|
| `python -m brd4_analysis.cli fetch` | Download PDBs + CCD ligand SDFs + entry metadata; write `data/provenance/provenance.json` (URLs, UTC timestamps, SHA-256, PDB REVDAT info) |
| `python -m brd4_analysis.cli analyse` | Parse structures, verify ligand identity, compute descriptors, write contact tables + figures to `results/` |
| `python -m brd4_analysis.cli all` | `fetch` then `analyse` |

Optional `--data-dir` / `--results-dir` override the defaults (`data/`, `results/`).

## Layout

```
data/raw/          original PDB files + entry metadata JSON (committed, fingerprinted)
data/ccd/          ligand CCD ideal-coordinate SDFs (identity reference)
data/provenance/   provenance.json
src/brd4_analysis/ fetch, parse, descriptors, contacts, visualize, pipeline, cli
tests/             pytest suite (runs offline against committed data)
results/tables/    ligand_descriptors.csv, contacts_summary_*, contact_pairs_*, threshold_counts.csv
results/figures/   contact bar charts + ligand 2D depictions
docs/REPORT.md     Phase 1 scientific report
```

## Conventions

- Distances in Angstrom; contact defined as `min heavy-atom distance <= threshold`.
- Contacts are reported as geometry only — no interaction-type classification.
- Descriptors are computed on CCD reference chemistry; the PDB-extracted
  ligand is verified against it by bond-order template matching.
- Experimental measurements and predictions are kept strictly separate
  (Phase 1 contains no predictions).
