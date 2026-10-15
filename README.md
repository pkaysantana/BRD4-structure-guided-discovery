# BRD4-structure-guided-discovery

Reproducible experimental baseline for BRD4 bromodomain-1 (BD1) ligand
complexes, inspired by Vidler et al. (2013), doi:10.1021/jm4011302.

**Phase 1** analyses PDB entries **4MEN** (ligand 25K) and **4MEO**
(ligand 25V). It covers ligand/solvent separation, CCD graph-consistency
checks, physicochemical descriptors, and heavy-atom protein–ligand contact
geometry at 3.5 / 4.0 / 5.0 Å. The scientific report, including the
unresolved 4MEO ligand-identity question, is in `docs/REPORT.md`.

## Quickstart

```bash
python -m venv .venv
.venv/Scripts/activate          # Windows; use .venv/bin/activate on Unix
pip install -r requirements.txt
pip install -e .
python -m brd4_analysis.cli all # verify inputs, analyse, re-render report numbers
pytest                          # offline; uses committed data/
```

## Commands

| Command | Effect |
|---|---|
| `python -m brd4_analysis.cli fetch [--force]` | Download missing inputs (PDBs, CCD ligand SDFs, entry metadata). Existing files are verified against `data/provenance/provenance.json` and never silently overwritten. `--force` re-downloads everything |
| `python -m brd4_analysis.cli analyse` | Parse structures, check ligand graph consistency, compute descriptors, write contact tables, figures and `results/summary.json` |
| `python -m brd4_analysis.cli report` | Re-render the generated number blocks in `docs/REPORT.md` from `results/summary.json` |
| `python -m brd4_analysis.cli all` | `fetch`, `analyse`, then `report` |

Global options `--data-dir`, `--results-dir` and `--report` override the
defaults.

## Layout

```
data/raw/          original PDB files + entry metadata JSON (committed, byte-preserved)
data/ccd/          ligand CCD ideal-coordinate SDFs (chemistry reference)
data/provenance/   provenance.json: URL, retrieval time, SHA-256 per input (paths relative to data/)
src/brd4_analysis/ fetch, parse, descriptors, contacts, visualize, pipeline, report, cli
tests/             pytest suite (offline; exact regression + determinism tests)
results/summary.json      deterministic headline values (byte-reproducible)
results/run_metadata.json run timestamp/platform (varies per run; not a result)
results/tables/    ligand_descriptors, contacts_summary_*, contact_pairs_*, vdw_overlaps_*, threshold_counts
results/figures/   contact bar charts + ligand 2D depictions
docs/REPORT.md     Phase 1 scientific report
```

## Conventions

- Distances are in Angstrom. A contact means `min heavy-atom distance <= threshold`.
- Contacts are reported as geometry only, with no interaction-type
  classification and no energies.
- Descriptors are computed on the deposited CCD definition. The PDB ligand
  is checked for *graph consistency* with it. This check does not show
  protonation, tautomer or hydrogen positions.
- HBA/HBD counts follow RDKit's version-dependent `CalcNumHBA`/`CalcNumHBD`.
  The version-stable Lipinski N+O / NH+OH counts are reported alongside
  (see `docs/REPORT.md` §3).
- Experimental measurements and predictions are kept strictly separate.
  Phase 1 contains no predictions.
