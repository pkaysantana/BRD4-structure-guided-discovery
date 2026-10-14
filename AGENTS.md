# Agent notes

## Environment

- Python 3.12.10, Windows. Virtualenv at `.venv/` (interpreter
  `.venv/Scripts/python.exe`).
- Install: `pip install -r requirements.txt && pip install -e .`
- Pinned deps (all ≥7 days old at pin time, 2026-10-14): rdkit 2026.3.6,
  biopython 1.88, numpy 2.5.3, pandas 3.0.6, matplotlib 3.11.2,
  pytest 9.1.1, requests 2.34.2.

## Verify

- Full pipeline: `.venv/Scripts/python.exe -m brd4_analysis.cli all`
- Tests: `.venv/Scripts/python.exe -m pytest` (14 tests, offline — they
  run against committed files in `data/raw/`).

## Gotchas

- RCSB ligand SDFs (`ligands/download/{ID}_ideal.sdf`) carry explicit Hs;
  `load_ccd_mol` removes them so the CCD molecule templates the
  heavy-atom-only PDB HETATM records.
- RDKit's Lipinski `NumHAcceptors` uses a SMARTS definition — caffeine
  gives 3, not the naive N+O count of 6. Test encodes this deliberately.
- In Biopython, water hetflag is `"W"`; crystallisation additives (EDO,
  DMS, ...) are filtered by `SOLVENT_COMPONENT_IDS` in `constants.py`.
- Ligand component IDs: 4MEN → 25K, 4MEO → 25V. Do not treat EDO/DMS/HOH
  as ligands.
