# Agent notes

## Environment

- Python 3.12.10, Windows. Virtualenv at `.venv/` (interpreter
  `.venv/Scripts/python.exe`).
- Install: `pip install -r requirements.txt && pip install -e .`
- Pinned deps (all ≥7 days old at pin time, 2026-10-14): rdkit 2026.3.6,
  biopython 1.88, numpy 2.5.3, pandas 3.0.6, matplotlib 3.11.2,
  pytest 9.1.1, requests 2.34.2.
- No git identity is configured on this machine. Commit with
  `git -c user.name=... -c user.email=...` and do not edit git config.

## Verify

- Full pipeline: `.venv/Scripts/python.exe -m brd4_analysis.cli all`
  (offline when inputs verify; it re-renders the numbers in `docs/REPORT.md`).
- Tests: `.venv/Scripts/python.exe -m pytest` (offline; uses committed `data/`).
  - `test_regression.py`: exact pinned values (checksums, identities,
    descriptors, Tanimoto 6/72, key distances, contact counts) and
    report/summary consistency.
  - `test_determinism.py`: re-running the analysis reproduces the committed
    `summary.json` and CSVs byte-for-byte.
- After changing any computation, run `cli analyse` and then `cli report`,
  and commit `results/` together with `docs/REPORT.md`.

## Gotchas

- Never type numbers into `docs/REPORT.md` by hand. Put them in a
  `<!-- BEGIN GENERATED -->` block rendered by `report.py`.
- RCSB ligand SDFs (`ligands/download/{ID}_ideal.sdf`) carry explicit Hs.
  `load_ccd_mol` removes them. Keeping them changes the Morgan Tanimoto
  (0.154 with H versus 0.0833 heavy-atom); this caused a past report error.
- `Lipinski.NumHAcceptors` is the C++ `CalcNumHBA`, not the Python
  `Lipinski.HAcceptorSmarts` pattern. Its definition changed between RDKit
  2025.03 and 2025.09 (caffeine 6 -> 3, 25K 5 -> 4).
- 25V atom names: `O1` is the acetamide carbonyl O and `O` is the N-oxide O.
- The 4MEO publication text (compound 5, no N-oxide detected by MS)
  conflicts with the deposited 25V N-oxide model. Keep that uncertainty in
  every 4MEO statement.
- In Biopython, the water hetflag is `"W"`. Crystallisation additives (EDO,
  DMS, ...) are filtered by `SOLVENT_COMPONENT_IDS` in `constants.py`.
- Ligand component IDs: 4MEN → 25K, 4MEO → 25V. Do not treat EDO/DMS/HOH
  as ligands.
- `data/**` and `results/**` are `-text` in `.gitattributes`, and generated
  JSON is written with `write_bytes` so it has LF line endings.
