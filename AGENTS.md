# Agent notes

## Active environment (Phase 4)

- Live development environment: WSL Ubuntu. Worktree:
  `/home/don/projects/BRD4-structure-guided-discovery-phase4`, branch
  `phase-4-adversarial-validation`.
- Interpreter: `/home/don/projects/BRD4-structure-guided-discovery-phase4/.venv/bin/python`
  (CPython 3.12.14 via uv).
- The separate `/home/don/projects/BRD4-structure-guided-discovery` checkout is
  the frozen Phase 3 evidence checkout on `phase3-4men-msa-20261017`.
- The Windows checkout is not the Phase 4 authority. The Windows/Python
  3.12.10 instructions below are retained only for legacy Phase 1 context.

## Verification permitted for Phase 4

- CPU-only tests in this worktree: `.venv/bin/python -m pytest`.
- Phase 4 battery regeneration in this worktree:
  `.venv/bin/python -m brd4_analysis.phase4.cli battery --out
  experiments/phase4/adversarial_screen_v1`.
- Phase 4 regression files: `tests/test_phase4_perturb.py` and
  `tests/test_phase4_battery.py`.
- Read-only frozen-bundle verification may run `phase2.workflow.verify_run`
  against `experiments/phase2/4men_msa_s20261019_v1`; it must not modify
  Phase 2/3 source, data, results or evidence.
- Do not run Phase 1 `analyse`/`report` regeneration for Phase 4 work unless
  explicitly asked: those commands intentionally update `results/` and
  `docs/REPORT.md` and are not the Phase 4 verification path.
- The published Phase 4 artifacts record the WSL environment that generated
  them. Cross-platform tests compare deterministic scientific fields while
  retaining that provenance; byte-identical regeneration is required in the
  recorded environment.

## Legacy Phase 1 reference (Windows)

These instructions describe the earlier Windows checkout and are not the
active Phase 4 setup:

- Python 3.12.10, Windows. Virtualenv at `.venv/` (interpreter
  `.venv/Scripts/python.exe`).
- Install: `pip install -r requirements.txt && pip install -e .`
- Full pipeline: `.venv/Scripts/python.exe -m brd4_analysis.cli all`
  (offline when inputs verify; it re-renders `docs/REPORT.md`).
- Tests: `.venv/Scripts/python.exe -m pytest` (offline; uses committed
  `data/`).
- After changing a Phase 1 computation, run `cli analyse` and then
  `cli report`, and commit `results/` together with `docs/REPORT.md`.

## Common environment details

- Pinned deps (all >=7 days old at pin time, 2026-10-14): rdkit 2026.3.6,
  biopython 1.88, numpy 2.5.3, pandas 3.0.6, matplotlib 3.11.2,
  pytest 9.1.1, requests 2.34.2.
- No git identity is configured on this machine. Commit with
  `git -c user.name=... -c user.email=...` and do not edit git config.

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
- Ligand component IDs: 4MEN -> 25K, 4MEO -> 25V. Do not treat EDO/DMS/HOH
  as ligands.
- `data/**` and `results/**` are `-text` in `.gitattributes`, and generated
  JSON is written with byte writes so it has LF line endings.

## Phase 4 source isolation

- Phase 4 code lives only in `src/brd4_analysis/phase4/`. Do not add or
  modify files under `src/brd4_analysis/phase2/` or the hashed top-level
  modules: `phase2.workflow.analysis_source_hashes` is recomputed by
  `verify_run`, and a change there invalidates the frozen primary bundle.
- Phase 4 writes only below `experiments/phase4`. The CLI rejects traversal,
  output-root aliases and symlinked destinations/artifacts; on supported
  POSIX systems it pins the destination through a directory descriptor and
  atomically replaces artifact entries.
