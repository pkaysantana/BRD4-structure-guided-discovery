# Phase 4 adversarial screen v1

CPU-only adversarial validation of the frozen Phase 3 acceptance gates for
the BRD4-BD1–25K reconstruction (PDB 4MEN, Boltz-2.2.1 run
`experiments/phase2/4men_msa_s20261019_v1`).

## Contents

- `results.json` — full post-audit (`post_audit_v1`) per-case record:
  registered-gate verdicts, all metrics (contacts at 3.5/4.0/5.0 A,
  geometry submetrics, pose and CA RMSD), proposed candidate-check values,
  case parameters, seeds, environment versions, input hashes and Phase 4
  source hashes.
- `FINDINGS.md` — generated pass/fail matrix. Do not edit by hand.
- `manifest.json` — generated SHA-256 record for the two output artifacts
  plus input/source provenance.

## Reproduce

```bash
.venv/bin/python -m brd4_analysis.phase4.cli battery \
    --out experiments/phase4/adversarial_screen_v1
```

Deterministic: fixed-seed noise, fixed transform parameters and no
timestamps. Inputs are pinned by SHA-256 (config, prediction, committed
metrics) in `results.json → provenance`. Scientific fields are portable;
`results.json → environment` truthfully records the WSL runtime that
generated these published artifacts, so exact byte reproduction is asserted
for that recorded environment rather than every CI platform. The corrected
post-audit battery is explicitly identified as `post_audit_v1`; it is not a
backdated preregistration.

## Interpretation

These are synthetic, fixed controls probing what the registered gates can
and cannot detect. An "accepted synthetic variant" means the scoring
harness cannot distinguish the transformed coordinates under its thresholds.
That is not automatically a chemically invalid structure, an empirical model
false positive, an estimate of how often Boltz produces such structures, or
a population-level failure rate. See `docs/phase4/ADVERSARIAL_SCREEN.md` for
the analysis.
