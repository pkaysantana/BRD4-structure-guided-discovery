# Experiment catalogue

All entries below are unexecuted. No GPU outputs or runtime manifests exist.

| Experiment ID | Role and status |
|---|---|
| `4men_single_sequence_s20261019` | Original preparation record from the first review. Superseded; files preserved unchanged. Its source hashes correspond to that earlier code version. |
| `4men_smoke_s20261019_v2` | Current optional no-MSA smoke test; prepared with an explicit `smoke_test` role. Not the primary scientific evaluation. |
| `4men_msa_s20261019_v1` | Primary full-MSA timing/structural evaluation, defined in `configs/phase2/4men_primary.json`. No experiment directory is created until a real MSA and provenance pass strict validation. |

The 32SC/9FWX comparison is an evidence audit of one duplicated public
complex, not two evaluation observations or two scheduled model runs.
Every future condition/run must have its own output directory. Existing
bundles are never silently resumed or overwritten.
