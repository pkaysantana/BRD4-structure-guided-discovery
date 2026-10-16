# Experiment catalogue

All entries below are unexecuted. No GPU outputs or runtime manifests exist.

| Experiment ID | Role and status |
|---|---|
| `4men_single_sequence_s20261019` | Original preparation record from the first review. Superseded; files preserved unchanged. Its source hashes correspond to that earlier code version. |
| `4men_smoke_s20261019_v2` | Preserved pre-review smoke preparation. Superseded after the terminal-OXT compatibility correction; its old analysis source hashes are intentional. |
| `4men_smoke_s20261019_v3` | Preserved atom-coverage correction checkpoint. Superseded by the final input/cache/provenance gate corrections. |
| `4men_smoke_s20261019_v4` | Current optional no-MSA preparation with strict informative-MSA validation and verified runtime provenance/cache gates. Not the primary scientific evaluation. |
| `4men_msa_s20261019_v1` | Primary full-MSA timing/structural evaluation, defined in `configs/phase2/4men_primary.json`. No experiment directory is created until a real MSA and provenance pass strict validation. |

The 32SC/9FWX comparison is an evidence audit of one duplicated public
complex, not two evaluation observations or two scheduled model runs.
Every future condition/run must have its own output directory. Existing
bundles are never silently resumed or overwritten.
