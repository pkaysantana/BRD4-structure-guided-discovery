# Experiment catalogue

Status updated 17 October 2026. One experiment, the primary full-MSA
evaluation `4men_msa_s20261019_v1`, has been executed on a GPU and
analysed. Every other bundle below is a preparation record that has never
been executed. Prepared files in every bundle are preserved unchanged;
their own `status` fields record the state at preparation time and are
not rewritten after execution.

| Experiment ID | Role | Status |
|---|---|---|
| `4men_single_sequence_s20261019` | Original preparation record from the first review | Superseded; never executed. Its source hashes correspond to that earlier code version. |
| `4men_smoke_s20261019_v2` | Pre-review no-MSA smoke preparation | Superseded after the terminal-OXT compatibility correction; never executed. Its old analysis source hashes are intentional. |
| `4men_smoke_s20261019_v3` | Atom-coverage correction checkpoint | Superseded by the final input/cache/provenance gate corrections; never executed. |
| `4men_smoke_s20261019_v4` | Optional no-MSA smoke preparation with strict informative-MSA validation and verified runtime provenance/cache gates | Prepared, never executed. Not the primary scientific evaluation and not needed to interpret it. |
| `4men_msa_s20261019_v1` | Primary full-MSA structural evaluation (Arm A), defined in [`configs/phase2/4men_primary.json`](../../configs/phase2/4men_primary.json) | **Executed and analysed.** One seed, one sample. All preregistered gates passed. See the [post-run report](../../docs/phase2/4MEN_PRIMARY_POST_RUN.md). |

## Primary experiment record (`4men_msa_s20261019_v1`)

The stages are kept separate. A later stage does not rewrite the files of
an earlier one.

| Stage | Evidence | State |
|---|---|---|
| Preregistration | `config.json`, `preregistration.json`, `command.json`, `validation.json`, `manifest.json`, inputs and maps; committed in `c9b09dc` before inference | Frozen. Fields such as `"executed": false` and `"status": "awaiting_user_review_no_prediction_generated"` describe the bundle as prepared, not its current state. |
| Runtime provenance | `runtime.json` with its `runtime.sha256` receipt, captured before inference | Complete; receipt matches. |
| Inference | `run_observation.json`, `logs/boltz_predict.log`, `outputs/boltz_results_4men/` (prediction, confidence, PAE/PDE/pLDDT arrays, processed inputs) | Boltz 2.2.1 on one RTX 5070 Ti Laptop GPU; exit code 0; 89.03 s wall time. |
| Scientific evaluation | `analysis/metrics.json`, written once by `brd4_analysis.phase2.cli analyse` | Harness success. Ligand RMSD 0.140 Å and protein Cα RMSD 0.191 Å after the protein fit; 4.0 Å contact precision and recall 1.0; geometry screen passed. |
| Independent numerical audit | Two externally generated reviews (Codex and a separate Claude Code session). No transcript or output file from either is present in this repository. | The reviews report that recomputation agreed with the stored metrics. They also report a 12% centroid-expansion test, which is planned as a Phase 4 unit test. Their limitations are recorded in the post-run report. |
| Execution provenance | Retrospective statement in the post-run report, section D | Run on the project owner's interactive instruction. No contemporaneous signed approval record exists, and none has been created or backdated. |
| Outstanding limitations | [Post-run report](../../docs/phase2/4MEN_PRIMARY_POST_RUN.md), sections G and H | Single case, single seed, unresolved training-set exposure, 3.5 Å contact sensitivity, permissive geometry screen. |

The result shows that the pipeline works end to end on one known,
pre-cutoff complex. It is not evidence of generalisation, of prospective
pose-prediction accuracy or of affinity.

## Future work

Phase 4 (robustness and memorisation probes) is in offline development
on a separate branch. No further GPU inference has been approved, and no
Phase 4 experiment is recorded in this catalogue. Phase 4 must not modify the files of `4men_msa_s20261019_v1`.

## Rules

The 32SC/9FWX comparison is an evidence audit of one duplicated public
complex, not two evaluation observations or two scheduled model runs.
Every future condition/run must have its own output directory and
experiment ID. Existing bundles are never silently resumed or
overwritten. A changed MSA, seed set, model or threshold requires a new
experiment ID.
