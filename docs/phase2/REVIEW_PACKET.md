# Phase 2 review packet

Prepared on 15 October 2026 for the Lilly AI for Chemistry interview on
Monday 19 October 2026 at 11:30 am. **No GPU job or Boltz prediction has
run; no model weights have been installed or downloaded.**

## A. Phase 1 merge recommendation

**Hold PR #1 at `5c06db2` for the CSV portability fix.** The existing 45
tests pass on Windows. The current head still lets pandas choose the OS
line ending, while every committed result table is CRLF. A prepared
[patch](phase1_csv_fix.patch) makes that explicit without changing any
existing result bytes. Apply it to the Phase 1 branch and verify on
Windows/Linux before merging. There are no GitHub status checks reported.

The deposited-model scientific scope is defensible. The 4MEO identity
uncertainty remains unresolved and need not block that baseline if it is
retained. See [the merge review](PHASE1_REVIEW.md) and
[the corrected research integration](../research/REVIEWED_INTEGRATION.md).
PR #1 was neither modified nor merged; Phase 1 source/results/report and
the supplied research originals were preserved.

Four of six research files were available and read. The Boltz protocol
and interview-question document are missing. The new protocol is
provisional, not an assertion that those deliverables were reviewed.

## B. Implementation plan and completed preparation

The offline implementation is in `src/brd4_analysis/phase2/`:

| Component | Completed behaviour |
|---|---|
| Configuration | Strict keys/types/ranges, one-seed/one-sample envelope, input hashes, fixed construct, explicit MSA policy, no structural templates/constraints or affinity inference |
| Data preparation | Verify construct and 25K chemistry; generate model input, residue/atom maps, frozen metrics/command/provenance; refuse an existing experiment directory |
| Interface audit | Inspect the exact checksum-verified Boltz 2.2.1 wheel's CLI, parser and writer without importing or running it |
| Pose analysis | Complete protein/ligand correspondence, proper protein fit, exact graph symmetry and no ligand-only fitting |
| Contacts | Residue precision/recall at the declared thresholds, including missing/extra contacts and empty-set behaviour |
| Geometry | Graph/coverage/finite-coordinate checks, bond lengths/angles, aromatic planarity, intra/intermolecular heavy-atom clash screens |
| Runtime provenance | Future host preflight checks model hashes, versions and atom names; records MSA/model/software/GPU information without launching Boltz |
| Research reconciliation | Preserve original CSV values; distinguish locally reproduced calculations from external claims; reproduce the 32SC/9FWX duplication |

Safe commands that do not run a model:

```text
.venv/Scripts/python.exe -m brd4_analysis.phase2.cli validate
.venv/Scripts/python.exe -m pytest
```

The prepared bundle already exists at
`experiments/phase2/4men_single_sequence_s20261019/`. `prepare` deliberately
refuses to overwrite it. A changed condition requires a new config ID and
a newly reviewed bundle. Research audit outputs are in `data/phase2/review/`.
New snapshots are separate from the Phase 1 input manifest.

Remaining stages are: receive/reconcile the two missing deliverables;
review this packet and the frozen single-run design; select and freeze a
separate Linux GPU environment; after approval, acquire/verify its assets
and execute the single run; analyse and report it, including any failure.
A multi-seed or MSA-comparison campaign is a later decision.

## C. Validation and limitations

The full offline suite passes **72 tests**: 45 existing Phase 1 tests and
27 Phase 2 tests. They include synthetic rigid-transform and symmetry
controls, a ligand displaced far from the pocket, incomplete protein
atoms, invalid geometry, configuration/hash errors, no-network preparation
and immutable output handling. The rigidly transformed reference scores
near-zero ligand/protein RMSD and contact precision/recall of 1; the
displaced ligand fails pose/contact acceptance. These are synthetic
controls, not Boltz results.

`cli analyse` and `cli report` were exercised using temporary Phase 1
outputs. The LF-default simulation changes the bytes of all eight CSVs;
the explicit-CRLF policy reproduces their bytes. Summary JSON remains
byte-identical, and rendered report text remains identical. This checks
the identified line-ending issue on Windows; it is not a native Linux
validation. Details:
[`phase1_line_endings_check.json`](../../data/phase2/review/phase1_line_endings_check.json).

The local 32SC check compares 1,081 protein/ligand heavy atoms, including
127 Cα and 19 ligand atoms: all compared coordinates equal 9FWX.
32SC's 2026 release date is not an independent test-case date; the earlier
complex was public on 7 August 2024. Both remain post-cutoff by the
reported experimental PDB date, but neither has verified exclusion from
every training source.

The prepared reference reproduces Phase 1's 5/11/13 contact-residue counts
and passes the declared geometry screen. The source inspection is not a
real Boltz parser/GPU integration test. Runtime dependency resolution,
memory use, elapsed time, stochastic variance, confidence calibration
and prediction quality are still unmeasured.

## D. Single-run GPU specification

See [PROTOCOL.md](PROTOCOL.md) for exact definitions, command and provenance.
The frozen condition is **4MEN, 127 residues including SM, neutral 25K
by SMILES, no MSA, one seed `20261019`, one sample, 3 recycles and 200
sampling steps**. No affinity head is requested. The single-sequence
condition tests the harness and is not Boltz's full-MSA performance claim.

Primary acceptance requires ligand RMSD ≤2 Å after the fixed protein fit,
protein Cα RMSD ≤2 Å, residue precision and recall each ≥0.70 at 4 Å,
and all geometry checks. Confidence is diagnostic only. Report the sole
pose even if it fails; do not retry or select a different seed without a
new review.

The result would support a narrow technical demonstration. It would not
establish generalisation, prospective discovery, affinity accuracy,
chemical sample identity, or resolution of the 4MEO dispute.

## E. Resources and explicit stop conditions

These are planning estimates, not measured performance or a provider quote:

| Resource | Proposed envelope |
|---|---|
| GPU | One CUDA NVIDIA GPU; prefer a 24 GB L4/A10/A5000-class device. Fit on smaller GPUs is unmeasured |
| Host | 4 CPU cores; at least 16 GB RAM, preferably 32 GB |
| Disk | Plan roughly 40 GB total; retain at least 25 GB free at runtime preflight |
| Model/cache download | Exactly 6,204,362,719 bytes (6.20 decimal GB) advertised for the three upstream assets; additional environment and extracted-CCD storage required |
| Setup | Rough allowance 15–45 minutes, depending on network, image and cached dependencies |
| Inference | Rough allowance 2–15 minutes for this small complex; unbenchmarked, especially with optimized kernels disabled |
| Hard inference limit | 30 minutes, one job, one sample; no automatic retries |
| Proposed paid window | At most one GPU-hour including setup; prefer completing downloads before reserving the GPU |
| Proposed total spend cap | US$5, subject to explicit approval and a host quote that fits compute, disk and transfer costs |

For a quoted GPU rate `r` dollars/hour, the 30-minute inference ceiling
alone costs at most `0.5 × r`; setup/storage/transfer may add costs.
No provider, current price, paid reservation or actual charge exists.
The advertised model asset sizes come from metadata; the files themselves
were not fetched.

Stop under any of these conditions:

1. User approval for environment/assets and the single paid run is absent,
   or the supplied protocol has not been reconciled with the chosen design.
2. The selected host, runtime freeze, funding limit or enforced timeout
   cannot be established before launch.
3. A source/model/MSA/config hash, construct, ligand graph or runtime atom
   mapping differs from the reviewed specification. Do not silently
   upgrade a package, change a tautomer, drop residues or repair the pose.
4. Model assets are missing, corrupt or from mutable/unverified URLs.
   Preflight must fail before Boltz can attempt its automatic downloads.
5. Resources fall below the approved allocation, CUDA initialization fails,
   an out-of-memory error occurs, or the wall-time/cost ceiling is reached.
   Terminate and preserve logs; do not change steps, seeds, precision or
   hardware and silently continue.
6. Output is incomplete, non-finite, chemically inconsistent or lacks
   provenance/confidence files. Record an invalid run rather than a
   selectively scored subset.
7. The single run completes, whether it passes or fails. Review the
   metrics and actual resource record before any further run.

For the interview, the defensible fallback already exists: demonstrate
the frozen construct/chemistry, show why a ligand-only fit would be wrong,
show the synthetic failure control, and explain the independently
reproduced 32SC duplication and the unresolved 4MEO identity.
Prefer to finish any approved pilot and freeze the narrative by Sunday
18 October rather than starting a campaign on interview morning.
