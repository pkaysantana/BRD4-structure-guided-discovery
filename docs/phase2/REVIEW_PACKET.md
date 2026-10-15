# Phase 2 review packet

**Historical preparation packet.** Phase 1 has since been verified as
merged on `main` at `d0b9f63`, including `cba7549` and its successful
Windows/Ubuntu CI. See [PR_REVIEW.md](PR_REVIEW.md) for the current Phase 2
review scope and validation. The dated recommendations below record the
earlier preparation stage; they do not request reverting the merge.

Prepared on 15 October 2026 for the Lilly AI for Chemistry interview on
Monday 19 October 2026 at 11:30 am. **No GPU job or Boltz prediction has
run; no model weights have been installed or downloaded.**
The final follow-up status and exact minimal primary configuration are in
[FOLLOW_UP.md](FOLLOW_UP.md).

All six supplied deliverables have now been read and reconciled.
See [PROTOCOL_RECONCILIATION.md](PROTOCOL_RECONCILIATION.md) for the
Downloads filenames, preserved hashes and explicit decisions.

## A. Phase 1 merge recommendation

**The CSV correction is ready locally; PR #1 must not be merged yet.**
The corrected writer explicitly uses CRLF at all five CSV write sites.
It is isolated on `phase-1-baseline` in the separate
`BRD4-phase1-csv-portability` worktree, together with a Windows/Linux CI
workflow. The 45 existing Phase 1 tests pass against that corrected code.
Actual writer checks under both Windows and simulated LF defaults preserve
all nine numeric files. The standalone [patch](phase1_csv_fix.patch)
remains the minimal writer/comment change against `5c06db2`.

Linux CI is configured, not run: this host has no installed WSL or Docker
runtime, and nothing has been pushed or dispatched. Native Linux CI and
user approval remain merge gates. The remote PR is unchanged.

The deposited-model scientific scope is defensible. The 4MEO identity
uncertainty remains unresolved and need not block that baseline if it is
retained. See [the merge review](PHASE1_REVIEW.md) and
[the corrected research integration](../research/REVIEWED_INTEGRATION.md).
PR #1 was neither modified nor merged. The isolated source fix changes
serialization policy only; Phase 1 input/result/report bytes and the
supplied research originals are preserved.

## B. Implementation plan and completed preparation

The offline implementation is in `src/brd4_analysis/phase2/`:

| Component | Completed behaviour |
|---|---|
| Configuration | Strict keys/types/ranges, one-seed/one-sample envelope, input hashes, fixed construct, primary cached-MSA versus smoke roles; draft MSA fields allowed only by the non-executing `plan` command |
| Data preparation | Verify construct and 25K chemistry; generate model input, residue/atom maps, frozen metrics/command/provenance; refuse an existing experiment directory |
| Interface audit | Inspect the exact checksum-verified Boltz 2.2.1 wheel's CLI, parser and writer without importing or running it |
| Pose analysis | Complete protein/ligand correspondence, proper protein fit, exact graph symmetry and no ligand-only fitting |
| Contacts | Residue precision/recall at the declared thresholds, including missing/extra contacts and empty-set behaviour |
| Geometry | Graph/coverage/finite-coordinate checks, bond lengths/angles, aromatic planarity, intra/intermolecular heavy-atom clash screens |
| Runtime provenance | Future host preflight checks model hashes, versions and atom names; records MSA/model/software/GPU information without launching Boltz |
| Research reconciliation | Preserve original CSV values; distinguish locally reproduced calculations from external claims; reproduce the 32SC/9FWX duplication |

Safe commands that do not run a model:

```text
.venv/Scripts/python.exe -m brd4_analysis.phase2.cli plan
.venv/Scripts/python.exe -m brd4_analysis.phase2.cli validate --config configs/phase2/4men_smoke.json
.venv/Scripts/python.exe -m pytest
```

The optional smoke bundle is
`experiments/phase2/4men_smoke_s20261019_v2/`. The old single-sequence bundle
is preserved as a superseded preparation record. The primary design is
`configs/phase2/4men_primary.json`, experiment `4men_msa_s20261019_v1`;
it cannot be prepared until its real MSA, hash and provenance are supplied.
No real A3M has been generated, downloaded or fabricated.

`prepare` refuses to overwrite an existing experiment. Latest research
audit outputs are in `data/phase2/review/followup/`, separate from the
Phase 1 input manifest.

Remaining stages are user review, approved acquisition/freezing of the
actual MSA and a separate Linux GPU environment, and approved acquisition/
verification of model assets. Only then prepare and run the single primary
pilot, analyse it and report any failure. A multi-seed or MSA-comparison
campaign is a later decision.

## C. Validation and limitations

All original **45 Phase 1 tests and 27 targeted Phase 2 tests are retained**.
Four additional tests cover pending primary-MSA gates, caching without
network access, rejection of a query-only primary MSA, and changed pinned
alignments. The targeted Phase 2 suite passes all 31 tests. Tests include synthetic rigid-transform and symmetry
controls, a ligand displaced far from the pocket, incomplete protein
atoms, invalid geometry, configuration/hash errors, no-network preparation
and immutable output handling. The rigidly transformed reference scores
near-zero ligand/protein RMSD and contact precision/recall of 1; the
displaced ligand fails pose/contact acceptance. These are synthetic
controls, not Boltz results.

After stacking the Phase 2 commits on the isolated Phase 1 correction,
the combined offline suite passes **76 tests**. All original baseline
input/result/report bytes, all six supplied-source hashes in Git and the
current smoke bundle were checked after the local rebase. The recorded
result and source fingerprints are in
[`validation.json`](../../data/phase2/review/followup/validation.json).

`cli analyse` and `cli report` were exercised using temporary Phase 1
outputs using the **actual corrected writer**. It reproduces the eight
CSVs and summary JSON under both default-EOL policies. The rendered report
bytes match its committed LF form and the original workspace report;
fresh Windows checkouts may use CRLF, and were not overwritten.
This is not a native Linux validation. Details:
[`phase1_validation.json`](../../data/phase2/review/followup/phase1_validation.json).

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
The primary design is **4MEN, 127 residues including SM, neutral 25K
by SMILES, cached full MSA (up to 8192 rows), one seed `20261019`, one
sample, 3 recycles and 200 sampling steps**. Step scale is 1.5; there are
no templates, constraints, steering potentials or affinity request.
The MSA bytes/hash remain an explicit preparation gate. The separate
no-MSA configuration is an optional smoke test, not the primary result.

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

1. User approval for environment/assets, MSA acquisition or the single paid
   run is absent, or the real primary MSA/provenance is not pinned.
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
