# Follow-up and final preparation

15 October 2026. Local work is ready for review. No push, merge, CI
dispatch, MSA search, model-weight download or GPU inference was performed.

## 1. Located files and reconciliation

Both files were found in `C:\Users\LLR User\Downloads\`, without the
`_claude_science` suffix in their actual filenames:

- `03_boltz2_protocol_4MEN.md` → archived as
  `docs/protocols/03_boltz2_protocol_4MEN_claude_science.md`.
- `04_interview_questions.md` → archived as
  `docs/research/04_interview_questions_claude_science.md`.

They were read in full and copied byte-for-byte. All six deliverables
are now accounted for. The
[reconciliation table](PROTOCOL_RECONCILIATION.md) records decisions and
the [reviewed interview guidance](../research/04_interview_questions_reviewed.md)
corrects the reused scientific overstatements.

Full MSA is the primary Arm A condition. The no-MSA run is an optional
smoke test, not the primary scientific result. 32SC and 9FWX remain one
duplicated public complex. Date eligibility is not verified training
membership; model maps, retained refinement models and predictions do not
identify the species in the 4MEO crystal.

## 2. Phase 1 correction

Local commit **`cba7549`**, on `phase-1-baseline`, contains exactly:

- Explicit CRLF at all five CSV write sites.
- The matching `.gitattributes` comment correction.
- A Windows 2022 / Ubuntu 24.04 CI workflow using Python 3.12.10,
  pinned official actions, the existing offline tests and result-drift
  verification.

The correction was developed in the separate
`C:\Users\LLR User\BRD4\BRD4-phase1-csv-portability` worktree.
All 45 existing Phase 1 tests pass. The actual patched writer reproduces
all eight CSVs and summary JSON byte-for-byte under both Windows and
simulated LF defaults. Temporary report rendering matches the committed
LF form and the original workspace's report bytes. Existing data/results/
report files were preserved.

Native Linux CI is **pending**, not reported as passed. WSL is not
installed and Docker is unavailable on this host. The workflow will run
only after an approved remote update or dispatch. The minimal
[writer patch](phase1_csv_fix.patch) remains available against `5c06db2`;
the local commit also includes CI. Do not merge without approval and a
successful native Linux check.

## 3. Git state and isolation

`phase-1-baseline` is one local commit ahead of its origin tracking branch.
`phase-2-controlled-pilot` has two separate Phase 2 commits stacked on that
corrected local base; it has no upstream branch. A local rebase established
that dependency without merging either PR/branch.

The Phase 2 review diff is:

```text
git diff phase-1-baseline..phase-2-controlled-pilot
```

It contains no Phase 1 writer, `.gitattributes`, CI-workflow or existing
result/report changes. The Phase 1 fix can be reviewed independently.
Both worktrees are left clean after recording this report. Remote refs
have not been updated.

## 4. Reviewed-push readiness

**Ready for a reviewed push of the prepared commits, subject to user
approval.** The combined offline suite passes **76 tests**: 45 existing
Phase 1 tests, all 27 original targeted Phase 2 tests, and four new MSA
tests. The original test function bodies were retained. The four new tests
cover the pending-primary-input gate, cached-MSA copying/provenance without
network calls, rejection of a query-only primary alignment and rejection
of a changed pinned alignment even during planning.

[Validation evidence](../../data/phase2/review/followup/validation.json)
records the result and source hashes. All original Phase 1 inputs/results/
report bytes and all six supplied-source hashes in Git were verified.
The smoke bundle also verifies after the rebase.

Inference remains blocked. There is no real A3M, no pinned MSA checksum,
and no selected/resolved GPU environment. Boltz and Torch are not installed
in the Phase 1 environment. Null MSA fields are explicit pending inputs;
they are not filled with a synthetic or query-only substitute.

## 5. Exact minimal GPU pilot

Primary design:
[`configs/phase2/4men_primary.json`](../../configs/phase2/4men_primary.json).
It is an unexecuted design until the actual alignment/provenance is bound.

| Parameter | Value |
|---|---|
| Experiment | `4men_msa_s20261019_v1` |
| Purpose | First Arm A timing and structural evaluation of the known 4MEN complex |
| Software | Boltz 2.2.1; RDKit 2026.3.6 for the verified ligand atom mapping |
| Runtime | Separate Linux / Python 3.12 GPU environment; Torch/CUDA, driver and full environment freeze pending |
| Protein | Chain A, 127 residues: SM expression remnant + O60885 residues 44–168 |
| Ligand | Chain B, neutral 25K, 20 heavy atoms, SMILES `Cc1ccc(CN(C)c2cc(C)nc3ncnn23)cc1` |
| MSA | Real cached full A3M at `inputs/protein.a3m`, exact query and SHA-256/provenance required |
| MSA limit | 8192 sequences; no depth-subsampling ablation |
| Model | `boltz2`; model repository revision `6fdef46d763fee7fbb83ca5501ccceff43b85607` |
| Structure checkpoint SHA-256 | `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1` |
| Seed | `20261019` |
| Recycling steps | 3 |
| Sampling steps | 200 |
| Diffusion samples / parallel samples | 1 / 1 |
| Step scale | 1.5 |
| Devices / accelerator | 1 / GPU |
| Data-loader workers / preprocessing threads | 0 / 1 |
| Optimized kernels | Disabled via `--no_kernels` |
| Templates / constraints / steering / affinity request | None |
| Output | mmCIF, sole `model_0`, confidence JSON; separate experiment directory |
| Resource proposal | Prefer 24 GB VRAM; ≥16 GB host RAM; ≥25 GB free disk at preflight |
| Limits requiring approval | 30-minute inference ceiling, one paid GPU-hour overall, US$5 proposed total cap; no automatic retry |

The exact prediction command is recorded in the [protocol](PROTOCOL.md).
It uses the local `msa:` file and does not query the server during
inference. The official immutable upstream guide and CLI/parser code were
re-fetched and match the saved v2.2.1 source. Static syntax verification
is not a claim of live runtime testing.

Safe review commands:

```text
.venv/Scripts/python.exe -m brd4_analysis.phase2.cli plan
.venv/Scripts/python.exe -m brd4_analysis.phase2.cli validate --config configs/phase2/4men_smoke.json
```

`prepare` with the default primary config must fail until its real MSA and
provenance are supplied. Its planned directory does not yet exist.
`4men_smoke_s20261019_v2` is separately prepared but remains optional and
unexecuted. No campaign is scheduled or implied by either configuration.
