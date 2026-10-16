# Phase 2 independent review

Phase 1 is merged on `main` at
`d0b9f63c53fc1aff5d29273191ee940a01e111be`. The merge contains
`cba754944b8ce80028fb9e0c2602fcbc7c35ff0c`, and its complete tree matches
that reviewed Phase 1 head. Both native Windows 2022 and Ubuntu 24.04 CI
jobs passed all 45 Phase 1 tests, numeric byte comparisons and
regeneration/drift checks.

This review is for **offline Phase 2 preparation and analysis utilities**.
There are no Boltz predictions or completed Legion hardware measurements.
The ThinkPad must not be used for GPU compatibility diagnostics.
Opening or approving this code PR does not authorise model downloads,
inference, MSA acquisition or GPU hosting.

## Contents and review order

1. `src/brd4_analysis/phase2/config.py`, `inputs.py` and `workflow.py`:
   validate the fixed 4MEN construct and 25K graph, generate versioned
   Boltz input, freeze provenance and metrics, and refuse existing output
   directories. Cached full MSA is the primary condition; no MSA is an
   explicitly labelled smoke test. Planning with pending MSA fields does
   not permit strict preparation or execution.
2. `metrics.py` and `tests/test_phase2.py`: complete protein/ligand
   correspondence, proper protein alignment, graph-symmetry ligand RMSD
   without a ligand fit, residue contact precision/recall and a limited
   heavy-atom geometry screen. These are operational metrics, not
   binding-affinity or physical-stability estimates.
   Atom-coverage diagnostics distinguish optional terminal OXT from valid
   predicted side-chain atoms that lack reference coordinates. The
   primary contact metric is unchanged; a separate reference-observed-atom
   diagnostic exposes the effect of incomplete experimental side chains.
3. `runtime.py`: a future, explicit Linux-host provenance check that
   verifies existing model assets and records the environment. It does
   not install packages, download assets or invoke Boltz.
4. `scripts/boltz_gpu_preflight.py` and `tests/test_gpu_preflight.py`:
   standalone Legion-only compatibility proposal, with a separate
   environment and synthetic CPU-only tests. The hardware collector and
   tiny CUDA operations have not been run. The intended laptop has
   user-reported 12 GB VRAM; neither fit nor timing has been measured.
5. `docs/research/`, `docs/protocols/` and `data/phase2/`: all six supplied
   deliverables, corrected scientific interpretations, hashed reference
   inputs and upstream source/metadata snapshots. Upstream Python excerpts
   are passive `.txt` evidence, not imported execution modules.
6. `configs/phase2/` and `experiments/phase2/`: the blocked primary design,
   an unexecuted smoke bundle, and clearly labelled superseded bundles.
   Existing experiment records are preserved rather than overwritten.

The current protocol is [PROTOCOL.md](PROTOCOL.md), with explicit
departures from the supplied draft in
[PROTOCOL_RECONCILIATION.md](PROTOCOL_RECONCILIATION.md).
The [Legion guide](LEGION_GPU_PREFLIGHT.md) is preserved byte-for-byte
from the earlier proposal; its introductory Phase 1 publishing statements
predate the merge. Its operational preflight instructions remain a
proposal for the powered-on Legion.

## Preservation and validation

The two original local Phase 2 commits and the preflight checkpoint are
retained under the local backup branch
`backup/phase2-pre-main-d0b9f63`. Rebase range comparison confirmed that
all three patches are unchanged:

| Before integration | Rebased commit |
|---|---|
| `d53ddbe` — offline harness | `296177a` |
| `f5ce0c7` — protocol/MSA reconciliation | `a3bae34` |
| `34277a2` — three Legion preflight files | `5424df3` |

The script, its tests and its guide match their original working-file
SHA-256 values. All 22 protected Phase 1 input/result/report files retain
their bytes. No Phase 1 implementation, existing test, dependency pin,
CI workflow or golden result is changed by this PR. The only existing
`main` file with edits is the README, adding Phase 2 navigation.
Historical Phase 1 review receipts and the old writer patch are retained
as evidence, not reapplied changes.

The integration checkpoint passed **84 tests**, and the atom-coverage
checkpoint passed **89 tests**. After the three final-review corrections,
the full local Windows suite passes **161 tests**, with one symlink test
skipped because the current Windows account lacks symlink privileges.
The suite contains 162 cases: 45 Phase 1, 44 Phase 2 preparation/scoring,
65 runtime-integrity/provenance and eight synthetic GPU-preflight cases.
The original 27 Phase 2 tests remain.
GPU-preflight tests block subprocess calls and Torch imports; they do not
call the hardware collector. No GPU diagnostics or inference were needed
for these results. The original integration receipt is
[`main_integration.json`](../../data/phase2/review/main_integration.json);
the earlier atom-coverage correction is recorded in
[`codex_atom_coverage.json`](../../data/phase2/review/codex_atom_coverage.json).
The final corrections and validation are recorded in
[`final_review_corrections.json`](../../data/phase2/review/final_review_corrections.json).
Runtime tests use tiny synthetic assets and mocked GPU services; they
neither import real Torch nor execute hardware commands.

## Final review corrections

The focused review of `5ebd95a` identified two P1 validation defects and
one P2 provenance defect. The corrections are limited to those boundaries:

- Extracted CCD files and directories must match the inventory and bytes
  in the checksum-verified `mols.tar`. Missing/altered/unexpected entries,
  unsafe paths, links, special files and duplicate normalized archive
  entries are refused before GPU probes.
- A primary MSA needs aligned standard amino-acid information in at least
  one non-query record within the declared 8192-record prefix. Query plus
  padding cannot qualify. Lowercase insertions and query/alignment lengths
  retain their checks; this does not independently establish homology.
- Scientific scoring requires a complete versioned runtime manifest,
  its capture-time checksum receipt, and consistency with frozen inputs,
  command, MSA, model assets/revision and relevant package requirements.
  Incomplete, changed or inconsistent records are refused before scoring.
  The receipt detects changes; it is not attestation against forging both
  the record and receipt, or proof of subsequent manual execution.

All existing scientific/reference artefacts and the original Legion
preflight files are preserved. Current smoke preparation is v4, which
freezes the corrected code hashes; earlier bundles remain historical.
Metric definitions, thresholds, ligand/construct inputs, model pins and
the primary full-MSA condition have not changed.

## Independent review correction

Codex identified that the reference contains terminal GLU168 OXT, while
Boltz 2.2.1's canonical standard-residue output omits it. A regression
reproduced the rejection before the fix. Upstream atom lists were then
checked against all 127 reference residues: Lys91 and Gln123 also have
unresolved distal atoms that a complete prediction includes.

The adapter now permits missing OXT only at the terminus, accepts additional
atoms only through the pinned reference's residue-type name/element
templates, and reports both coverage differences. The templates were
checked against the pinned upstream canonical sets. No experimental
coordinates are invented. All other observed heavy atoms, including
terminal backbone O, remain required.

Primary contacts still use all predicted heavy atoms, with the original
cutoffs and thresholds. The additional contact diagnostic is descriptive
and cannot replace the primary acceptance score. Geometry uses all
predicted atoms. Five added regression cases cover OXT omission, retained
terminal-O requirements, canonical additional atoms/contact accounting,
and rejection of unknown or wrongly typed extras.

That source policy was frozen in
`experiments/phase2/4men_smoke_s20261019_v3/`, now preserved as a superseded
checkpoint. V4 retains it along with the final validation corrections. The
primary full-MSA design remains blocked on real MSA input/provenance.
Phase 1 and the three original Legion preflight files remain unchanged.

## Scientific and execution boundaries

- 4MEN is a retrospective known-complex control. Pre-cutoff release means
  training eligibility, not verified membership or demonstrated recall.
  Related complexes, distillation and affinity-label exposure remain risks.
- 32SC and 9FWX are one duplicated public complex, not two independent
  evaluation observations. The earliest established public complex date
  is used in the assessment.
- The 4MEO sample identity remains unresolved. Deposited models,
  external map/refinement computations, reproduced geometry, predictions
  and hypotheses are kept distinct.
- The primary cached A3M, its checksum and provenance are still absent.
  The query FASTA and synthetic unit-test alignments are not substitutes.
- Boltz 2.2.1 syntax has been checked against its versioned upstream
  source. An installed Boltz/GPU end-to-end run has not been validated.
- The proposed Legion environment remains separate from Phase 1.
  `--no_kernels` does not cure an incompatible CUDA build, and a tiny
  CUDA probe would not establish model memory fit.
- A later first primary run would use one seed, one diffusion sample,
  one parallel sample, no affinity and the full cached MSA. A depth
  reduction, new seed or remote GPU requires a separately reviewed
  experiment; none is scheduled here.

Independent review should focus on the MSA preparation gate, atom/residue
mapping, avoidance of ligand refitting, output/provenance integrity and
the absence of unintended execution paths. The existing scientific and
hardware limitations must remain visible even when synthetic tests pass.
