# Phase 4 — adversarial validation of the Phase 3 acceptance gates

Status: implemented and executed (CPU-only, offline). The post-audit
`post_audit_v1` battery contains 47 deterministic cases and is reproduced
byte-for-byte in the recorded WSL environment by
`.venv/bin/python -m brd4_analysis.phase4.cli battery --out
experiments/phase4/adversarial_screen_v1`; the pass/fail matrix is
`experiments/phase4/adversarial_screen_v1/FINDINGS.md` (generated) and full
values are in `results.json`.

## What was tested

The five registered Phase 3 gates evaluated by `score_prediction`
(`src/brd4_analysis/phase2/workflow.py`):

| gate | registered limit |
|---|---|
| ligand RMSD (protein-frame, symmetry-aware, no ligand fit) | ≤ 2.0 A |
| protein CA RMSD (125 fit CAs, proper-rotation Kabsch) | ≤ 2.0 A |
| contact precision (residue-level, 4.0 A, primary cutoff) | ≥ 0.7 |
| contact recall (same) | ≥ 0.7 |
| geometry screen | bond ≤ 0.25 A, angle ≤ 30 deg, aromatic plane ≤ 0.15 A, intraligand vdW ratio ≥ 0.75, zero protein-ligand overlaps ≥ 0.8 A |

Baseline (unmodified committed prediction): ligand RMSD 0.1405 A, CA RMSD
0.1914 A, 4.0 A contacts P/R = 1.0/1.0 — all gates pass.

The battery applies fixed transforms to this prediction (or replays the
reference) and scores each through a replica of `score_prediction` that is
pinned field-for-field against the committed `analysis/metrics.json` and the
real scorer. Permanent end-to-end tests also serialize faithful predicted
complexes — retaining actual predicted protein coordinates and atom names —
through the unmodified scorer. Covered cases include the 12% ligand
expansion, the signed +20° angle opening, and a protein expansion that must
fail the protein-CA gate.

## Post-audit battery revision

The originally drafted battery had 45 cases. The corrected `post_audit_v1`
battery has 47 cases because the audit exposed two coverage gaps:

- the old `angle_bend` sign convention treated a positive argument as angle
  closing; `post_audit_v1` defines `degrees` as the intended signed valence
  change and explicitly tests −35°, −20°, +20° and +35°;
- `protein_scale_1.25` was added as a permanent protein-side case that fails
  the registered CA-RMSD gate.

This is a corrected exploratory battery, not a backdated preregistration:
the post-audit additions were designed after the earlier results had been
reviewed.

## Reproduced finding

A uniform **+12% expansion of the predicted ligand about its centroid**
passes all five registered gates (`centroid_scale_1.12`):

| metric | original prediction | 12% expansion |
|---|---:|---:|
| ligand RMSD | 0.1405 A | 0.3972 A |
| protein CA RMSD | 0.1914 A | 0.1914 A |
| 4.0 A contact precision / recall | 1.000 / 1.000 | 0.9167 / 1.000 |
| maximum bond error vs CCD | 0.0541 A | 0.2048 A |
| minimum intraligand vdW ratio | 0.7912 | 0.8862 |
| severe protein-ligand overlaps | 0 | 0 |

Why it works: uniform scaling preserves angles and molecular handedness and
scales every pairwise distance, including every covalent bond and nonbonded
separation. Exact coplanarity is preserved, while nonzero absolute planarity
deviations scale with the molecule. Here intraligand clashes move farther
apart (minimum vdW ratio improves), while correlated bond-length drift of
~12% remains below the absolute 0.25 A bond gate. The residue-level contact
set is quantised and does not respond continuously to sub-Å coordinate
changes.

The continuous accepted uniform-scaling interval for **this prediction and
these thresholds** is approximately **−5.21% to +15.33%** (scale factors
0.947888–1.153270). The lower boundary is the intraligand vdW-ratio gate;
the upper boundary is the worst-bond absolute-error gate. The sampled
positive factors 1.02–1.15 pass, but that sampled range is not the complete
acceptance interval.

Caveat found during implementation: the same +12% expansion applied to the
**deposited** ligand coordinates is caught by the bond gate, because the
deposited model's own bond lengths already deviate up to 0.118 A from the
CCD template. The finding is specific to expanding a near-template
prediction.

## Accepted synthetic variants

Acceptance means the registered gates admit the transformed coordinates. It
is not an empirical model error rate, and a rigid displacement is not
automatically a chemically invalid structure.

| case | interpretation |
|---|---|
| `centroid_scale_1.02`–`1.15` | sampled positive inflations inside the wider verified acceptance interval |
| `ligand_translate_0.25`–`1.00A` | rigid ligand displacements remain under the 2.0 A pose tolerance and retain the contact set |
| `ligand_rotate_15deg` | a rigid 15° reorientation gives pose RMSD 0.5943 A and remains accepted |
| `bond_stretch_C5-N_+0.20A` | one covalent bond is 0.20 A too long (1.4587→1.6587 A), below the 0.25 A absolute gate; many cross-fragment nonbonded distances also change |
| `angle_bend_C5-N-C7_+20deg` | the corrected signed +20° opening passes all five gates: pose RMSD 0.5525 A, max angle error 23.2019° (< 30°), contacts 1.0/0.9091, zero severe overlaps |
| `protein_scale_1.05` | 5% predicted-protein inflation is accepted: CA RMSD 0.7399 A (< 2.0), pocket CA RMSD 0.6919 A |
| `jitter_sigma_0.05` | fixed-seed σ = 0.05 A uncorrelated noise is accepted; this may be a legitimate tolerance rather than an invalid structure |

## Cases correctly rejected

Compression at or beyond the tested 15% cases, sampled expansions ≥20%,
translations ≥1.5 A, rotations ≥45°, 90°/180° torsion flips about either
linker bond, bond changes ≥ |0.30| A, angle changes −35°, −20° and +35°,
anisotropic 1.25× scaling, jitter σ ≥ 0.20 A, clash-directed translations,
25% protein inflation, and the ligand mirror image are rejected by at least
one registered gate. Rejection directions are not symmetric: the +20° angle
opening is accepted, while the −20° closing move creates a severe
protein-ligand overlap and fails geometry.

The mirror image is rejected through pose and clash consequences, not
through a chirality-aware internal check; the pilot is declared achiral-only.

## Candidate checks (exploratory — not registered, not adopted)

Seven stronger checks are implemented in `src/brd4_analysis/phase4/checks.py`
and evaluated on every case. Their limits were chosen on the same committed
complex/control set being screened, so this battery provides no independent
sensitivity, specificity or legitimate-structure false-rejection estimate.
They are diagnostics for future preregistered evaluation, not calibrated
general-purpose validity criteria.

| proposed check | scope classification | exploratory limit | catches in this screen | limitation / notable miss |
|---|---|---:|---|---|
| `bond_signed_bias` | conditional: reference-free prospective when the predicted graph has a trusted CCD template; otherwise requires a suitable reference | abs(mean) ≤ 0.04 A | correlated inflation/deflation and larger jitter | opposing bond errors cancel; misses the isolated +0.20 A bond stretch |
| `bond_relative_error` | conditional as above | ≤ 10% | ≥ ~8–9% correlated inflation; isolated bond errors | sensitive to single-coordinate/template error; flags σ = 0.05 jitter |
| `similarity_scale` | reference-dependent retrospective | abs(ln s) ≤ 0.03 | scaling, large conformational changes and reflection | requires deposited coordinates; conformational changes also alter fitted scale |
| `internal_distance_rmsd` | reference-dependent retrospective | ≤ 0.15 A | conformational and scaling changes | measures agreement with one conformer; can falsely reject a legitimate alternate conformer; misses +0.20 A bond stretch at 0.1425 A |
| `contact_distance_bias` | reference-dependent retrospective | ≤ 0.10 A | systematic radial/contact shifts and protein inflation | requires the reference complex; signed shifts and changing nearest pairs complicate interpretation |
| `pose_rmsd_tightened` | reference-dependent retrospective | ≤ 0.5 A | larger pose displacements | reconstruction criterion only; tightening may reject legitimate tolerated poses; misses +0.20 A stretch |
| `pocket_ca_rmsd` | reference-dependent retrospective protein-fit check | ≤ 1.0 A | the added 25% protein-inflation control | receptor-side measure only; misses 5% protein inflation and ligand-side chemistry by construction |

The proposed four-check combination — `similarity_scale`,
`internal_distance_rmsd`, `bond_signed_bias` and `pose_rmsd_tightened` — is
not comprehensive. It accepts both `bond_stretch_C5-N_+0.20A` and
`protein_scale_1.05`. It also uses deposited-coordinate checks, so it is a
retrospective reconstruction diagnostic rather than a fully
reference-independent prospective validator.

## Scientific limitations

- **Synthetic screen, not a failure-rate estimate.** Every non-control is a
  deterministic coordinate transform of one prediction on one target.
  Acceptance of a synthetic variant says what this scoring harness tolerates;
  it does not estimate how often Boltz-2.2.1 or another model produces such
  structures.
- **Accepted does not mean chemically invalid.** Rigid translations and
  rotations preserve all internal interatomic distances exactly. Low-amplitude
  independent coordinate jitter changes bond lengths, bond angles and signed
  aggregate statistics, although those changes may remain within registered
  tolerance. The +12% expansion, +0.20 A bond stretch and systematic protein
  inflation are clearer implausible-geometry examples; the outcome labels in
  `results.json` keep this distinction.
- **Single-complex scope.** The battery runs against 4MEN/25K only. The
  numerical acceptance interval is specific to this prediction, ligand size,
  symmetry and pocket and must not be generalised.
- **CCD template is a 3D conformer.** The pinned ideal 25K coordinates are
  three-dimensional, not a flat 2D embedding. Whole-molecule comparisons to
  that particular conformer conflate conformational choice with geometry;
  bond lengths and valence angles remain more directly interpretable.
- **Automorphism subtlety.** The registered geometry screen uses the
  pose-optimal graph automorphism; for distant poses this can affect which
  atoms are compared. The candidate checks report the best statistic over
  exact graph automorphisms and disclose that convention.
- **Reflection handling is incidental.** The mirror image was caught by pose
  and clash consequences, not by a chirality-aware check.
- **Protein-side coverage remains narrow.** Two protein-scale cases now
  cover one accepted and one rejected magnitude; they do not cover side-chain
  moves, missing atoms or residue-level adversarial edits.
- **Shared implementation assumptions.** The replica reuses the frozen
  metrics/input code so parity tests cannot independently prove those frozen
  calculations. The end-to-end tests verify that the in-memory evaluator
  invokes equivalent registered logic, not that the frozen metrics are
  universally optimal.
- **Detection floor.** `centroid_scale_1.02` passes all proposed checks as
  well as the registered gates. Correlated inflation below roughly 2–4% is
  below the noise floor set by this ligand's deposited-vs-template bond
  differences.

## Reproducibility and provenance

`results.json` records the battery revision, design status, per-case
parameters, the jitter seed, Python/platform/package versions, input hashes
and LF-normalized Phase 4 source hashes. `manifest.json` records SHA-256
hashes of the generated `results.json` and `FINDINGS.md`. Environment
provenance is excluded from cross-platform scientific-result comparison but
remains recorded truthfully in the published artifact. The CLI confines
output to `experiments/phase4`, rejects aliased roots and symlinked
destinations/artifacts, replaces artifacts atomically without truncating a
hard-linked inode, pins the destination through a directory descriptor on
supported POSIX systems, and returns a nonzero status if a required control
fails. On Windows it retains lexical confinement and atomic inode replacement,
but Python does not expose descriptor-relative rename, so a detected
post-validation directory swap is rejected without claiming a descriptor-level
race guarantee.
