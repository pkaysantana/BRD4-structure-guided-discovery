# Reviewed scientific integration

Review date: 15 October 2026. Baseline: PR #1, commit
`5c06db21811a29e9925ab008ca555129c71c0a35`.
The initial engineering review covered four supplied files. The follow-up
located and read the protocol and interview questions in Downloads,
completing all six. This is not a repeat literature review.
Published results cited by the supplied brief remain
attributed to that source unless a local reproduction is identified below.

The [protocol reconciliation](../phase2/PROTOCOL_RECONCILIATION.md)
records the subsequent decisions: cached full MSA is the primary
condition, the no-MSA run is an optional smoke test, and the first
primary timing run is one seed/one sample. The same evidence corrections
below apply to the protocol and to the
[reviewed interview answers](04_interview_questions_reviewed.md).

## Evidence categories

| Category | Meaning in this project |
|---|---|
| Published report | What the cited authors say they measured or interpreted; assay, sample and method context must travel with the claim |
| Deposited model/metadata | Coordinates, assigned component, occupancies, B-factors and archive metadata; an interpretation of experimental data |
| Locally reproduced calculation | A calculation run here on identified, hashed inputs; independent implementation does not make the underlying experimental data independent |
| External computation | A supplied map sample, PDB-REDO result or validation statistic not reproduced here |
| Prediction | A future Boltz output, with model/input provenance; none has been generated |
| Hypothesis | An explanation to test, with no promotion to a measured result |

## Corrections that govern downstream use

**Model-phased maps.** Withdraw the claim that density sampled at the
putative oxide oxygen is an “upper bound” on unbiased evidence. Model
phases, map weighting and bulk-solvent treatment can affect density in
either direction; there is no demonstrated bound. A small difference-map
residual does not establish an atom's existence or full occupancy.
The supplied map samples and profiles have not been reproduced because
the maps and sampling code were not supplied.

Omit and polder maps can reduce particular sources of model or solvent-mask
bias. They do not provide experimentally independent phases or guarantee
an unbiased, decisive chemical identification. Omitting one atom while
retaining the surrounding refined model leaves dependencies. Competing
models would need consistent processing, refinement/restraints, solvent
treatment, occupancy/B-factor handling and the same held-out reflections,
with local fit and residuals assessed together. No such calculation was
performed here.

Do not treat approximately 0.8 Å map-grid sampling as the experiment's
resolving power or as proof that alternatives cannot be distinguished.
Interpolation and grid spacing differ from diffraction resolution,
signal quality and local identifiability. The supplied profiles neither
prove discrimination nor establish that all public diffraction evidence
is incapable of discrimination.

**PDB-REDO.** Re-refinement of the same diffraction data from a related
starting model is a useful sensitivity check. Retention of the same ligand
restraints and assignment is not independent chemical confirmation or a
comparison against a neutral quinoline/water model. A changed B-factor
contrast does not by itself prove that the original contrast was a
refinement artefact; occupancy, TLS and restraint parameterization are
coupled. Whether occupancy was actually optimized must be checked before
interpreting a retained value of 1.00.

Withdraw the rule that a change in R-free smaller than 0.008 is
“uninterpretable.” An estimated uncertainty on a single R-free is not the
uncertainty on a paired difference between models evaluated on the same
reflections. A small ligand also contributes little to a global score.
The supplied PDB-REDO numbers remain external computations, not local
re-refinement results.

**Sample identification.** The safe statement is: the paper reports no
N-oxide detected in the material it analysed; the deposited 4MEO ligand
model is the N-oxide 25V. Neither statement identifies the chemical species
in the particular crystal. The chain connecting the assayed, MS-tested,
resynthesised, soaked and crystallised material, including lot/vial IDs
and analytical detection limits, has not been established here.
“ChemBridge 50866471” is a supplied attribution, not a verified identifier
of the crystal's sample.

A trace impurity, partial occupancy, oxidation and an incorrect ligand
assignment remain hypotheses. Do not rank impurity selection as the
“most parsimonious” explanation from the present evidence, say it predicts
a mixed site, or assert that crystallization selects by affinity alone.
Concentrations, equilibria, solubility and crystal conditions also matter.
The reported biochemical IC50 of authentic 14a is not a binding constant
under crystallization conditions. Negative MS does not establish purity
below an unspecified detection limit or refute a different crystal sample.

For compounds 3 and 5, retain the supplied statement that a quantitative
IC50 was not determined within the accessible solubility range. Correct
the identity assessment's internally inconsistent “IC50 was lower than
its solubility” wording. Do not invent an IC50 or use either compound as
a precise affinity label. The reported 42 µM value belongs to compound
14a in the publication, not to a validated identification of the 4MEO
crystal ligand.

**Training membership.** The
[Boltz-2 report](https://jeremywohlwend.com/assets/boltz2.pdf) describes a
1 June 2023 cutoff for experimental PDB structural data and additional
distillation/MD and affinity supervision. This is not a universal cutoff
for every training source or a per-checkpoint membership manifest.
4MEN/4MEO are eligible by date and strongly exposed to leakage risk.
“Memorised,” “known training member” and definite affinity-label leakage
are unsupported without exact dataset/checkpoint evidence. The supplied
ChEMBL document attribution is a risk lead, not proof of label inclusion.
The supplied Boltz-1 cutoff is not needed here and is not propagated as a
verified fact.

A 4MEN result can demonstrate a working, measured pipeline on a known
complex. It cannot alone establish generalisation or distinguish recall
from transferable learning. Conversely, calling it incapable of providing
any useful capability evidence is too broad. SMILES avoids requesting the
specific `25K` CCD component as input; it does not erase knowledge of the
graph, related chemotypes, protein fold or binding pocket.

**Later-release evaluation.** 32SC really is recorded as released on
5 August 2026. However, its complete compared protein and ligand
heavy-atom coordinates are identical to 9FWX, released on 7 August 2024.
All 127 Cα atoms and all 19 A1IGX ligand atoms match. The exact-graph RCSB
search returns both entries. See the
[reproduced assessment](../../data/phase2/review/32SC_assessment.json).

Do not present 32SC as a newly independent 2026 case. Treat this complex
once, using its earliest established public date. That date is after the
reported PDB structural cutoff, but before the selected checkpoint's
public release; exclusion from distillation, affinity data and related
complexes is unverified. It is a conditional retrospective case, not a
clean prospective or novel-chemotype benchmark. The pairwise Morgan
similarity to 25K in the assessment is not a search for the maximum
similarity against all pre-cutoff BET ligands.

32SC also has a citation inconsistency: the RCSB/mmCIF DOI differs from
PDBe's publication record, while the latter agrees with 9FWX. Preserve
the discrepancy; do not silently substitute one record for another.

**Controls and calibration.** A single prediction cannot calibrate
confidence or estimate seed variability. Five seeds is a proposed design
choice, not a statistically justified minimum. Different non-binders
without experimental bound poses do not have a meaningful “correct”
ligand pose RMSD against 25K. Use rigidly transformed and displaced
copies of the same graph to test the metric implementation; a future
binding-classification study needs its own controls and assay labels.

## Reconciliation with Phase 1

The reviewed head already corrects the stale 0.154 similarity to 6/72 =
0.0833, identifies `25V:O1` as the acetamide carbonyl and `25V:O` as the
N-oxide, and limits the template check to CCD graph consistency. Its
geometric contact thresholds are 3.5, 4.0 and 5.0 Å. The brief's 4.2 Å
lists use a different definition and should not replace those tables.

The measurements audit retains all 27 supplied rows. Nine rows have local
checks: archive resolution metadata, CCD formula, deposited ligand mean
B and occupancy, the principal protein and water distances, the additional
Asp88 distance, and the deposited N1–O distance/B-factors. Numeric values
are in the machine-readable reconciliation. Checks of a scalar do not
verify every interpretation or annotation in that row. All remaining
rows retain an explicit external/unreproduced status.

Distances near hydrogen-bond ranges remain geometric observations, not
hydrogen-bond assignments without hydrogen/angle and chemical context.
Weak density, model B-factors and local fit do not uniquely imply affinity,
occupancy or chemical identity. The RSCC ordering quoted for the four
Vidler models is an external validation comparison, not a new local map
analysis or proof of an assignment error.

4MEO remains suitable for a carefully labelled deposited-model baseline.
It is excluded from the initial Boltz pose-accuracy test and from any
claim that a model can adjudicate the sample identity. This unresolved
scientific question is not itself a blocker to merging a baseline that
states its scope. The actual outstanding merge issue is documented in
[the Phase 1 review](../phase2/PHASE1_REVIEW.md).
