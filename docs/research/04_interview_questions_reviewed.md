# Five interview questions — reviewed answer guidance

Engineering review of the supplied question document, 15 October 2026.
Published Vidler observations below are attributed to the supplied
research and Phase 1 report. The new local work is the software,
reproduced calculations and evidence audit. No Boltz result exists.

## 1. How strong is the published screening hit-rate comparison?

Describe it as Vidler's result: six hits from 143 compounds in the
structure-based arm, against zero from 97 controls, as reported in the
supplied brief. Six events support a limited comparison; they do not
separate the effects of compound selection, docking and expert triage.
A zero-hit control does not establish a zero underlying hit probability.
Do not claim a significance level or confidence interval that has not
been calculated.

An additional docking-without-manual-triage arm could help study triage,
provided compound selection and assay conditions were controlled.
Define hits and treatment of solubility-limited measurements in advance.
Human inspection is a plausible contributor, not a proven sufficient
explanation. Do not call this project's coordinate analysis a virtual
screen performed by you.

## 2. How would you investigate the 4MEO ligand assignment?

Separate the evidence levels. The deposited CCD model 25V is an N-oxide
graph; the coordinates contain its modelled N1–O bond. The publication
describes compound 5 and reports no N-oxide detected by MS. This is a
documented conflict, not a chemical identification of the crystal's
contents. Atom labels, restraints and a short modelled bond do not
independently identify the species.

The supplied PDBe map samples and PDB-REDO calculations were not reproduced
here. Model-phased density is not independent proof of the contested atom.
PDB-REDO reuses the diffraction data and can retain the starting chemical
assignment. Omit/polder maps can reduce particular biases, not guarantee
unbiased phases. Compare alternative models under consistent refinement
and held-out reflection choices, considering local fit and residuals as
well as global R-free. Do not turn an uncertainty on one R-free into a
universal cutoff for a paired difference.

Establish which lot/vial was assayed, MS-tested, soaked and crystallised,
with analytical limits stated. An impurity, mixed occupancy, oxidation
or an incorrect assignment remain hypotheses. The measured IC50 of
authentic 14a is not a crystal-condition binding constant, and
crystallisation does not select on affinity alone. A co-folder's preference
cannot adjudicate this identity conflict.

## 3. What would success on 4MEN tell you?

It would show that the recorded pipeline can reproduce aspects of this
known complex under a specified condition. 4MEN predates the reported
PDB structural cutoff, so agreement is recall-compatible. The exact
checkpoint's membership, related-complex exposure and affinity-label
inclusion have not been established; do not claim proven memorisation.

The primary condition is sequence + a cached full MSA, without a structural
template. Its first planned run has one seed and one sample. The no-MSA
run is only an optional smoke test. One run does not establish seed
stability, confidence calibration or a population success rate.

The proposed later case needed an additional check: locally, 32SC and
9FWX have identical compared protein/ligand coordinates. They represent
one public complex, first released in 2024, rather than an independent
2026 observation. Post-cutoff release alone also does not prove absence
from every training source or chemotype family. A credible later study
would deduplicate cases and audit protein/pocket and ligand relatives.

## 4. How do you compute ligand pose RMSD?

First declare the residue mapping and fit set. Here, the proper Kabsch fit
uses all 125 native Cα positions, author residues 44–168, excluding the
two expression-tag residues. Apply that transform to the ligand; never
refit the ligand itself.

Require the complete heavy-atom graph and enumerate graph isomorphisms
that preserve chemistry. Minimise the coordinate RMSD over those mappings
in the fixed protein frame. Pre-transforming coordinates does not disable
the fitting done by RDKit `GetBestRMS`; our implementation explicitly
enumerates matches and calculates distances without a ligand fit.

Missing atoms or incomplete mappings fail validation. A synthetic rigid
transform and symmetry permutation should score near zero, whereas moving
only the ligand should remain a large error. Those controls are tested.
The implemented pocket diagnostic is an 8 Å reference-pocket Cα error
under the same global fit; it is not a second fitted ligand RMSD.

Report contact precision/recall and the limited geometry screen alongside
RMSD. The preregistered ≤2 Å gate is operational. A higher RMSD does not
automatically identify a different binding site, and a low RMSD alone does
not establish physical validity.

## 5. How does the water network limit the interpretation?

The deposited 4MEN model places 25K N2 near Asn140 ND2 and N1 near HOH316.
The local distances reproduce the supplied approximate 2.84 Å and 2.82 Å
values. These are geometric observations, not a complete hydrogen-bond or
energetic analysis. A low model B-factor does not independently determine
occupancy, affinity or a binding mechanism.

This Boltz input/output setup does not predict the reference structural
waters. Declare their exclusion from the protein-residue contact metric
instead of crediting or penalising absent water-mediated contacts. Good
geometry may be consistent with learned effects of solvent, but it does
not reveal how the model represented those effects.

The publication's explicit retained waters and fixed holo receptor make
its docking setup different from sequence/MSA co-folding or a protein-only
template ablation. For 4MEO, the water-versus-oxide dispute still needs
experimental-data and sample evidence regardless of a predicted pose.
