# Deliverable 5 — Independent scientific assessment and prioritised uncertainties

Part A reviews the assumptions underlying the project as a whole (Task 4). Part B is the prioritised
uncertainty register. I have reviewed the *science*, not the Phase 1 code — the repository is private and I
did not access it.

---

## Part A — Independent assessment

### A1. Unsupported or over-reachable claims to remove from the project's framing

| Claim | Status | Correct statement |
|---|---|---|
| "4MEO is the BRD4(1) complex of compound 5" | **Unsupported.** The deposited ligand 25V is C18H16N2O2, the *N*-oxide. | "4MEO contains a ligand modelled as the quinoline *N*-oxide; the paper cites it as compound 5. The identity is unresolved." |
| "Compound 5 inhibits BRD4 with IC50 X" | **No such value exists.** Its IC50 exceeded its solubility. | "Compound 5 showed detectable inhibition; IC50 not determinable. The *N*-oxide 14a has IC50 42 µM." |
| "Vidler's docking achieved 0.30 Å RMSD, so docking works" | **Context-stripped.** | "Cross-docking into 3MXF, a holo receptor differing by 0.37 Å from the target structure, with six crystallographic waters retained." |
| "Boltz-2 reproducing 4MEN validates co-folding" | **Unsupported — leakage.** Released 2013-09-25; cutoff 2023-06-01. | "Validates the pipeline. Generalisation requires a post-cutoff arm." |
| "Boltz-2's affinity prediction for compound 3 is a prediction" | **Unsupported — label leakage.** This paper's data is ChEMBL document CHEMBL2434916. | "Plausibly in-distribution; and compound 3 has no measured IC50 to compare against." |
| "Deposited coordinates establish chemical identity" | **Category error.** | "A deposited model is an interpretation of electron density." |
| "The virtual screen's 4.2% hit rate proves the KAc-mimetic-first strategy" | **Confounded, n = 6.** | "Consistent with; the design cannot separate it from expert visual triage." |

### A2. Missing controls

1. **No held-out (post-cutoff) arm** — the single most important gap. Addressed as Arm D in Deliverable 3;
   32SC recommended.
2. **No seed-variance quantification.** A single Boltz-2 run is a draw from a stochastic sampler. Five
   seeds minimum, median and range reported.
3. **No negative control for the metrics.** Run a non-binder or a decoy through the identical pipeline. If
   RMSD, contact recovery and confidence do not separate it from 25K, the metrics are not measuring what
   you think.
4. **No confidence calibration.** Correlate `confidence_score` against actual RMSD. On a memorised target it
   will likely be high and flat — demonstrating that is a result.
5. **No MSA provenance control.** `--use_msa_server` is a moving target; save and reuse the `.a3m` or the
   work is not reproducible.
6. **No apo/template ablation separating fold recall from pose prediction** — the A/B/C arm structure in
   Deliverable 3 supplies this.
7. **For 4MEO: no occupancy-refinement control.** The competing models were never refined against each other.

### A3. Alternative interpretations the project should carry explicitly

- **4MEO identity:** four live hypotheses (pure *N*-oxide / compound 5 + short-contact water / mixed
  occupancy / oxidation during crystallisation or data collection). Deliverable 2 §5 ranks them; none is
  excluded.
- **4MEO's weak ligand density** is equally explained by a genuinely weak, mobile, partially occupied ligand
  as by a species-assignment error — compound 5's unmeasurable IC50 predicts exactly that.
- **Good pocket reproduction by a co-folder** is explained by memorisation at least as well as by learned
  physics.
- **My own corrected inference, recorded deliberately:** I initially read B(O) = 34.5 vs B(N1) = 13.6 Å² as
  evidence of partial occupancy. PDB-REDO's re-refinement of the same data removes most of that disparity
  (28.7 vs 21.7). The anomaly was substantially a refinement artefact. I downgraded the inference rather
  than keeping it because it was convenient.

### A4. Design improvements worth making

1. **Specify ligands by SMILES, not CCD code.** `ccd: 25K` hands the model the component template for a
   ligand unique to the entry being reproduced.
2. **Fix the residue-index convention once** (`seq[i-42]`, verified: 81 W, 88 D, 97 Y, 140 N) and use it
   everywhere. Off-by-one in residue mapping is the most common silent error in this kind of work.
3. **Hash the reference coordinates.** 4MEN is at revision 1.3 (2023-09-20); a future revision would
   silently change your RMSDs.
4. **Pre-register thresholds** before inspecting output.
5. **Report the four Vidler structures as a set.** 4MEP/4MEN/4MEQ/4MEO with their RSCC values is what makes
   4MEO's 0.872 interpretable — an isolated number is not evidence of an outlier.
6. **Exclude water-mediated contacts from recovery metrics, declared.**
7. **Report precision alongside contact recovery**, or a promiscuous prediction scores well.
8. **Treat 4MEO as forward-prediction only.** Co-folding both candidate species is hypothesis-generating; a
   model preference is not crystallographic evidence and must never be used to adjudicate identity.
9. **State provenance for every number** — measured / computed-by-me / predicted / inferred. The project's
   credibility rests more on this than on any single result.

### A5. What is genuinely strong about the project

The target choice is excellent: a small, well-characterised system with four sibling co-crystal structures
from one paper, one of which contains a real, documented chemical-identity ambiguity. That ambiguity is a
better interview asset than a clean reproduction, because it demonstrates the thing that is actually hard —
distinguishing what the data establish from what the model asserts. Keep it central.

---

## Part B — Prioritised uncertainty register

Priority = (effect on conclusions) × (tractability). **P1 = resolve before claiming anything.**

### P1 — Resolve first

1. **Chemical identity of the species in the 4MEO crystal.** *Unresolved.* Blocks every claim about compound
   5's binding mode. Needs: polder/omit map at the oxide oxygen from `4MEO-sf.cif`; competing-model
   refinement with occupancy refined, scored on ΔR-free against σ(R-free) = 0.0077; LC-MS of ChemBridge
   50866471 against an authentic 14a standard. *Blocked tonight: requires gemmi/CCP4/Phenix.*
2. **Whether any Boltz-2 result on 4MEN can support a capability claim.** *Resolved as "no".* Mitigation is
   Arm D (32SC). Must be settled in the framing before results are generated, not after.
3. **Residue-numbering and construct convention.** *Resolvable now, and cheap.* `seq[i-42]`; decide on
   including the `SM` cloning remnant and keep it fixed. An error here invalidates every contact metric.

### P2 — Resolve before quantitative claims

4. **Boltz-2 compute envelope for this system.** *Estimated only, not measured.* One timing run settles it.
   Everything I wrote about VRAM and wall-clock is inference.
5. **Hosted-API cost.** *Unknown.* `boltz_estimate_structure_and_binding` returns an exact figure for free,
   but the connector returned **401 Bearer token authentication failed** — reconnect it under
   Customize → Connectors. I have not guessed a price.
6. **Seed-to-seed variance of Boltz-2 on this complex.** *Unknown.* Five seeds per arm.
7. **Whether 32SC's ligand is chemotype-novel** relative to pre-cutoff BET inhibitors. Compute max Tanimoto;
   a post-cutoff entry with a familiar chemotype still leaks.
8. **Affinity-head interpretability for this series.** Labels are in ChEMBL (CHEMBL2434916) *and* compound 3
   has no measured IC50. Report the number only with both caveats, or omit it.

### P3 — Known limitations to state, not resolve

9. **No structural waters in co-folding output.** Structural, not fixable. Affects the KAc network directly.
10. **Static prediction vs a single 100 K crystal snapshot** in P2₁2₁2₁ with DMSO and ethylene glycol in the
    lattice; no crystal-packing term.
11. **Confidence calibration on a memorised target** is expected to be uninformative; measure and report it.
12. **4MEO's weak ligand density** has at least two sufficient explanations (§A3); do not assert one.

### P4 — Access limitations, reported not worked around

13. **Supporting Information** (Figures S1–S3, crystallographic Table S1) not retrievable open-access; the
    substructure library and original docking poses are not independently verifiable here. Reachable with an
    institutional ACS login.
14. **Figure 11** — the published density figure central to the 4MEO question — not downloadable. I analysed
    the deposited maps instead, which is better evidence, but the published figure could not be compared
    against the final coordinates.
15. **GitHub repository is private** and no authorised integration is connected. I accessed nothing, assumed
    nothing about Phase 1, and modified nothing.
16. **Package installation unavailable this session** by your instruction. Consequences: no omit map, no
    occupancy refinement, no RDKit symmetry-RMSD implementation, no Boltz-2 run. All input data
    (`4MEO-sf.cif`, both map pairs, PDB-REDO model, coordinates, CCDs) is cached in `data/` for when the
    environment is available.

### Not assessed

Crystallographic twinning/pseudosymmetry in either entry; anisotropy of the diffraction data; alternative
space-group assignments; the radiation-damage hypothesis for 4MEO (listed as a live explanation in
Deliverable 2 §5 but with no evidence sought either way).
