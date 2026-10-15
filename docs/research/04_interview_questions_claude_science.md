# Deliverable 4 — Five technical questions you should be able to answer

Framed as an interviewer would ask them, each with the evidence you need and the trap to avoid.
Every number quoted here is traceable to Deliverables 1–3 and `measurements_4MEN_4MEO.csv`.

---

### Q1. "Your virtual screen reported a 4.2% hit rate against 0% for the control. How strong is that evidence,
and what would you have done differently?"

**What to say.** It is a real signal — 6/143 against 0/97, and all six actives came from the structure-based
arm. But it rests on **six events in total**, so any per-branch rate (2.6% literature-substructure vs 6.0%
similarity) has a confidence interval wide enough to swallow the difference; the authors themselves note the
branch ordering came out opposite to their expectation. More importantly the two arms differ in **several
ways simultaneously** — KAc-mimetic substructure enumeration, docking, *and* expert visual triage of ≥500
poses per run. The design cannot attribute the result to the KAc-mimetic-first hypothesis specifically;
human inspection alone is a sufficient alternative explanation.

**What you'd change:** add a third arm — docking plus automated scoring, no visual inspection — to isolate the
contribution of human triage. And pre-register the hit definition, since four of the six "hits" had IC50
above their solubility limit and are inhibition-detected rather than potency-measured.

**Trap:** don't call 0/97 "proof the method works". Say what it is: consistent with, underpowered for
attribution.

---

### Q2. "PDB 4MEO is cited as compound 5, but the deposited ligand is a quinoline *N*-oxide. Walk me through
how you'd resolve that."

**What to say.** Establish the facts first, in order. The deposited component **25V** is C18H16N2O2 with an
explicit `O–N1` single bond in the CCD and `[n+]…[O-]` in the SMILES — it is chemically identical to the
paper's synthesised *N*-oxide **14a**, and different from compound **5** (C18H16N2O, HRMS found 277.1333).
In the coordinates N1–O is **1.358 Å**: a covalent bond, with no water modelled. So the deposition models
the oxide, unambiguously, while the paper's text reports compound 5 with a bridging water at 2.14 Å and
explicitly records that **mass spectrometry did not detect the *N*-oxide**.

Then say what the crystallography can and cannot settle. 25V is the worst-fitting ligand of the four Vidler
depositions — RSCC 0.872 and RSR 0.097, against 0.953/0.055 for 25K in 4MEN, fit percentile 0.51 vs 0.84.
In the PDBe maps the oxide oxygen sits in only ~1.6σ of 2Fo−Fc density against 3.44σ at the nitrogen it is
bonded to, and the density along the N1→O vector is a broad shoulder rather than a peak at either 1.358 Å
or 2.14 Å. **But** those maps are phased on the deposited model, so they are biased toward showing an oxygen
there, and PDB-REDO's independent re-refinement keeps the *N*-oxide at full occupancy and improves R-free
from 0.2126 to 0.2072. The public evidence does not discriminate.

**The resolution:** a polder/omit map at that atom, plus competing-model refinement with occupancy refined
(oxide / quinoline+water / mixture), scored on ΔR-free against σ(R-free) = 0.0077 — and, decisively, LC-MS
of the actual commercial sample against an authentic 14a standard. The chemically coherent hypothesis is
that a minor *N*-oxide impurity was preferentially selected into the crystal, because 14a has a measured
IC50 of 42 µM while 5 has none — crystallisation selects on affinity, not abundance.

**Trap:** do not say "the PDB says it's an *N*-oxide, so it's an *N*-oxide." A deposited model is an
interpretation of density, not a measurement of chemical identity. Equally, do not assert the impurity
hypothesis as fact — it is the best explanation, not a result.

---

### Q3. "You want to test Boltz-2 on 4MEN. What does success actually tell you?"

**What to say.** By itself, almost nothing about capability. 4MEN was released **2013-09-25**; Boltz-2's
documented structural training cutoff is **2023-06-01**. The entry is in the training set, BRD4-BD1 is one
of the most heavily represented pockets in the PDB, and the affinity labels leak too — this paper's assay
data is in ChEMBL as document **CHEMBL2434916** (20 activities, 10 IC50s against BRD4-BD1), and Boltz-2's
affinity head trains on ChEMBL. So a good result validates **my harness** — sequence handling, residue
numbering, protein superposition, symmetry-corrected RMSD, contact definitions — which is genuinely worth
establishing, and that is the claim I would make.

The informative experiment is the held-out arm: the same pipeline on a BRD4-BD1 complex released after
2023-06-01. There are 62 such X-ray entries better than 2.0 Å; **32SC** (2026-08-05, 1.25 Å, ligand A1IGX,
RSCC 0.933) is a clean single-ligand case. And even then I'd report the new ligand's maximum Tanimoto to
known BET chemotypes, because a post-cutoff *entry* with a pre-cutoff *chemotype* still leaks.

**Trap:** never present a sub-2 Å RMSD on a pre-cutoff structure as evidence of generalisation. Interviewers
at a company doing this seriously will be listening for exactly that.

---

### Q4. "How do you compute ligand RMSD properly, and why does it so often get reported wrong?"

**What to say.** Three failure modes, all common:

1. **Atom-order RMSD on a symmetric molecule.** 25K's *p*-tolyl ring has a topological two-fold, so a naive
   RMSD can report several Å of error for a pose that is correct. You must minimise over the molecular
   graph's automorphisms — RDKit `GetBestRMS` or `spyrmsd`.
2. **Aligning the ligand instead of the protein.** The protein is superposed and the resulting transform is
   *applied* to the ligand; the ligand is never fitted to itself. `GetBestRMS` aligns by default, which
   silently turns a pose-accuracy metric into a conformer-similarity metric. Pass pre-transformed
   coordinates and suppress further alignment.
3. **Quoting an RMSD without naming its superposition.** Global Cα versus binding-site-only superposition
   give different numbers; a global fit can look excellent while the pocket is subtly wrong, and the ligand
   RMSD depends on the pocket. Report both.

Add that RMSD alone is insufficient: a PoseBusters-style validity check (bond geometry, planarity, internal
and protein–ligand clashes) should accompany it, because a 1.5 Å RMSD with distorted geometry is not a
success. Threshold convention: < 2 Å reproduced, 2–5 Å right site/wrong pose, > 5 Å wrong site —
pre-registered, with the distribution over seeds reported, not the best run.

---

### Q5. "Bromodomain recognition depends on a conserved water network. What does that mean for co-folding
methods?"

**What to say.** It is a structural limitation, and it should be stated rather than scored around. In 4MEN
the ligand makes two polar anchors: **25K N2···Asn140 ND2 = 2.84 Å**, the acetyl-carbonyl mimic, and
**25K N1···HOH316 = 2.82 Å** — a water with B = 6.2 Å², the best-ordered solvent atom in the site, part of
the structurally conserved network that also contacts Tyr97. Vidler's own docking **retained six
crystallographic waters as part of the receptor** (UIDs 9, 12, 15, 23, 33, 209) and the authors attribute
their residual pose discrepancies for compounds 1 and 5 precisely to bridging waters the docking could not
model.

Boltz-2 does not output structural waters. So a water-mediated contact cannot be recovered, and must be
excluded from the contact-recovery denominator with that exclusion declared — otherwise you are penalising
the model for something it cannot produce, or worse, silently crediting it. The deeper point: the model may
have learned the *geometric consequence* of the water network implicitly, since it fits thousands of
bromodomain structures, which is why pocket geometry can look right while the physical mechanism is absent
from the output. That gap is exactly where a co-folding result stops being mechanistically interpretable.

**Bonus credit:** this is also why 4MEO is hard — the contested density is a putative *bridging water versus
a covalent oxide*, and a method that cannot represent water has no way to adjudicate it.
