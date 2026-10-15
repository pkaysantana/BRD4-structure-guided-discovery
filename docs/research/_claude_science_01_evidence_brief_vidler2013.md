# Deliverable 1 — Scientific evidence brief: Vidler et al. (2013), BRD4 structure-based virtual screening

**Primary source.** Vidler LR, Filippakopoulos P, Fedorov O, Picaud S, Martin S, Tomsett M, Woodward H,
Brown N, Knapp S, Hoelder S. *Discovery of Novel Small-Molecule Inhibitors of BRD4 Using Structure-Based
Virtual Screening.* J Med Chem 2013;56(20):8073–8088. DOI 10.1021/jm4011302; PMID 24090311; PMC3807807.
Full text retrieved from PubMed Central (ACS AuthorChoice, CC-BY) and read in full for this brief.
The Supporting Information PDF (`jm4011302_si_001.pdf`, containing Figures S1–S3 and crystallographic
Table S1) was **not** retrievable through the open-access route — see "Access limitations" below.

Evidence labels used throughout: **[M]** measured/reported experimentally, **[C]** computed by me in this
session from public primary data, **[I]** inferred — my interpretation, not a claim of the source.

---

## 1. Why bromodomains, and BRD4 in particular, were drug-discovery targets

**[M]** Bromodomains are acetyl-lysine (KAc) reader modules — protein–protein interaction domains, not
enzymes. The paper frames the BET subfamily (BRD2, BRD3, BRD4, BRDT), each with two bromodomains of high
sequence similarity and highly similar binding sites, as the tractable subset.

**[M]** The structural rationale given is specific and geometric: the KAc site is a central cavity formed by
a left-handed four-helix bundle flanked by the ZA and BC loops; it is primarily hydrophobic, and the polar
recognition is carried by a family-conserved asparagine (Asn140 in BRD4-BD1) plus a structurally conserved
water. Occupying two adjacent regions — the "hydrophobic shelf" and the "ZA channel" — is what converts a
KAc mimetic into a nanomolar, selective inhibitor; (+)-JQ1 (PDB 3MXF) is the exemplar.

**[M]** The disease rationale cited is target-validation genetics rather than chemistry: RNAi/shRNA screens
identifying BRD4 as a therapeutic target in AML (Zuber 2011, Nature 478:524; Blobel 2011, Cancer Cell
20:287), plus lung adenocarcinoma, glioblastoma, neuroblastoma, and inflammation (Nicodeme 2010, Nature
468:1119).

**[I]** Worth being precise in interview: the 2013 justification for BRD4 is *chemically* "this pocket is
small, enclosed and has a defined H-bond pattern, therefore druggable", and *biologically* "knockdown kills
AML cells". Those are two independent arguments, and only the first is what the virtual screen actually
exploits. The authors' own prior work (Vidler 2012, J Med Chem 55:7346) is the druggability analysis that
underpins the first argument.

---

## 2. How the virtual screening campaign actually worked

The design intent was to make docking tractable by pre-enriching for KAc mimetics, rather than docking
everything. **[M]** Two branches fed a common docking-and-triage funnel:

1. **Literature Substructures branch.** KAc-mimetic substructures were extracted from published
   bromodomain inhibitors and deliberately mutated into pharmacophore-preserving relatives (the worked
   example is 1,2,4-triazole → isoxazole, starting from (+)-JQ1). These became substructure queries.
2. **Similarity Searching branch.** Three structure-based pharmacophore searches, two shape searches
   (vROCS/ROCS) and one ECFP4 2D-fingerprint search (Pipeline Pilot), all templated on the bound
   conformation of (+)-JQ1 from 3MXF, with the receptor as excluded volume. Searches ran over the ICR
   Vendor Collection, 2.4 million commercially available compounds, with ≤100 conformers per molecule
   (OMEGA) and single protonation states (FILTER).

**[M]** Reported pharmacophore hit counts are a clean illustration of constraint relaxation: the five-point
pharmacophore (two triazole acceptors, a hydrophobic at the pocket base, plus ZA-channel and shelf
hydrophobics) matched 6,540 ligands; relaxing to one of the two acceptors gave 268,474; removing the shelf
point instead gave 41,890. Shape and 2D searches each carried their top 100,000 forward.

**[M]** Docking: Glide against 3MXF, prepared with Protonate3D (MOE) and the Maestro protein-preparation
wizard, 20 Å box around (+)-JQ1, ligands prepared with LigPrep/Epik. Six crystallographic waters (UIDs 9,
12, 15, 23, 33, 209) were retained as part of the receptor and treated as rigid. Sets >5,000 compounds were
funnelled HTVS → top 5,000 → SP; smaller sets went straight to SP.

**[M]** Triage was explicitly manual and criterion-driven, motivated by the docking false-positive
literature (Ferreira 2010, J Med Chem 53:4891): reject high-energy conformations (validated against
observed CSD dihedral distributions), require heteroatoms to be satisfied, reject clashes or unfillable
vacancies, reject PAINS/reactive groups (Baell & Holloway 2010). At least the top 500 poses per docking run
were inspected by eye, yielding 150 compounds.

**[M]** **The control is the most important methodological feature of the paper.** A fully automated
ligand-based arm — the same 2D-fingerprint and shape searches on (+)-JQ1, top 50 unique Murcko scaffolds
each, no docking, no inspection — selected 100 compounds with zero overlap with the virtual-screen set.

**[M]** Outcome: 250 purchased, 10 DMSO-insoluble, 240 screened by AlphaScreen at up to 250 µM. Six
inhibitors, all from the virtual-screen arm, none from the control arm.

| Arm | Tested | Hits | Hit rate |
|---|---|---|---|
| All compounds | 240 | 6 | 2.5% |
| Virtual screen (both branches) | 143 | 6 | 4.2% |
| — Literature substructures | 76 | 2 | 2.6% |
| — Similarity searching | 67 | 4 | 6.0% |
| Control (ligand-based only) | 97 | 0 | 0% |

**[I]** Critical reading of the control: 0/97 versus 6/143 is a real signal, but with six events in total the
confidence interval on any per-branch hit rate is very wide, and the two arms differ in *several* ways at
once (docking, human inspection, substructure enumeration, scaffold-diversity selection). The paper
attributes the difference to the KAc-mimetic-first design; the experiment as run cannot separate that from
the effect of expert visual triage alone. The authors' own observation that the Literature Substructures
branch *under*-performed the Similarity branch (2.6% vs 6.0%), contrary to their prior expectation,
is honestly reported and further argues against a simple single-cause explanation.

---

## 3. How acetyl-lysine mimetics were identified and evaluated

**[M]** Identification was substructure-driven: a chemotype library assembled from published inhibitors,
"intuitive design", and similarity searching, converted into substructure queries and run against
eMolecules to retrieve every commercially available instance. The authors contrast this with similarity
search, which they note misses chemotype representatives because of tautomer, protomer and conformer
differences in vendor libraries.

**[M]** Evaluation used three nested tests: (i) biochemical inhibition (AlphaScreen, His-BRD4(1) 250 nM
versus a tetra-acetylated biotinylated H4 peptide, 24-point 1:2 dilutions from 150 µM); (ii) co-crystal
structure determination; (iii) SAR by purchase and synthesis of analogues.

**[M]** Four previously unprecedented KAc mimetics resulted. Potency of the primary hits was modest —
compound **1** (chloropyridone) IC50 4.7 µM, compound **2** (isoxazole) 80.9 µM; the other four hits had
IC50 above their solubility limit and could not be quantified. SAR recovered the gap: triazolopyrimidine
**8b** IC50 24 µM (LE 0.33), **8c** improved over **3** and **4**, pyrrolopyrimidones **15a/15b** reached
LE > 0.3.

**[M]** Two mechanistic results stand out as structurally informative rather than merely potent:
- In compound **1**, a **chlorine** occupies the base of the KAc pocket where the KAc methyl normally sits —
  described as the first such example, and proposed as a transferable replacement for the usual methyl.
- The **triazolopyrimidine** of **3** and **4** reproduces the acetyl-carbonyl interaction through *two ring
  nitrogens* substituting for the two carbonyl lone pairs, and the 2-amino group of **4** completes the
  acceptor–acceptor–donor motif seen in dihydroquinazolinones.

**[C]** I verified the deposited chemistry of the primary triazolopyrimidine: CCD **25K** in 4MEN is
C15H17N5, `Cc1ccc(cc1)CN(C)c2cc(nc3n2ncn3)C`, 20 heavy atoms, no stereocentres — matching compound **3**,
*N*,5-dimethyl-*N*-(4-methylbenzyl)-[1,2,4]triazolo[1,5-a]pyrimidin-7-amine, as reported.

---

## 4. What the experimental structures establish about ligand recognition

**[M]** Four co-crystal structures were deposited, all BRD4-BD1, all in-house Rigaku FRE rotating-anode data
at 1.52 Å wavelength, MOSFLM/SCALA, molecular replacement from 2OSS with PHASER, built in ARP/wARP + COOT,
refined in REFMAC5 with TLS: **4MEP** (compound 1), **4MEN** (compound 3), **4MEQ** (compound 4),
**4MEO** (compound 5).

**[C]** Independent validation retrieved this session (RCSB ligand quality scores, PDBe validation API,
PDB-REDO) — note that the ligand of 4MEO is the clear outlier on density fit:

| Entry | Compound | Resolution | R-work/R-free | Ligand | RSCC | RSR | Fit percentile |
|---|---|---|---|---|---|---|---|
| 4MEP | 1 | 1.85 Å | — | 24Y | 0.949 | 0.077 | 0.741 |
| 4MEN | 3 | 1.81 Å | 0.163 / 0.219 | 25K | 0.953 | 0.055 | 0.838 |
| 4MEQ | 4 | 1.77 Å | — | 25O | 0.965 | 0.049 | 0.890 |
| 4MEO | 5 (see Deliverable 2) | 1.72 Å | 0.1653 / 0.2126 | **25V** | **0.872** | **0.097** | **0.506** |

**[C]** Construct, measured from the deposition: a single chain, UniProt O60885 (BRD4_HUMAN), author
residues 42–168, 127 residues, **no chain breaks**, expressed in *E. coli*. Author numbering equals UniProt
numbering from residue 44; residues 42–43 (Ser-Met) are a cloning remnant. Space group P2₁2₁2₁, cell
37.519 × 44.229 × 78.020 Å (4MEN), monomeric assembly.

**[C]** Binding-site contacts I computed from 4MEN coordinates (heavy-atom distances; residues within
4.2 Å of 25K): Trp81, Pro82, Phe83, Val87, Leu92, Leu94, Tyr97, Tyr139, Asn140, Asp145, Ile146, Met149 —
i.e. the WPF shelf (Trp81-Pro82-Phe83), the gatekeeper Ile146, and the conserved Asn140/Tyr97 pair. The two
polar anchors are **25K N2···Asn140 ND2 = 2.84 Å** and **25K N1···HOH316 = 2.82 Å**, the latter being the
structurally conserved water (B = 6.2 Å², the best-ordered water in the site).

**[M]** What the structures establish, stated conservatively:
- The triazolopyrimidine *is* a genuine KAc mimetic — this is the hypothesis the structures were obtained to
  test, and they confirm it. Direct evidence.
- Compound **5** binds in an orientation where the **acetyl group is not the KAc mimetic**; instead it
  contacts the **backbone NH of Asp88**, an interaction the authors state had not been seen in any published
  bromodomain structure. **[C]** I confirm the geometry in the deposition: acetamide O1···Asp88 N = 2.67 Å.
- Docking reproduced the crystallographic pose well for **3** and **4** (rmsd ≤ 0.30 Å, "within experimental
  error") and acceptably for **1** and **5** (rmsd < 2.0 Å), with discrepancies attributed to bridging waters
  the docking could not model.

**[I]** The ≤0.30 Å docking result should be quoted with its caveat, which the authors themselves supply:
the protein barely moved between the 3MXF template and the triazolopyrimidine-bound structures (0.37 Å
rmsd). This is near-ideal *cross-docking into an essentially pre-organised holo receptor*, not a prospective
test of induced fit. It is a strong result for pose prediction and a weak one for pocket prediction.

---

## 5. Supported versus inferred, as the paper itself stands

**Experimentally supported:** the six hit identities and their AlphaScreen inhibition; IC50 values for
compounds **1**, **2**, **8b–e**, **14a**, **15a/b**; the four co-crystal binding modes at 1.72–1.85 Å;
the 0/97 control-arm result; the identity and purity of synthesised compounds (HRMS, ¹H NMR); the
triazolopyrimidine as a KAc mimetic; the Asp88 backbone-NH contact of compound **5**.

**Inferred or model-dependent in the paper:**
- The *causal* attribution of the hit-rate difference to the KAc-mimetic-first strategy (confounded, as
  above; n = 6 events).
- "Four unprecedented KAc mimetics" — a literature-coverage claim, true as of 2013, not an experiment.
- Potency rationalisations ("the amino group completes the acceptor–acceptor–donor motif"; "8c improves over
  4 because the benzylamino reaches the hydrophobic shelf") are structure-guided interpretations; the
  matched-pair comparisons support them but do not isolate the mechanism.
- Compounds whose IC50 exceeded solubility (including **5**) have *inhibition detected*, not potency
  measured. Any downstream use of **5** as a potency datapoint is unsupported.
- The chemical identity of the species in the 4MEO crystal — the subject of Deliverable 2.

---

## 6. Relationship to modern protein–ligand co-folding

**[I] Framing for interview.** Vidler 2013 is a *pose-prediction-into-a-known-receptor* campaign. Modern
co-folders (Boltz-2, Chai-1, AlphaFold3, OpenFold3) predict the receptor and the ligand *jointly* from
sequence plus a chemical graph. The comparison is therefore not like-for-like in three specific ways:

1. **Receptor provenance.** Glide was given 3MXF, a holo BRD4-BD1 structure with (+)-JQ1 bound, including
   six hand-selected ordered waters. A co-folder is given a sequence. Reproducing the 4MEN pose from
   sequence is a strictly harder task than the docking the paper reports — and a correspondingly weaker
   claim if it succeeds, because of point 3.
2. **Explicit water.** The conserved water network is load-bearing for bromodomain recognition; the paper
   both retains waters in the receptor and attributes its docking discrepancies to waters it could not
   model. Current co-folders do not place structural waters, so any water-bridged contact is, at best,
   implicitly absorbed into the learned geometry. **[C]** In 4MEN the single most informative solvent atom
   (HOH316, 2.82 Å to the ligand) is exactly the kind of feature a co-folder cannot output.
3. **Training-set overlap.** **[C]** 4MEN and 4MEO were released 2013-09-25. Boltz-2's documented structural
   training cutoff is 2023-06-01, and Boltz-1's was 2021-09-30; both entries fall well inside either. Any
   agreement between a Boltz-2 prediction and 4MEN is *recall-compatible* and cannot by itself demonstrate
   generalisation. **[C]** The leakage is not only structural: this paper's assay data is in ChEMBL as
   document CHEMBL2434916, carrying 20 activities (10 IC50s) against BRD4-BD1 (CHEMBL1163125) — and Boltz-2's
   affinity head is trained on ChEMBL-derived data. For this compound series, both the structure labels and
   the affinity labels are plausibly in-distribution.

**[I]** The honest scientific position: a 4MEN co-folding experiment is a **positive control that validates
the harness**, not evidence about model capability. It is still worth running — you need to know your
pipeline, metrics and alignment are correct on a case where the answer is known — but the claim must be
scoped to that. This is the single most important framing point in the project, and Deliverable 3 is built
around it.

---

## Access limitations (reported, not worked around)

- **Supporting Information** (Figures S1–S3; crystallographic data-collection/refinement Table S1) could not
  be retrieved open-access. Consequences: the full substructure library (S1) and the original docking poses
  (S2) are not independently verifiable here, and the authors' own data-collection statistics were replaced
  by values I read from the depositions and PDBe. The ACS page is reachable with an institutional login.
- **Figures 1–15** of the main text, including **Figure 11** (the 2Fo−Fc density around compound 5 that is
  central to Deliverable 2), were not downloadable by the full-text route. I have therefore analysed the
  *deposited maps themselves* rather than the published figure.
- **GitHub repository** `pkaysantana/BRD4-structure-guided-discovery` is private and no authorised GitHub
  integration is connected to this session. I have not accessed it and have assumed nothing about the Phase 1
  codebase; nothing in these deliverables duplicates or modifies it.
