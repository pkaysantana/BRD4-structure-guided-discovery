# Phase 1 report — Experimental baseline for BRD4-BD1 ligand complexes

**Status:** analysis of deposited experimental models only. No predicted
structures, no docking, no prospective screening.

**How numbers get here.** Every computed value in this report sits inside a
`<!-- BEGIN GENERATED -->` block that is re-rendered from
`results/summary.json` and `data/provenance/provenance.json` by
`python -m brd4_analysis.cli report` (run automatically by `cli all`).
`tests/test_regression.py` fails if the committed report differs from a
fresh render. Values quoted from the publication or the wwPDB validation
report are external citations, labelled as such, and are not reproduced by
this pipeline.

**Reference:** Vidler LR, Filippakopoulos P, Fedorov O, Picaud S, Martin S,
Tomsett M, Woodward H, Brown N, Knapp S, Hoelder S. *Discovery of Novel
Small-Molecule Inhibitors of BRD4 Using Structure-Based Virtual Screening.*
J. Med. Chem. 2013, 56(20), 8073–8088. doi:10.1021/jm4011302;
PMID 24090311; PMCID PMC3807807.

## 0. Study context (from the publication)

Vidler et al. ran a **structure-based virtual screen**. They used
substructure queries to assemble commercially available compounds carrying
putative acetyl-lysine (KAc) mimetics, then docked that set against BRD4
BD1 to pick compounds. They tested 143 compounds experimentally and
reported six hits, four of them with KAc mimetics new at the time. Crystal
structures were solved for four hits, and potency was improved by
synthesising and purchasing analogues. In the publication's numbering,
**4MEN is BRD4(1) with compound 3** (a triazolopyrimidine hit) and **4MEO
is BRD4(1) with compound 5** (a 2-methylquinoline hit). This is not a
fragment-screening campaign.

## 1. Data and provenance

<!-- BEGIN GENERATED: provenance -->
| Entry | Ligand (CCD) | Method | Resolution (Å) | Deposited | Latest REVDAT | Retrieved | PDB file SHA-256 |
|---|---|---|---|---|---|---|---|
| 4MEN | 25K | X-RAY DIFFRACTION | 1.81 | 2013-08-27 | #4, 20-SEP-23 | 2026-10-14 | `a1b21bdc30365fb0ca905a6f9f00cc5955e3b382e6b726c02e4a2162be55b6a0` |
| 4MEO | 25V | X-RAY DIFFRACTION | 1.72 | 2013-08-27 | #4, 20-SEP-23 | 2026-10-14 | `16f7df00cb65a5658c5e3bbb97458327ce73a27993b7304f035753215fa1b3de` |

Input files recorded in the manifest: 6.
<!-- END GENERATED: provenance -->

Source URLs, per-file UTC retrieval timestamps and checksums are in
`data/provenance/provenance.json`. Manifest paths are relative to `data/`.
Raw files are committed unmodified (`.gitattributes` marks `data/**` as
`-text`). `cli fetch` does not overwrite an existing file. It re-verifies
the file against its recorded SHA-256 and refuses on mismatch unless
`--force` is given.

## 2. Ligand identity: what is and is not established

**Separating the ligand from other HET groups.** The bound ligand is
identified by its CCD component ID: 25K in 4MEN, 25V in 4MEO. These IDs
were taken from RCSB non-polymer entity metadata. Ethylene glycol (EDO),
dimethyl sulfoxide (DMS) and water are classified as crystallisation
additives or solvent.

**Graph consistency with the CCD definition.** For both entries, assigning
bond orders from the CCD template to the deposited HETATM records gives
exactly the CCD canonical isomeric SMILES
(`identity_check.level = "smiles_match"`). This shows that the deposited
heavy-atom connectivity agrees with the CCD definition of the component.
It **does not** show experimental protonation state, tautomeric form or
hydrogen positions. The template supplies the bond orders, formal charges
and hydrogen counts, and X-ray data at this resolution carry no hydrogen
information. Neither ligand has a stereocentre, so there is no
stereochemistry to verify.

**Unresolved identity question for 4MEO (25V).** The deposited component,
25V, is the quinoline **N-oxide**: CCD name
*N-[3-(2-methyl-1-oxidoquinolin-4-yl)phenyl]acetamide*, formula
C18H16N2O2. The publication, however, describes the 4MEO complex as
compound 5, the parent quinoline without the N-oxide oxygen. Its account,
quoted from the paper and not re-derived here, is as follows. Electron
density was seen very close to the quinoline nitrogen and was interpreted
as a bridging water at an unusually short 2.14 Å. The authors speculated
that the active species might be the N-oxide, but mass spectrometry did
not detect it. They then synthesised the N-oxide (compound 14a), which was
more potent (reported IC50 42 μM). The deposited model therefore puts a
covalently bonded N-oxide oxygen (atom `O`) where the paper's narrative
describes density attributed to water. **Whether the crystal contained
compound 5, the N-oxide, or a mixture cannot be settled from the
coordinates.** Answering it would need the structure factors and electron
density, or the deposition correspondence. This baseline reports 25V as
deposited and carries this uncertainty into every 4MEO statement below,
including the 25V descriptors, which describe the N-oxide.

## 3. Ligand properties (CCD reference chemistry, deposited components)

<!-- BEGIN GENERATED: descriptors -->
| Descriptor | 25K (4MEN) | 25V (4MEO) |
|---|---|---|
| Formula | C15H17N5 | C18H16N2O2 |
| Isomeric SMILES | `Cc1ccc(CN(C)c2cc(C)nc3ncnn23)cc1` | `CC(=O)Nc1cccc(-c2cc(C)[n+]([O-])c3ccccc23)c1` |
| InChIKey | JAEVUMCKPBRDAJ-UHFFFAOYSA-N | YSGOTUAMXVUKOD-UHFFFAOYSA-N |
| MW (Da) | 267.34 | 292.34 |
| cLogP (Crippen) | 2.38 | 3.41 |
| HBD (RDKit CalcNumHBD) | 0 | 1 |
| HBA (RDKit CalcNumHBA) | 4 | 2 |
| HBD, Lipinski N-H + O-H | 0 | 1 |
| HBA, Lipinski N + O | 5 | 4 |
| TPSA (Å²) | 46.32 | 56.04 |
| Rotatable bonds | 3 | 2 |
| Heavy atoms | 20 | 22 |
| Aromatic rings | 3 | 3 |
| Fraction Csp3 | 0.27 | 0.11 |
| Formal charge | 0 | 0 |
<!-- END GENERATED: descriptors -->

**Descriptor definitions (RDKit, pinned version in `summary.json`).**
MW uses average isotopic masses. cLogP is Wildman–Crippen. TPSA is
Ertl (N and O only). Rotatable bonds use the RDKit default definition. Two
H-bond count conventions are reported:

- **`CalcNumHBA`/`CalcNumHBD`** (behind `Lipinski.NumHAcceptors`/
  `NumHDonors`) are SMARTS-based and **version-dependent**. Between RDKit
  2025.03 and 2025.09, `CalcNumHBA` stopped counting N-substituted
  aromatic (pyrrole-type) nitrogens. In our tests this changed caffeine
  from 6 to 3 and 25K from 5 to 4; 25V is 2 in every version tested
  (2023.09.6, 2024.09.6, 2025.03.6, 2025.09.6, 2026.03.6, 2026.09.1).
  MW, cLogP, TPSA, rotatable bonds and the Morgan Tanimoto were identical
  across all of those versions.
- **Lipinski N+O / NH+OH counts** follow the original rule-of-five
  definition and do not depend on the version.

**Fingerprint similarity.**

<!-- BEGIN GENERATED: similarity -->
Tanimoto (4MEN vs 4MEO) = **0.0833** (6 common / 72 union on-bits; 38 and 40 on-bits per ligand). Fingerprint: Morgan radius 2, 2048 bits, heavy-atom graph (explicit H removed).
<!-- END GENERATED: similarity -->

The fingerprint is computed on the heavy-atom graph. An earlier draft of
this report quoted 0.154. That value came from a superseded pipeline state
in which the CCD SDF was loaded **with explicit hydrogens**, which adds
H-environment bits. It never matched the committed `summary.json` and has
been withdrawn. The similarity is low, so by this fingerprint the two
ligands are chemically distinct. That is consistent with the publication
presenting them as different KAc-mimetic chemotypes.

## 4. Experimentally observed contacts (heavy-atom geometry only)

A residue is "in contact" at threshold *t* iff its minimum heavy-atom
distance to the ligand is ≤ *t*. Only polymer residues are counted. No
hydrogen-bond or interaction-type classification is applied.

<!-- BEGIN GENERATED: thresholds -->
| Threshold | 4MEN residues | 4MEO residues |
|---|---|---|
| ≤ 3.5 Å | 5 | 5 |
| ≤ 4.0 Å | 11 | 11 |
| ≤ 5.0 Å | 13 | 14 |
<!-- END GENERATED: thresholds -->

Shortest protein–ligand heavy-atom pairs (`contact_pairs_*.csv`):

<!-- BEGIN GENERATED: closest-pairs -->
| Entry | Protein atom | Ligand atom | Distance (Å) |
|---|---|---|---|
| 4MEN | A:ASN140:ND2 | 25K:N2 | 2.844 |
| 4MEN | A:PHE83:CE1 | 25K:C10 | 3.353 |
| 4MEN | A:ASN140:OD1 | 25K:C12 | 3.373 |
| 4MEN | A:PHE83:CD1 | 25K:C10 | 3.378 |
| 4MEN | A:TRP81:CH2 | 25K:C2 | 3.420 |
| 4MEN | A:ASN140:OD1 | 25K:N2 | 3.427 |
| 4MEO | A:ASP88:N | 25V:O1 | 2.670 |
| 4MEO | A:ASN140:ND2 | 25V:O | 2.881 |
| 4MEO | A:LEU92:CD1 | 25V:O1 | 2.910 |
| 4MEO | A:LEU92:CD1 | 25V:C3 | 2.945 |
| 4MEO | A:ASN140:OD1 | 25V:C10 | 3.289 |
| 4MEO | A:LEU92:CD1 | 25V:C4 | 3.335 |
<!-- END GENERATED: closest-pairs -->

**Atom roles in 25V.** These are checked by a test, which confirms them
from bond lengths and template-assigned bond orders. `25V:O1` is the
**acetamide carbonyl oxygen** (double-bonded to C1). `25V:O` is the
**N-oxide oxygen** (bonded to N1). An earlier draft of this report had
these two swapped. The swap affected the interpretation only, not the
tables.

### Interpretation (geometric; no energies inferred)

- **4MEN (compound 3).** The closest pair is Asn140 ND2 to the
  triazolopyrimidine ring nitrogen N2. The distance is typical of a
  hydrogen bond, but the label is unverified: no angles or hydrogens were
  used. This places a ring nitrogen where the KAc carbonyl usually sits.
  The ligand also packs against Phe83, Trp81, Leu92 and Pro82, which
  includes the WPF shelf (Trp81/Pro82/Phe83).
- **4MEO (compound 5 / 25V).** The closest pair is the **acetamide
  carbonyl O1 to the Asp88 backbone N**. This matches the publication's
  description that the acetyl group "was not acting as the KAc mimetic but
  formed an interaction with the backbone NH donor of Asp88". The
  **N-oxide oxygen `O`** sits at hydrogen-bond-compatible distance from
  Asn140 ND2. Section 2 explains why this contact needs particular care:
  in the paper's account, the density at this position was attributed to
  a water.
- **Asn140 side-chain orientation.** At these resolutions, X-ray density
  cannot reliably tell the Asn amide OD1 and ND2 apart. Calling Asn140
  ND2 the donor rests on the depositors' assignment plus chemical
  plausibility, since both ligand atoms are acceptors. It is not measured.
- **Waters.** Conserved waters in the KAc pocket are well documented for
  bromodomains. Several waters sit within 3.5 Å of each ligand (table
  below). Water-mediated contacts are not part of the residue counts.

### Model quality and close-contact screen

<!-- BEGIN GENERATED: model-quality -->
| Entry | Ligand occupancy | Ligand B min / mean / max (Å²) | Waters ≤ 3.5 Å of ligand | Other HET groups, min distance (Å) |
|---|---|---|---|---|
| 4MEN | 1.00 | 9.6 / 14.3 / 18.3 | 6 | EDO201 10.91 |
| 4MEO | 1.00 | 13.6 / 26.1 / 37.2 | 4 | EDO201 10.63, DMS203 3.06 |
<!-- END GENERATED: model-quality -->

<!-- BEGIN GENERATED: vdw-overlaps -->
| Entry | Protein atom | Ligand atom | Distance (Å) | Bondi sum (Å) | Overlap (Å) |
|---|---|---|---|---|---|
| 4MEO | A:LEU92:CD1 | 25V:C3 | 2.945 | 3.40 | 0.455 |
| 4MEO | A:LEU92:CD1 | 25V:O1 | 2.910 | 3.22 | 0.310 |

Screen: carbon-containing heavy-atom pairs with Bondi overlap ≥ 0.25 Å. No pairs flagged in: 4MEN.
<!-- END GENERATED: vdw-overlaps -->

**4MEO Leu92 CD1 close approach.** The Leu92 CD1 methyl carbon sits
closer than the sum of Bondi van der Waals radii to both the acetamide
carbonyl O1 and the phenyl carbon C3. The C···O distance is short for a
methyl–carbonyl contact. Possible explanations, which the coordinates
alone cannot separate:

1. a genuine tight packing contact (for example a C–H···O contact);
2. coordinate error in a less well-ordered ligand: the 25V B-factors are
   higher than those of 25K. Per-atom values are in `summary.json` under
   `ligand_model_quality.b_factor_by_atom`. The highest-B atoms are the
   acetamide N, the N-oxide O, phenyl C2/C3 and O1, i.e. the region facing
   Leu92 plus the disputed N-oxide oxygen;
3. restraint or model bias during refinement, or a subtly misplaced
   ligand or side chain.

*External corroboration, not reproduced by this pipeline:* the wwPDB
validation report for 4MEO (`4meo_validation.xml.gz`, report generated
9 Mar 2026, retrieved 2026-10-14 from
`files.rcsb.org/pub/pdb/validation_reports/me/4meo/`) lists all-atom
clashes between Leu92 HD11 and 25V C3 (0.79 Å overlap) and C2 (0.48 Å).
It also lists a 0.94 Å clash between 25V H7 and DMS203, and a ligand
real-space correlation coefficient (RSCC) of 0.872. The 4MEN validation
report gives an RSCC of 0.953 for 25K. These figures support treating the
4MEO ligand geometry with more caution than 4MEN. No interaction energy
is inferred from any of them.

## 5. Limitations

1. **Heavy-atom geometry only.** There are no hydrogens, angle criteria or
   protonation model. Distances "compatible with" hydrogen bonds are not
   hydrogen-bond assignments. Coordinate uncertainty (roughly 0.1–0.2 Å at
   this resolution; not estimated per atom) is not propagated into contact
   flags, so pairs within that margin of a threshold could change class.
2. **Ligand chemical identity for 4MEO is uncertain** (Section 2). The
   CCD graph check cannot resolve it.
3. **Waters and additives are excluded from the residue counts.** Their
   proximity is reported in the model-quality table. In 4MEO, DMS203 lies
   close to the ligand.
4. **Single model and single pose.** Ligand occupancy is shown in the
   model-quality table. Protein alternate conformers occur only at
   residues at least 9 Å from the ligand (4MEN Gln78; 4MEO Gln59, Leu71,
   Thr166), and Biopython's default altloc selection is used.
5. **Descriptors** use the deposited CCD definitions with RDKit's default
   charge and tautomer handling. They are not experimentally measured
   properties.
6. **Two complexes are not statistics.** The baseline describes these two
   structures and is a reference point for Phase 2, not a binding model.

## 6. Relevance to early hit identification

These two complexes provide **observed, quantitative reference
geometries** for testing whether a computational method can reproduce
experiment. That includes recovering the observed ligand placement, the
Asn140-proximal acceptor in each complex, the WPF-shelf packing in 4MEN,
and the unusual acetamide–Asp88 contact in 4MEO. The publication itself
reports that docking reproduced the crystal poses of compound 3 closely,
and of compound 5 less exactly, partly because of bridging waters.

These features **describe the two observed complexes. They are not
universal requirements for a BRD4 hit.** Other BRD4 ligands bind with
different anchoring patterns, water networks and property profiles. 4MEO
itself shows a binding feature that departs from the usual pattern. The
descriptor values in Section 3 characterise these two ligands and should
not be read as acceptance criteria. Using them as filters in a screen
would be a separate, untested modelling choice.
