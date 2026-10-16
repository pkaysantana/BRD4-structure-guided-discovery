# Deliverable 2 — Assessment of the 4MEO chemical-identity ambiguity

**Question.** PDB 4MEO is cited by Vidler et al. (2013) as the BRD4(1) complex of **compound 5**,
*N*-(3-(2-methylquinolin-4-yl)phenyl)acetamide. The deposited ligand, CCD **25V**, is a **quinoline
*N*-oxide**. The paper separately discusses density near the quinoline nitrogen interpreted as a bridging
water. Which chemical species was actually in the crystal?

**Answer in one line:** the deposition unambiguously *models* the *N*-oxide, that model is
crystallographically self-consistent and survives independent re-refinement, but the paper's own mass
spectrometry argues the crystallised material was not the *N*-oxide — and the publicly available,
model-phased maps **cannot discriminate** between the two interpretations. The chemical identity of the
bound species remains genuinely unresolved.

Labels: **[M]** reported by the source, **[C]** computed by me this session from public primary data,
**[I]** my inference.

---

## 1. What the original paper actually reports

**[M]** Direct quotations of substance, paraphrased closely:

- Compound **5** "demonstrated inhibition of BRD4(1), but its IC50 was lower than its solubility and could
  not be determined." So **5 has no measured potency.**
- The quinoline nitrogen "does not appear to be forming a direct interaction with the protein."
- "From the crystal structure, electron density could be observed in close proximity to the quinoline
  nitrogen (Figure 11), suggesting an interaction through a bridging water molecule."
- "the distance between the oxygen of the water molecule and the ligand appeared to be very short (2.14 Å),
  and we speculated that the active ingredient of this compound was in fact the quinoline *N*-oxide;
  **however, mass spectrometry did not detect the *N*-oxide**."
- They nevertheless synthesised the *N*-oxide **14a** (mCPBA, DCM, reflux, 3 h; 14% yield, 7.2 mg) alongside
  a resynthesis of **5**. Resynthesised **5** behaved like the original sample; **pure 14a gave
  IC50 = 42 µM** — the paper's stated conclusion being "a single-atom change… resulting in the displacement
  of a single water molecule and achieving an increase in activity."
- The deposition statement names **"4MEO (BRD4(1)/compound 5)"**.

**[I]** Three things follow immediately. First, the paper's *text* commits to compound **5** plus a bridging
water, and explicitly records evidence *against* the *N*-oxide (negative MS). Second, the 2.14 Å figure is
itself anomalous: a quinoline-N···O(water) hydrogen bond should be ≈2.7–3.0 Å, and 2.14 Å is far too short —
the authors flag this themselves as the reason for suspecting the *N*-oxide. Third, the 42 µM potency of
**14a** is a *measured* value for an authentic sample of the *N*-oxide and tells us the *N*-oxide binds
better than anything they could quantify for **5**. That potency ordering is the hinge of the whole problem,
because it means even a minor *N*-oxide impurity could preferentially occupy the site.

---

## 2. What the RCSB deposition specifies

**[C]** Retrieved from the RCSB chemical component dictionary and the deposited coordinates:

| Property | CCD 25V (deposited in 4MEO) | Compound 5 in the paper | *N*-oxide 14a in the paper |
|---|---|---|---|
| Name | *N*-[3-(2-methyl-1-oxidoquinolin-4-yl)phenyl]acetamide | *N*-(3-(2-methylquinolin-4-yl)phenyl)acetamide | 4-(3-acetamidophenyl)-2-methylquinoline 1-oxide |
| Formula | **C18 H16 N2 O2** | C18 H16 N2 O | C18 H16 N2 O2 |
| HRMS [M+H]⁺ | — | calcd 277.1335, found 277.1333 | calcd 293.1285, found 293.1287 |
| Canonical SMILES | `Cc1cc(c2ccccc2[n+]1[O-])c3cccc(c3)NC(=O)C` | — | — |
| InChIKey | YSGOTUAMXVUKOD-UHFFFAOYSA-N | — | — |

**[C]** The *N*-oxide is explicit, not an artefact of nomenclature: the CCD bond block for 25V contains
`25V O N1 SING N N 12`, i.e. a formal single bond from oxygen `O` to the quinoline nitrogen `N1`, with
`N1–C7` double and aromatic and `N1–C9` single and aromatic. The SMILES carries `[n+]…[O-]`.

**[C]** Therefore **25V is 14a, not 5.** The deposited ligand is chemically identical to the synthesised
*N*-oxide, and chemically different from the compound the paper's text says the crystal contains.
25V appears in exactly one PDB entry (4MEO); 25K appears in exactly one (4MEN).

**[I]** This is a genuine, documented inconsistency between a paper and its own deposition, not a
misreading. The paper says "4MEO = compound 5"; the PDB says 4MEO contains the *N*-oxide. Both were
deposited/submitted in the same window (deposition 2013-08-27, release 2013-09-25, paper accepted
2013-10-28), so the most likely explanation is that the depositors resolved the ambiguity in favour of the
*N*-oxide at the coordinate level while the manuscript text retained the bridging-water reading. I can
evidence the discrepancy; I cannot evidence the intent.

---

## 3. What the deposited structure and crystallographic validation establish

**[C]** Geometry of the contested atom, measured from the deposited coordinates:

| Quantity | Deposited 4MEO | PDB-REDO re-refinement |
|---|---|---|
| N1–O distance | **1.358 Å** | 1.382 Å |
| B(N1) / B(O) | 13.6 / **34.5** Å² | 21.7 / 28.7 Å² |
| Occupancy of O | 1.00 | 1.00 |
| Alternate conformations of ligand | none | none |
| R-free | 0.2126 | 0.2072 |

A 1.358 Å N–O separation is a covalent bond. The deposited model is unambiguously an *N*-oxide and contains
**no water** at the contested position; the nearest water to the oxide oxygen is HOH327 at 3.03 Å.

**[C]** In the deposited model the oxide oxygen takes over the KAc-mimetic role: **O···Asn140 ND2 = 2.88 Å**
and **O···HOH327 = 3.03 Å**, where HOH327 is itself 2.46 Å from Tyr97 OH — i.e. the conserved structural
water. The remaining polar contacts are the acetamide carbonyl **O1···Asp88 N = 2.67 Å** (the novel backbone
contact the paper highlights) and **O1···HOH324 = 2.91 Å**, plus amide **N···Gln85 NE2 = 3.55 Å** and
**N···HOH324 = 2.79 Å**. Residues within 4.2 Å: Trp81, Pro82, Phe83, Gln85, Pro86, Val87, Asp88, Leu92,
Leu94, Tyr97, Tyr139, Asn140, Ile146 (plus a DMSO molecule at 3.06 Å).

**[I]** Note the structural consequence, which is more interesting than the identity question alone: in the
deposited model compound **5**'s quinoline *is* engaged with Asn140 — through its oxide oxygen — so the
paper's statement that the quinoline nitrogen makes no direct protein interaction describes a *different
model* from the one deposited. The two readings are not cosmetic variants; they assign different
pharmacophores to the same ligand.

**[C]** Validation metrics single 25V out as the weakest-fitting ligand of the four Vidler structures:

| | 4MEP/24Y | 4MEN/25K | 4MEQ/25O | **4MEO/25V** |
|---|---|---|---|---|
| RSCC | 0.949 | 0.953 | 0.965 | **0.872** |
| RSR | 0.077 | 0.055 | 0.049 | **0.097** |
| Fit percentile | 0.741 | 0.838 | 0.890 | **0.506** |
| Intermolecular clashes | 0 | 2 | — | **4** |
| Mean ligand B | — | 14.3 Å² | — | **26.1 Å²** |

4MEO's entry-level data-quality percentile (56.8) is also the lower of the pair against 4MEN (74.3), and its
RSRZ outlier fraction is twice as high (3.15% vs 1.57%) — despite 4MEO having the *better* resolution
(1.72 vs 1.81 Å) and better geometry percentile (82.0 vs 67.2).

**[I]** The ligand model in 4MEO is measurably poorer than its siblings. That is consistent with — but does
not prove — a species-assignment problem; it is equally consistent with a more mobile, lower-occupancy
ligand of weaker affinity, which is exactly what a compound with unmeasurably weak IC50 would give.

---

## 4. Can public electron-density evidence discriminate? — I tested this directly

Public data **is** available and I used it rather than speculating about it:
structure factors `files.rcsb.org/download/4MEO-sf.cif` (209 kB) and PDBe EDS maps
`www.ebi.ac.uk/pdbe/entry-files/4meo.ccp4` (2Fo−Fc) and `4meo_diff.ccp4` (Fo−Fc), 997 kB each,
grid 48 × 54 × 96 over the P2₁2₁2₁ cell (≈0.78–0.82 Å spacing), σ = 0.3162 and 0.1110 respectively.
I parsed the CCP4 maps with numpy and sampled by trilinear interpolation.

**[C] Density at the contested atoms (σ units):**

| Atom | B (Å²) | 2Fo−Fc at centre | 2Fo−Fc local max | Fo−Fc at centre |
|---|---|---|---|---|
| 25V **N1** (quinoline N) | 13.6 | **3.44** | 3.70 | −0.31 |
| 25V **O** (oxide O) | 34.5 | **1.55** | 2.47 | −0.13 |
| HOH327 (conserved water) | 14.2 | 3.42 | 3.36 | +0.49 |
| Asn140 ND2 | 19.2 | 2.33 | 2.70 | +0.25 |
| *4MEN control:* HOH316 | 6.2 | 3.31 | 4.04 | +0.12 |

**[C] Density profile along the N1 → O vector** (the direct test of 1.358 Å vs 2.14 Å):

| distance from N1 (Å) | 1.0 | 1.2 | **1.36** | 1.6 | 1.8 | 2.0 | **2.14** | 2.4 | 2.6 | 3.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| 2Fo−Fc (σ) | 1.91 | 1.65 | **1.56** | 1.55 | 1.60 | 1.70 | **1.73** | 1.39 | 1.13 | 0.71 |
| Fo−Fc (σ) | −0.63 | −0.42 | **−0.13** | +0.57 | +1.22 | +1.92 | **+2.37** | +2.71 | +2.93 | +3.39 |

Centroid of positive 2Fo−Fc density over the 0.9–2.8 Å window: **1.73 Å from N1** — between the deposited
oxide position (1.358 Å) and the paper's water position (2.14 Å).

**[C] What this does and does not show.**
- The 2Fo−Fc density along this vector is a **broad, low shoulder (1.5–1.7σ), not a resolved peak at either
  candidate distance.** For comparison, a genuinely well-ordered fully-occupied water in the same maps
  (HOH327, 3.42σ; 4MEN's HOH316, 3.31σ) is twice as strong. Whatever sits there is partially occupied,
  mobile, or both.
- Positive Fo−Fc rises steadily beyond the oxide oxygen, reaching +3.4σ at 3.0 Å from N1. **I checked what
  lies there before interpreting it**: the points at 2.6–3.0 Å are only 2.1–2.5 Å from the already-modelled
  Asn140 ND2 and 2.5–2.7 Å from HOH327. This feature is therefore *not* a cleanly separable missing-water
  site; at ~0.8 Å grid spacing it cannot be distinguished from slight mis-positioning or under-modelling of
  that Asn140/HOH327 pair.
- Fo−Fc at the oxide oxygen itself is **−0.13σ — not significantly negative.** The data therefore do **not**
  reject a fully occupied oxide oxygen.

**[C] Decisive caveat — model phase bias.** PDBe EDS maps are computed from the deposited structure factors
*with the deposited model's phases*. At the one atom whose existence is in question, the map is biased
toward showing density. The 1.55σ at the oxide position is therefore an **upper bound** on the unbiased
evidence, and the apparent centroid shift to 1.73 Å is not a reliable discriminator. **An unbiased test
requires an omit or polder map**, which I could not compute (see "Blocked", below).

**[C] Independent re-refinement weakens my first inference, and I am recording that.** My initial reading was
that B(O) = 34.5 Å² against B(N1) = 13.6 Å² — a 2.5× disparity between covalently bonded neighbours, with
oxygen the heavier atom — pointed to partial occupancy of the oxide oxygen, i.e. a mixture of **5** and
**14a**. PDB-REDO's independent re-refinement of the same data retains the *N*-oxide at full occupancy,
places N1–O at 1.382 Å, and **largely removes the disparity** (N1 13.6 → 21.7, O 34.5 → 28.7 Å²) while
improving R-free from 0.2126 to 0.2072. So a substantial part of the original B anomaly was a feature of the
depositors' REFMAC/TLS run, not of the data. The partial-occupancy hypothesis is **not supported** by that
check and I have downgraded it from "indicated" to "open".

---

## 5. Chemical reasoning that the crystallography cannot settle

**[I]** The most parsimonious reconciliation of all the evidence — offered as a hypothesis, flagged as such:

*The commercial sample of 5 (ChemBridge 50866471, 5 mg) contained a minor quinoline N-oxide impurity which,
being the more potent species (IC50 42 µM for authentic 14a versus unmeasurable for 5), was
preferentially selected into the crystal from a 5 mM soak at pH 8.0, 4 °C.*

This hypothesis is consistent with every observation: the too-short 2.14 Å contact (a covalent N–O seen as a
"water"), the negative bulk MS for the *N*-oxide (a few per cent impurity is easily missed in a bulk sample
while dominating occupancy in a crystal, because crystallisation selects on affinity, not abundance), the
depositors' decision to model an oxide, the weaker RSCC and higher ligand B, and the broad rather than
peaked density. It predicts a *mixed* site, which is also why neither pure model fits perfectly.

**[I]** Competing explanations that remain live:
1. **Pure *N*-oxide** in the crystal, exactly as deposited — supported by the covalent-distance geometry,
   PDB-REDO's agreement, and the non-negative Fo−Fc at the oxide O; contradicted by the negative MS.
2. **Compound 5 plus a genuinely short-contact water** — the paper's reading; requires explaining a 2.14 Å
   N···O distance, which is hard on physical grounds alone.
3. **Mixed occupancy of 5 and 14a** — explains the weak/broad density best; not currently supported by
   occupancy refinement, and 1.72 Å data with ~0.8 Å map sampling may simply be unable to resolve it.
4. **Radiolytic or soak-time oxidation** of the quinoline during crystallisation or data collection —
   raised for completeness; I have no evidence for or against it.

---

## 6. What remains uncertain, and what would settle it

**Unresolved (honest status):**
- The chemical identity of the species bound in the 4MEO crystal. **Not established** by the deposition, the
  paper, or any public evidence I could retrieve.
- Whether 4MEO's site is single-species or mixed.
- Whether the paper's Figure 11 shows a map consistent with the finally deposited coordinates, since the
  figure was not retrievable.

**Settled:**
- CCD 25V *is* the *N*-oxide (14a), chemically distinct from compound **5**. **Established.**
- The paper's text and its deposition disagree about which compound 4MEO contains. **Established.**
- The deposited and PDB-REDO models both place a covalent N–O bond with no bridging water. **Established.**
- 25V is the poorest-fitting ligand of the four Vidler depositions. **Established.**

**Decisive experiments, in priority order (none require new crystallography):**
1. **Polder / omit map at the oxide oxygen** from `4MEO-sf.cif` — remove `O` of 25V, refine, and compute an
   unbiased difference map. This is the single measurement that removes phase bias. Phenix `polder_map`, or
   CCP4 with Servalcat/REFMAC.
2. **Competing-model refinement with occupancy refinement**, scored on R-free, local RSCC/RSR and residual
   Fo−Fc: (a) *N*-oxide, occ 1.0; (b) neutral quinoline + water at refined occupancy; (c) *N*-oxide at
   refined occupancy; (d) two-component mixture. Report ΔR-free with its σ (PDB-REDO gives σ(R-free) = 0.0077
   for this entry, so differences below ~0.008 are not interpretable).
3. **Re-analyse the actual sample.** LC-MS/NMR of ChemBridge 50866471, ideally the same lot, with an
   authentic 14a standard and a detection limit stated. This is the only experiment that tests the impurity
   hypothesis directly, and it is cheap.
4. Contact the depositors (Knapp/Filippakopoulos, SGC Oxford) — the ligand assignment rationale may simply be
   recorded in their notes.

---

## Blocked tonight (environment, not evidence)

Items 1 and 2 above require `gemmi`/`CCP4`/`Phenix`/`servalcat` for map calculation and refinement. Package
installation was unavailable this session by your instruction, so **no omit map, no occupancy refinement and
no re-refinement were performed by me.** What I did instead — parsing the PDBe maps with numpy and sampling
them — is reported above with its phase-bias limitation stated. The structure factors, maps and PDB-REDO
model are all downloaded and cached in `data/` for whenever the environment is available.
