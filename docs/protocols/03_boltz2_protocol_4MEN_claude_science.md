# Deliverable 3 — Minimal Boltz-2 experimental protocol (4MEN as primary reference)

**Status: specification only. No GPU inference was run, no paid job was submitted, no compute was incurred.**

**Purpose, stated honestly.** This experiment asks whether a Boltz-2 pipeline, run from sequence, reproduces
the experimentally observed BRD4-BD1/compound-3 complex (PDB 4MEN). Because 4MEN was released 2013-09-25 and
Boltz-2's documented structural training cutoff is 2023-06-01, **this is a positive control that validates
the harness, not a test of model capability.** Arm D below is the only arm that carries generalisation
evidence. Design the write-up around that distinction from the start.

---

## 1. Protein sequence and construct boundaries

**[C] Measured from the 4MEN deposition** (UniProt O60885, BRD4_HUMAN), single chain, author residues
42–168, **127 residues, no chain breaks**, expressed in *E. coli*, monomeric biological assembly:

```
SMNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEE
```

- Author numbering equals UniProt numbering from residue 44. Residues 42–43 (`SM`) are a **cloning remnant**,
  not native BRD4.
- Index mapping for all downstream analysis: `seq[i - 42]` for author residue *i*. Verified: residue
  81 = W, 88 = D, 97 = Y, 140 = N. **Use this mapping, not a hand count** — Boltz output is 0-indexed and
  renumbered from 1, so every residue reference must be translated explicitly.
- **Construct decision (make it once, document it):** run the crystallised construct *including* `SM`, so
  that residue indices map 1:1 onto the reference without an offset table. Record the alternative
  (44–168, native-only) as a sensitivity check; do not mix the two between arms.

## 2. Ligand structure and stereochemistry

**[C] CCD 25K**, verified against the paper's compound **3**
(*N*,5-dimethyl-*N*-(4-methylbenzyl)-[1,2,4]triazolo[1,5-a]pyrimidin-7-amine):

- Formula **C15H17N5**, MW 267.33, **20 heavy atoms**
- SMILES: `Cc1ccc(cc1)CN(C)c2cc(nc3n2ncn3)C`
- InChIKey: `JAEVUMCKPBRDAJ-UHFFFAOYSA-N`
- **No stereocentres** — 25K is achiral, so no stereochemical specification is needed and no enantiomer
  control is required. State this explicitly rather than leaving it implicit.
- Tautomer/protomer: neutral form as above. The triazolopyrimidine has no ionisable centre near
  physiological pH; no protomer enumeration needed. (Contrast with compound **1**, a pyridone, where it
  would be.)

**Specify the ligand by SMILES, not by CCD code.** Passing `ccd: 25K` hands the model the deposited
chemical-component template for a ligand that exists in exactly one PDB entry — the very entry you are
trying to reproduce. SMILES input is the weaker, and therefore more honest, prior. Run the CCD variant only
as a labelled sensitivity arm.

## 3. Model version and inference settings

Two execution routes; pick one and keep it fixed across arms.

**Route A — self-hosted open-source `boltz` (recommended for this project).** MIT weights, full control,
no per-job cost, and the settings are inspectable and quotable in a write-up.

```yaml
# 4men_cofold.yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: SMNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEE
  - ligand:
      id: L
      smiles: 'Cc1ccc(cc1)CN(C)c2cc(nc3n2ncn3)C'
properties:
  - affinity:
      binder: L
```

```bash
boltz predict 4men_cofold.yaml --use_msa_server --out_dir out/ \
     --recycling_steps 3 --diffusion_samples 5 --output_format pdb
```

**Route B — hosted Boltz API via the connector.** The server exposes **`boltz-2.1`** (not plain Boltz-2 —
record the version difference in your methods) with documented defaults: `recycling_steps` 3,
`sampling_steps` 200, `step_scale` 1.638, `num_samples` 1–10, up to 4 templates, and
`binding: ligand_protein_binding` for the affinity head (≤2048 heavy atoms; the open-source path documents a
tighter 128-atom affinity cap in v2.2.x — 25K's 20 atoms clear both).

> **Blocked:** `boltz_get_account_context` returned **HTTP 401, "Bearer token authentication failed"**, so I
> could **not** obtain a live cost estimate from `boltz_estimate_structure_and_binding`. Reconnect the Boltz
> API connector (Customize → Connectors) and the estimate call is free and returns an exact figure before any
> paid run. I have not guessed a price.

## 4. Template and MSA treatment

This is the axis that determines what the experiment means, so it is the factor to vary:

| Arm | Template | MSA | What it tests |
|---|---|---|---|
| **A (primary)** | none | full (`--use_msa_server`) | Sequence-only reproduction of 4MEN |
| **B** | none | `msa: empty` | MSA dependence for a well-represented fold |
| **C** | 3MXF, protein chain only | full | The Vidler-equivalent setup: the same holo receptor Glide was given |
| **D (generalisation)** | none | full | Same pipeline on a **post-cutoff** complex — see §9 |

Notes that matter:
- **Never template on 4MEN itself.** That is the answer.
- Arm C is the scientifically interesting comparison to the paper: Vidler's docking achieved rmsd ≤ 0.30 Å
  *because* it was given 3MXF, whose protein differs from the triazolopyrimidine-bound structure by only
  0.37 Å. Arm C reproduces that information advantage; Arm A withholds it. The A-vs-C difference is the
  quantitative cost of not having a holo structure.
- `msa: empty` is an **accuracy** manipulation, not a memory saving — the MSA search is CPU-side and peak
  VRAM is unchanged. Do not use it to fit a smaller GPU.
- Record the MSA actually used. `--use_msa_server` queries `api.colabfold.com`, which is a **moving target**:
  the same command six months apart can return different alignments. Save the `.a3m` and re-use it via the
  `msa:` field for every arm, or the experiment is not reproducible.

## 5. Reproducibility controls

1. **Seeds.** 5 independent seeds (e.g. 0–4) per arm, reported as median and full range — never a single
   run. Note the API's own caveat: *results can vary across hardware and model versions even with the same
   seed*, so a seed is a weak guarantee, not a strong one.
2. **Pin and record:** `boltz` package version, model-weight checksum, CUDA/driver version, GPU model,
   `--no_kernels` on/off (the reference path is ~2× slower but numerically identical), and the saved `.a3m`
   hash.
3. **Freeze the reference.** Hash the 4MEN coordinate file you compare against; a later PDB revision would
   silently change your RMSDs. 4MEN is at revision 1.3 (2023-09-20).
4. **Pre-register the thresholds** in §7 before looking at any output.
5. **One writer per output directory**; never overwrite an arm.

## 6. Protein structural alignment

Alignment must be done on the protein and the resulting transform **applied** to the ligand — the ligand is
never aligned to itself.

1. Superpose predicted chain A onto 4MEN chain A on **Cα atoms**, residues 42–168 (all 127, no gaps).
   Report backbone Cα RMSD. Expectation for a well-represented fold like this: sub-Å; treat >1.5 Å as a
   pipeline fault, not a finding.
2. Superpose again on **binding-site residues only** — the 12 residues within 4.2 Å of 25K in 4MEN:
   **Trp81, Pro82, Phe83, Val87, Leu92, Leu94, Tyr97, Tyr139, Asn140, Asp145, Ile146, Met149** (all-atom or
   Cα; state which). Report this separately. A global fit can look excellent while the pocket is subtly
   wrong, and it is the pocket that the ligand RMSD depends on.
3. Report **both** alignments with every ligand RMSD. A ligand RMSD quoted without naming its superposition
   is uninterpretable.

## 7. Symmetry-aware ligand RMSD

25K contains topologically equivalent atoms (the *p*-tolyl ring has a two-fold; the 5-methyl and the
*N*-methyl are internally equivalent to nothing, but the aromatic ring permutation is real), so a naive
atom-order RMSD will overstate the error.

- Compute **symmetry-corrected heavy-atom RMSD** by minimising over the molecular graph automorphisms —
  RDKit `GetBestRMS` (which enumerates substructure matches) or `spyrmsd`'s symmetry-corrected RMSD.
  Hydrogens excluded.
- **Do not let the RMSD routine re-align the ligand.** Apply the protein-derived transform, then compute RMSD
  in place. RDKit's `GetBestRMS` *does* align by default — use the symmetry matching but supply the
  pre-transformed coordinates and suppress further alignment, or compute the automorphism-minimised RMSD
  directly.
- **Pre-registered thresholds** (PoseBusters convention): **< 2.0 Å = pose reproduced**; 2–5 Å = right site,
  wrong pose; > 5 Å = wrong site. Report the full distribution across seeds and the per-arm best and median.
- Report a **PoseBusters-style validity check** alongside RMSD (bond lengths/angles, planarity, internal and
  protein–ligand clashes). A 1.5 Å RMSD with distorted geometry is not a success.

## 8. Binding-site contact recovery

**[C] Reference contact set, computed by me from 4MEN coordinates** (heavy-atom distance < 4.2 Å):

| Residue | Min. distance to 25K (Å) |
|---|---|
| Trp81 | 3.42 | 
| Pro82 | 3.48 |
| Phe83 | 3.35 |
| Val87 | 3.78 |
| Leu92 | 3.47 |
| Leu94 | 3.51 |
| Tyr97 | 4.11 |
| Tyr139 | 3.76 |
| Asn140 | **2.84** |
| Asp145 | 3.75 |
| Ile146 | 3.60 |
| Met149 | 3.64 |

Metrics:
- **Contact recovery** = |predicted ∩ reference| / |reference|, same 4.2 Å cutoff. Report precision too, so
  that a model which contacts everything cannot score well.
- **Specific polar anchor:** the KAc-mimetic hydrogen bond **25K N2···Asn140 ND2 = 2.84 Å**. Report whether
  it is present and its predicted distance. This single contact is the chemical content of the paper's
  claim that the triazolopyrimidine is a KAc mimetic; a prediction that recovers the pocket but misses this
  has missed the point.
- **Acknowledge the water explicitly.** In 4MEN, 25K N1 hydrogen-bonds the conserved structural water
  (HOH316, 2.82 Å, B = 6.2 Å² — the best-ordered water in the site). Boltz-2 does not output structural
  waters, so this contact **cannot** be recovered. Exclude it from the denominator and say so, rather than
  scoring the model against something it cannot produce.

## 9. Training-data overlap — risks and the one arm that controls for them

**[C] Verified facts, not assumptions:**
- 4MEN and 4MEO: deposited 2013-08-27, **released 2013-09-25**.
- Boltz-2 structural training cutoff: **2023-06-01** (Boltz-2 preprint, Passaro & Wohlwend et al. 2025);
  Boltz-1's was 2021-09-30. Both entries sit well inside either.
- **Affinity-label leakage is also real and verifiable:** this paper's assay data is in ChEMBL as document
  **CHEMBL2434916** — 20 activities including 10 IC50 values against BRD4-BD1 (target CHEMBL1163125). Boltz-2's
  affinity head trains on ChEMBL-derived data. Any affinity prediction for this compound series must be
  reported as **plausibly in-distribution**, and the `affinity_pred_value` for 25K should not be presented as
  a prospective prediction.
- BRD4-BD1 is one of the most heavily represented pockets in the PDB, so even a held-out ligand sits on a
  thoroughly memorised receptor. Pocket reproduction is near-guaranteed and carries little information.

**Arm D — the generalisation control.** Run the identical pipeline on a BRD4-BD1 complex released *after*
2023-06-01. **[C]** I queried RCSB: 62 X-ray BRD4 (O60885) entries better than 2.0 Å were released after the
cutoff. Concrete candidates, all single-chain BD1 with one ligand of interest:

| PDB | Released | Resolution | Ligand | Ligand RSCC |
|---|---|---|---|---|
| **32SC** | 2026-08-05 | 1.25 Å | A1IGX | 0.933 |
| 24KM | 2026-07-29 | 1.86 Å | NIM | 0.939 |
| 9VCZ | 2026-06-10 | 1.50 Å | A1ERJ | — |

**32SC is the recommended held-out case** (highest resolution, single ligand, well-fitted). Before using it,
check the ligand's maximum Tanimoto similarity to known BET chemotypes — a post-cutoff *entry* with a
pre-cutoff *chemotype* still leaks. Report that similarity alongside the result.

**Interpretation rule, fixed in advance:** success on 4MEN (Arms A–C) licenses the statement *"the pipeline
is correct and the model reproduces a memorisable complex"*. Only Arm D licenses any statement about
predictive capability, and only with the chemotype-similarity caveat attached.

## 10. Model confidence and limitations

Read `confidence_*_model_0.json`. For this **single-protein-chain** input, `iptm` is still defined for the
protein–ligand interface, but note the skill's warning that a strictly single-chain input yields no ipTM —
verify which you have rather than assuming, and read `ptm` and `complex_plddt` regardless.

- Community pass lines: `iptm` > 0.5 for an interface, `complex_plddt` > 0.7 for the fold; `confidence_score`
  is the aggregate used for ranking.
- **Calibrate, don't trust.** Record the rank correlation between `confidence_score` and the actual
  symmetry-corrected RMSD across all seeds and arms. On a memorised target confidence will be uniformly high
  and therefore uninformative; demonstrating that is itself a worthwhile result.
- Known limitations to state in the write-up: no structural waters (§8); a single static output rather than
  an ensemble, against a crystal structure that is itself one ordered snapshot at 100 K in P2₁2₁2₁ with a
  DMSO and an ethylene glycol in the lattice; no crystal-packing term, so lattice-influenced side-chain
  rotamers cannot be reproduced; and the affinity head returns log₁₀(IC50/µM), which for compound **3** has
  no measured ground truth to compare against — the paper could not determine its IC50.

## 11. Compute requirements and cost

**[C] Documented facts:**
- System size: 127 residues + 20 heavy atoms ≈ **~147 tokens**, far below the 512-token/4608-atom training
  crop. This is at the small end of what Boltz-2 is built for.
- The open-source guidance addresses OOM by lowering `--diffusion_samples` / `--max_parallel_samples` or
  moving to an 80 GB tier — i.e. 80 GB is the remedy for large systems, not the requirement for small ones.
- `--no_kernels` (missing `cuequivariance_ops_torch`) is ~**2× slower**, numerically identical.
- `--use_msa_server` adds a documented **30–90 s pause per chain**, CPU-side.

**[I] Estimates — to be measured, not quoted as results.** A ~147-token complex should fit comfortably in a
**single 24 GB GPU** (A10G/L4/RTX 4090 class) at `diffusion_samples 5`; an 80 GB card is not needed. Expect
**order-of-minutes per sample**, so a 5-seed × 5-sample arm is plausibly well under an hour of GPU time, and
the whole four-arm design a few GPU-hours. **Treat every number in this paragraph as a hypothesis.**

**Required first step before any campaign:** a single timing run — one seed, one sample, Arm A — logging
wall-clock, peak VRAM (`torch.cuda.max_memory_allocated`) and MSA time separately. Scale from the measured
value; do not scale from my estimate.

**Cost:**
- *Self-hosted:* cost = (measured GPU-hours) × (your hourly rate). The only unknown is the first factor, and
  the timing run supplies it.
- *Hosted API:* call `boltz_estimate_structure_and_binding` — free, no job started, exact figure. Currently
  blocked by the 401 above. **I have not estimated a price for the hosted route.**

**Approval gate:** nothing in §1–§11 has been executed. No GPU inference, no paid job, no repository change.
Running Arms A–D requires your explicit go-ahead, and I would recommend the single timing run first.

## 12. Minimal run order

1. Timing run (1 seed, 1 sample, Arm A) → measured VRAM and wall-clock. **Gate: your approval.**
2. Arm A, 5 seeds → primary result.
3. Arm C (3MXF template) → the Vidler-equivalent comparison.
4. Arm B (`msa: empty`) → MSA dependence.
5. **Arm D (32SC)** → the only arm that speaks to generalisation.
6. Analysis: alignment (§6) → symmetry-corrected RMSD (§7) → contact recovery (§8) → confidence calibration
   (§10). Thresholds pre-registered at step 0.

## 13. What this protocol deliberately does not do

- It does not treat 4MEO as a Boltz-2 target. **4MEO is excluded from the primary experiment**: its ligand
  identity is unresolved (Deliverable 2), so there is no trustworthy ground truth to score against. If you
  want it as a case study, the only defensible framing is a *forward* one — co-fold **both** candidate
  species (compound **5**, `Cc1cc(-c2cccc(NC(C)=O)c2)c2ccccc2n1`, and the *N*-oxide **14a**, CCD 25V,
  `Cc1cc(c2ccccc2[n+]1[O-])c3cccc(c3)NC(=O)C`) and ask whether the model distinguishes them. Report that as
  a hypothesis-generating observation only: a co-folder's preference is **not** crystallographic evidence and
  must not be used to adjudicate the identity question.
- It does not benchmark Boltz-2 against Vidler's docking rmsd values. Those were obtained by docking into a
  given holo receptor; a sequence-based co-fold is a different task, and a head-to-head number would be
  misleading in both directions.
