# 4MEN primary full-MSA experiment: post-run report

Experiment: `4men_msa_s20261019_v1`
Branch: `phase3-4men-msa-20261017`
Report date: 17 October 2026

The experiment was prepared under the Phase 2 harness and executed as the
Phase 3 milestone, so its files are under `experiments/phase2/` and
`docs/phase2/`. This report describes those files. It does not change them.

Unless stated otherwise, every value below is read from a committed file
in the experiment directory
[`experiments/phase2/4men_msa_s20261019_v1/`](../../experiments/phase2/4men_msa_s20261019_v1/),
which is abbreviated to `EXP/`. Values marked **external** come from audit
reports supplied to the project and are not reproduced in this repository.

## A. Executive summary

Boltz-2 (v2.2.1) was given the BRD4 BD1 sequence, a pinned full MSA and
the SMILES of ligand 25K. No structure, template or constraint was
supplied. It reconstructed the deposited 4MEN complex closely. After a
Cα fit of the protein onto the crystal structure, the ligand heavy-atom
RMSD was **0.140 Å** and the protein Cα RMSD was **0.191 Å**. At the
primary 4.0 Å cutoff, the predicted ligand-contacting residues were the
same 11 residues as in the crystal structure. The ligand passed the
heavy-atom geometry screen. Every preregistered gate passed, so the
preregistered **harness success is True**.

The run was a single retrospective reconstruction of one known complex,
with one seed and one sample. 4MEN was public about ten years before the
reported Boltz-2 structural training cutoff. The result shows that the
pipeline runs end to end with complete provenance and that its scoring
behaves sensibly on a known case. It does **not** show that Boltz-2 would
predict new BRD4 ligands prospectively, that it generalises beyond its
training data, or anything about affinity. The 3.5 Å contact sensitivity
analysis (precision 0.667, recall 0.40) shows that agreement at the
residue level does not extend to sub-ångström agreement on every contact.

## B. Scientific background

**Target.** BRD4 is a BET-family acetyl-lysine (KAc) reader. Its first
bromodomain (BD1, UniProt O60885) has a largely hydrophobic KAc pocket.
The pocket lies between the ZA and BC loops, with polar recognition by the
conserved Asn140 and a structurally conserved water. The WPF shelf
(Trp81–Pro82–Phe83) and the ZA channel line the pocket entrance. See the
[evidence brief](../research/_claude_science_01_evidence_brief_vidler2013.md)
and the [reviewed integration memo](../research/REVIEWED_INTEGRATION.md),
whose corrections take precedence.

**Study.** Vidler et al. (2013, *J Med Chem* 56:8073, doi:10.1021/jm4011302)
ran a structure-based virtual screen for KAc mimetics. They combined
substructure, pharmacophore, shape and 2D-similarity searches templated on
(+)-JQ1 from 3MXF, followed by Glide docking into 3MXF and expert triage.
Hits were confirmed crystallographically. Their published results are
attributed to that paper and were not reproduced here.

**Complex.** PDB 4MEN is BRD4 BD1 bound to ligand **25K**, compound 3 of
Vidler et al.
(*N*,5-dimethyl-*N*-(4-methylbenzyl)-[1,2,4]triazolo[1,5-*a*]pyrimidin-7-amine;
C15H17N5; 20 heavy atoms; achiral). The pinned file
[`data/raw/4MEN.pdb`](../../data/raw/4MEN.pdb) records X-ray diffraction
at 1.81 Å, deposition on 27 August 2013 and initial release on
25 September 2013. The construct is the `SM` cloning remnant followed by
O60885 residues 44–168 (127 residues, all Cα present).

**Why this complex was chosen as a control.** 25K is small, neutral and
achiral. External validation statistics quoted in the evidence brief
(not reproduced locally) indicate that it fits its density well. The construct is complete and
already characterised in Phase 1. These properties make 4MEN a clean
control for the harness, with known ground truth, before any held-out
case is attempted. 4MEO is excluded because its ligand identity is
unresolved (see the Phase 1 report). Choosing a well-known, pre-cutoff
complex also means that success is expected whether the model learned
transferable features or recalls the training data. The experiment
therefore tests the pipeline, not the model's capability.

## C. Experimental design (frozen before inference)

The configuration is
[`configs/phase2/4men_primary.json`](../../configs/phase2/4men_primary.json).
The protocol is [`PROTOCOL.md`](PROTOCOL.md). The frozen bundle
(`EXP/config.json`, `EXP/preregistration.json`, `EXP/command.json`,
`EXP/manifest.json`, inputs and maps) was committed in `c9b09dc` at
00:59 UTC on 17 October 2026. That commit precedes the inference start
recorded at 01:34:33 UTC.

| Item | Frozen value |
|---|---|
| Protein (chain A) | 127 residues: `SM` remnant + O60885 44–168 |
| Ligand (chain B) | SMILES `Cc1ccc(CN(C)c2cc(C)nc3ncnn23)cc1`; InChIKey `JAEVUMCKPBRDAJ-UHFFFAOYSA-N`; neutral |
| MSA | Local A3M `data/phase2/msa/4men_full.a3m`, SHA-256 `8d0374e9…c9d`; 10,877 records including the query; ColabFold MMseqs2 env mode via `api.colabfold.com`, Boltz 2.2.1 client; retrieved 2026-10-17T00:45:46Z in 2.5 s; **database release unknown** |
| MSA row limit | `--max_msa_seqs 8192` (upstream default; the full file is archived even though only the first 8,192 rows can be used) |
| Model | `boltz2`, Boltz 2.2.1; structure checkpoint `boltz2_conf.ckpt` from `boltz-community/boltz-2` at revision `6fdef46d…5607` |
| Seed | 20261019 |
| Sampling | 1 diffusion sample, 1 parallel sample, 3 recycling steps, 200 sampling steps, step scale 1.5 |
| Off | Templates, constraints, steering potentials, affinity head, optimised kernels (`--no_kernels`) |
| Pose selection | The sole `model_0`; no reranking or selection against the reference |
| Protein alignment | Unweighted proper Kabsch fit on all 125 native Cα (author residues 44–168); the transform is applied to the ligand; no separate ligand fit |
| Ligand mapping | Minimum RMSD over all exact heavy-atom graph isomorphisms (2 evaluated); atom names reproduced from Boltz's SMILES parser and recorded in `EXP/ligand_atom_names.json` |
| Contacts | Unique protein residues with any heavy atom ≤ cutoff from any ligand heavy atom; no waters, additives or inferred H-bonds; primary 4.0 Å, sensitivity 3.5 and 5.0 Å |

Preregistered acceptance gates (all must pass for harness success):

| Gate | Threshold |
|---|---|
| Ligand RMSD after the protein fit | ≤ 2.0 Å |
| Protein Cα RMSD | ≤ 2.0 Å |
| Contact precision / recall at 4.0 Å | ≥ 0.70 / ≥ 0.70 |
| Bond-length deviation from CCD ideal | ≤ 0.25 Å |
| Bond-angle deviation | ≤ 30° |
| Aromatic ring planarity | ≤ 0.15 Å |
| Intraligand nonbonded distance / Bondi sum | ≥ 0.75 |
| Protein–ligand Bondi overlap | < 0.80 Å for every pair |

Confidence outputs were preregistered as diagnostics only, with no
threshold.

## D. Execution and infrastructure

**Measured** (from `EXP/run_observation.json`, `EXP/runtime.json` and
`EXP/logs/boltz_predict.log`):

| Item | Value |
|---|---|
| Runtime provenance captured | 2026-10-17T01:26:03Z (`runtime.json`; receipt `runtime.sha256` = `2172a3dd…c9c`, matching the file) |
| Inference start / finish | 2026-10-17T01:34:33Z / 01:36:02Z |
| Wall time | 89.03 s (whole `boltz predict` process, including start-up, preprocessing and model loading) |
| Exit code | 0 |
| GPU | NVIDIA GeForce RTX 5070 Ti Laptop GPU, 12,227 MiB, driver 582.05, UUID `GPU-7ca5cd76-…` |
| Host | WSL2 Linux 6.18.33.2, glibc 2.43, Python 3.12.14 |
| Torch / CUDA / cuDNN | torch 2.7.1+cu128, CUDA build 12.8, cuDNN 9.7.1 (90701) |
| Key packages (inference env) | boltz 2.2.1, rdkit 2026.3.6, numpy 1.26.4, biopython 1.84, pytorch-lightning 2.5.0 |
| Checkpoint SHA-256 | `boltz2_conf.ckpt` `090e82ac…28e1`; `boltz2_aff.ckpt` `dcc5cd37…ca9e` (present, unused); `mols.tar` `39e076d9…1d7` |
| CCD cache | 45,227 extracted files verified against `mols.tar` (tree digest `d11de775…0ca2`) |
| Precision | bfloat16 automatic mixed precision (log) |
| Completion evidence | Log: `GPU available: True (cuda), used: True`, one prediction batch at 100%, `Number of failed examples: 0`; prediction, confidence and PAE/PDE/pLDDT files present |

**Sampled or approximate, not exact:**

- **GPU memory.** The maximum of 26 `nvidia-smi` samples taken every 3 s
  was **6,900 MiB**. This is device-wide usage, and the true instantaneous
  peak may have been higher. It is not a PyTorch allocator peak.
- **Model compute time.** The log's progress bar shows the single
  prediction step at about 20 s. This figure is coarse and comes from the
  log only. The rest of the 89 s wall time is process start-up, input
  processing, checkpoint loading and output writing. This was not timed
  separately.

**Recorded departures and gaps relative to the protocol:**

- The 12 GB laptop GPU is below the configuration's
  `recommended_vram_GB: 24`. That value was a recommendation, not a gate,
  and the run completed.
- The protocol asked for inference-time network access to be blocked.
  The recorded policy is `HTTP proxy refusal and HF offline flags; not full
  isolation`. The cache was verified before the run, but the run was not
  network-isolated.
- The protocol's analysis command names the Windows `.venv`. Analysis
  actually ran on the WSL2 Linux host (`analysis_platform` in
  `metrics.json`) with the same pinned analysis packages (biopython 1.88,
  numpy 2.5.3, rdkit 2026.3.6) and the source hashes in the manifest.
- `run_observation.json` has no cost field. The run used local hardware
  with no billing, so there is no billed duration or cost.
- Determinism settings (`CUBLAS_WORKSPACE_CONFIG`, `PYTHONHASHSEED` and
  others) were unset. Together with bf16 mixed precision, this means the
  seed does not guarantee bit-identical coordinates on rerun, as the
  protocol already states.
- Lightning warned that the checkpoint was written with v2.5.0.post0 and
  loaded with v2.5.0. No effect on outputs is known. It is recorded here
  for completeness.
- `runtime.json` states that user approval "must be recorded separately".
  No separate approval record is archived. See the execution-provenance
  statement below.

**Execution provenance (retrospective statement, written 17 October
2026).**

- The project owner gave an interactive instruction to run this
  experiment, and the run followed that instruction.
- No separate, contemporaneous, signed approval record exists in this
  repository. This report is not a substitute for one, and no such record
  has been created or backdated.
- The frozen fields `"execution_authorized": false` (in `manifest.json`
  and `validation.json`) and `"executed": false` (in `command.json`)
  describe the bundle as it was prepared. They were deliberately left
  unchanged after the run.
- The evidence that inference actually took place is the material listed
  in this section: `run_observation.json`, the log, and the output hashes
  bound in `metrics.json`. It is not an approval artefact.

## E. Structural results

All values are from `EXP/analysis/metrics.json`. The prediction scored is
`model_0.cif` (SHA-256 `2d06e475…b28b`).

### Primary metrics and gates

| Metric | Value | Gate | Pass |
|---|---|---|---|
| Ligand heavy-atom RMSD after the protein Cα fit | 0.140485 Å | ≤ 2.0 Å | Yes |
| Protein Cα RMSD (125 Cα) | 0.191367 Å | ≤ 2.0 Å | Yes |
| Contact precision at 4.0 Å | 1.0 (11/11) | ≥ 0.70 | Yes |
| Contact recall at 4.0 Å | 1.0 (11/11) | ≥ 0.70 | Yes |
| Geometry screen | valid; no failed checks | all subchecks | Yes |
| **Harness success** | **True** | all gates | — |

Pocket Cα RMSD (reference residues within 8 Å of the ligand, same global
fit, diagnostic only): 0.128328 Å.

### Geometry screen detail

| Check | Prediction | Threshold | Deposited 4MEN (pre-run, `validation.json`) |
|---|---|---|---|
| Max bond-length error | 0.0541 Å | ≤ 0.25 Å | 0.1176 Å |
| Max bond-angle error | 3.50° | ≤ 30° | 7.18° |
| Max aromatic plane deviation | 0.0073 Å | ≤ 0.15 Å | 0.0158 Å |
| Min nonbonded / Bondi ratio | 0.791 | ≥ 0.75 | 0.812 |
| Max protein–ligand overlap | 0.023 Å | < 0.80 Å | 0.256 Å |
| Severe overlap pairs | 0 | 0 | 0 |

The predicted ligand is closer to the CCD ideal geometry than the
deposited model is. A plausible explanation, not tested here, is that a
generative model trained on refined and idealised structures will tend to
produce idealised geometry. It does not show that the predicted geometry
is more correct than the experimental one. The screen covers heavy atoms
only: no force-field energy, hydrogens, torsion-strain or
hydrogen-bond assessment.

### Contact sensitivity analysis

| Cutoff | Reference residues | Predicted residues | TP | FP | FN | Precision | Recall |
|---|---|---|---|---|---|---|---|
| 3.5 Å | 5 | 3 | 2 | 1 | 3 | 0.667 | 0.40 |
| 4.0 Å (primary) | 11 | 11 | 11 | 0 | 0 | 1.0 | 1.0 |
| 5.0 Å | 13 | 13 | 13 | 0 | 0 | 1.0 | 1.0 |

At 3.5 Å the reference set is {Trp81, Pro82, Phe83, Leu92, Asn140} and the
predicted set is {Phe83, Met149, Asn140}. Since every one of these
residues is in both 4.0 Å sets, the membership data alone imply the
following:

- Trp81, Pro82 and Leu92 have a minimum ligand distance ≤ 3.5 Å in the
  crystal structure but in (3.5, 4.0] Å in the prediction.
- Met149 has a minimum ligand distance in (3.5, 4.0] Å in the crystal
  structure but ≤ 3.5 Å in the prediction.

Each of these four disagreements is a distance shift across the 3.5 Å
boundary of less than 0.5 Å. The actual per-residue distances are not
stored in `metrics.json`, and this report does not compute them.

The diagnostic that restricts the prediction to atoms resolved in the
reference (`contacts_reference_observed_atoms`) gives identical sets at all
three cutoffs. The Lys91 and Gln123 side-chain atoms that are absent from
the deposited model therefore do not affect any contact result.

**Why 4.0 Å set agreement does not mean identical interactions.** The
metric compares the sets of residues that have at least one heavy atom
within the cutoff. It ignores which ligand atom makes each contact, the
exact distance below the cutoff, the angles and the number of contacting
atoms per residue. Two poses can have identical residue sets while
differing in specific atom pairs. The metric also cannot assess hydrogen
bonds: no hydrogens or angles are evaluated. In addition, the prediction
contains no waters, so it cannot reproduce water-mediated recognition.
This includes the conserved-water contact to the ligand that is described
for 4MEN in the evidence brief. Atom-level recovery of the Asn140 contact
was not measured in this experiment.

### Atom coverage

All ligand heavy atoms and all required reference protein heavy atoms are
present. The terminal `A:GLU168:OXT` is absent from the prediction, as
the protocol permits. Six predicted side-chain atoms have no reference
coordinates (Lys91 CD/CE/NZ, Gln123 CD/NE2/OE1). These were accepted
through the residue-template check and are listed in `atom_coverage`.

### Confidence outputs (diagnostic only)

From `confidence_4men_model_0.json`: confidence score 0.983, pTM 0.983,
ipTM / ligand ipTM 0.987, complex pLDDT 0.981, complex ipLDDT 0.986,
complex PDE 0.277, complex iPDE 0.303. `protein_iptm` is 0.0 because there
is only one protein chain. The pairwise chain ipTM matrix is asymmetric:
the [0][1] entry is 0.767 and the [1][0] entry is 0.987. These scores are
uncalibrated for this use, come from a single sample and a single case,
and are expected to be high on a complex that is probably in the training
distribution. They are not used as evidence of accuracy.

## F. Independent audit findings

Two independent reviews of the completed run were generated outside this
repository: one by Codex and one by a separate Claude Code session. These
reviews are external. **No transcript or output file from either review
is present in this repository.** The findings below are summarised as the
project owner supplied them. They are attributed to the reviewers, and
none of them has been reproduced by this report.

**Verification reported by the auditors (external):**

- Codex recomputed the ligand RMSD independently with a Horn
  unit-quaternion superposition instead of the harness's Kabsch routine.
  It reported agreement with the stored values.
- The auditors' numerical recomputation of the stored metrics agreed with
  `metrics.json`.

**Limitations and interpretations raised by the auditors:**

1. The 3.5 Å contact sensitivity result is scientifically important. It
   limits claims of interaction-level fidelity (section E).
2. Training-set exposure is unresolved (section G).
3. The geometry screen is limited and permissive. In an adversarial test
   that was part of the independent audit, the ligand was expanded by 12%
   about its centroid, and the distorted structure still passed the
   original gates (section H). The test came from the audit, not from this
   experiment. It is not in the repository's test suite, and its code and
   output are not archived. It should be re-implemented as a reproducible
   unit test in Phase 4.
4. Ligand-conformer preparation may be nondeterministic. According to the
   audit, Boltz builds the ligand's starting conformer from SMILES with
   RDKit, and that step may not be governed by the run seed. Neither point
   has been verified here. The
   processed ligand input actually used (`processed/mols/4men.pkl`) is
   archived, so this run's starting point is preserved even if a rerun
   would differ.
5. All original model outputs must be preserved, including the PAE, PDE
   and pLDDT arrays and the processed inputs. This has been done (section
   I).
6. One seed on one case cannot demonstrate generalisation or estimate
   variance.

**Checked in this documentation pass (local):**

- The current branch is `phase3-4men-msa-20261017` at `d2e96b8`. All 26
  experiment files are tracked, and none are untracked or ignored.
- Every value quoted in this report matches the committed JSON.
- `runtime.sha256` matches the SHA-256 of `runtime.json`, and both match
  `runtime_manifest_sha256` in `metrics.json`.
- `prediction_sha256` and `confidence_sha256` in `metrics.json` match the
  committed output files.
- `prepared_manifest_sha256` matches `manifest.json`, and the hashes of
  the frozen files listed in `manifest.json` match the committed files.
- The command recorded in `run_observation.json` and `runtime.json` is
  identical to the frozen `command.json`.

This pass did not recompute any structural metric.

## G. Training-data leakage and interpretation

Two different routes of leakage need to be kept apart.

**1. Direct leakage of experimental coordinates into inference inputs: not
found.** The input YAML (`EXP/inputs/4men.yaml`) contains only the
sequence, the SMILES string and a path to the local MSA. The processed
Boltz record (`processed/records/4men.json`) shows `templates: []`, no
pocket or contact constraints and no affinity request. An A3M alignment
carries sequences, not coordinates. The ligand was specified by SMILES,
not by the `25K` component ID. Reference coordinates were used only after
inference, for scoring.

**2. Learned exposure during model training: plausible, but neither
confirmed nor excluded.** 4MEN was released on 25 September 2013. The
Boltz-2 report describes a 1 June 2023 cutoff for experimental PDB
structural training data, plus additional distillation, MD and affinity
supervision. 4MEN is therefore eligible by date. Other BRD4 and BET
complexes, including 3MXF, 2OSS and the rest of the Vidler series, are
also eligible. No per-checkpoint training manifest has been examined, so
the membership of 4MEN in the training set of this checkpoint is
**unverified**.

Both of the following claims are unsupported by the present evidence:
"Boltz-2 memorised 4MEN" and "this result is independent of the training
data". A 0.14 Å ligand RMSD is consistent with recall of a training
example, with transferable learning of a well-sampled pocket, or with
both. A single pre-cutoff case cannot tell these explanations apart.

## H. Limitations and failed assumptions

1. **3.5 Å sensitivity.** Contact agreement is perfect at 4.0 and 5.0 Å
   but not at 3.5 Å (precision 0.667, recall 0.40). Four residues move
   across the 3.5 Å boundary. The headline 1.0/1.0 holds for the
   preregistered primary definition only.
2. **Single seed, single sample, single case.** The run gives no estimate
   of seed variance, of the sampling distribution or of failure rate. It
   cannot calibrate confidence and cannot support generalisation.
3. **The geometry screen has limited sensitivity.** In the independent
   audit's adversarial test, a uniform 12% expansion about the ligand
   centroid passed the original gates. That result is external and has not
   been reproduced here. The reasoning below is an analytical explanation,
   not a local measurement:
   - Uniform scaling leaves every bond angle unchanged and keeps planar
     rings planar, so the angle and planarity checks cannot detect it.
   - It lengthens bonds of about 1.3–1.5 Å by roughly 0.16–0.18 Å, which
     is within the 0.25 Å tolerance.
   - It moves atoms apart, so the intraligand nonbonded check cannot fail.
   - It moves each atom by 12% of its distance from the centroid, which
     adds well under the 2.0 Å RMSD gate for a 20-heavy-atom ligand.

   **The synthetic structure is a deliberately distorted test object. It
   is not the experimental prediction.** The actual prediction's largest
   bond-length error is 0.054 Å and its largest angle error is 3.5°. Its
   geometry is well formed by every check that was run. The finding is
   about the validator: it would not have caught this class of error.
   Passing the screen is therefore weak evidence of physical validity.
   A reproducible version of this perturbation is planned as a Phase 4
   unit test. Any stricter gate developed from it will apply only to
   future experiments, not retroactively to this one.
4. **Unknown MSA database release.** The ColabFold service did not
   expose a database release, and the provenance records it as `unknown`.
   The archived bytes and hash make this run reproducible, but the same
   query today may return a different alignment.
5. **Possible training overlap** (section G).
6. **Limited meaning of the confidence scores** (section E).
7. **No waters or hydrogens** in the prediction or in the contact metric.
   Water-mediated and hydrogen-bond recognition were not assessed.
8. **Incomplete inference isolation and non-deterministic settings**
   (section D).
9. **Assumption reconsidered.** The configuration recommended 24 GB of
   VRAM. A 12 GB laptop GPU sufficed, with 6,900 MiB sampled usage. This
   applies to this 127-residue, single-ligand, single-sample run only.

## I. Reproducibility and evidence register

All files below are committed on `phase3-4men-msa-20261017`. None of the
evidence exists only locally. Model weights and the CCD cache are not
committed. They are identified by immutable Hugging Face revision and
SHA-256 in `runtime.json`.

Commits:
- `c9b09dc`: MSA pinned and bundle frozen; precedes inference.
- `c9b6679`: `metrics.json`.
- `d2e96b8`: raw outputs, processed inputs, log and runtime files.

`metrics.json` was committed about 18 hours before the raw outputs.
Because it records the prediction, confidence and runtime hashes, the
later commit can be checked against it, and every hash matches.

| File (under `EXP/` unless a full path is given) | Purpose | SHA-256 |
|---|---|---|
| [`configs/phase2/4men_primary.json`](../../configs/phase2/4men_primary.json) | Primary design | `735234212086da69b67105136e668cda70fa82dc3d5552bd984c6330e77398c0` |
| [`data/phase2/msa/4men_full.a3m`](../../data/phase2/msa/4men_full.a3m) | Pinned MSA (also copied to `inputs/protein.a3m`) | `8d0374e94722d690f26a96b6539e133f55dfb11cb508083d8c12d3abac2baf53` |
| [`data/phase2/msa/4men_acquisition_20261017.json`](../../data/phase2/msa/4men_acquisition_20261017.json) | MSA acquisition record | `df4d8f6c931487df0ca0240633d5e349308078dc2955d9a471f0c97ed23673f2` |
| [`data/raw/4MEN.pdb`](../../data/raw/4MEN.pdb) | Reference structure | `a1b21bdc30365fb0ca905a6f9f00cc5955e3b382e6b726c02e4a2162be55b6a0` |
| `config.json` | Frozen experiment config | `8c044ee909a1a04b10b11426564a44c338011ec501e0e9b0b90742ffbfb6a165` |
| `preregistration.json` | Metrics, gates, fit/pocket residues | `7ac3d2bc63d30e4187c82ce4ca5ad20901a894ae8ab82f0888a0d14af60450e8` |
| `command.json` | Frozen command | `6a9b5076e2253eb356852cbb0974774fe966a006ef28aefacd0a97c5688d4308` |
| `validation.json` | Offline validation, reference contacts and geometry | `aeecf2e563ebfc9761a7fed165f73700ab4432abd5933414425090a9fbb9e4fe` |
| `manifest.json` | Hashes of frozen files and analysis source | `1dc10a6d798868280137a4663a6d088b9ae63281ebb85796bec04d4ef88c2ed2` |
| `inputs/4men.yaml` | Boltz input | `5667b80740a209ccb099f69ed82fd705f621db9442230b03321309eb68fa868f` |
| `ligand_atom_names.json` | SMILES ↔ CCD atom names | `0dffb08f3edade1f72b871aab96d8bd351a5aa1b434486ce8ced5759fc72d040` |
| `residue_mapping.csv` | Input index ↔ author/UniProt numbering | `b56c707f1b2a11b28ef5a6ae22064287da67f12e44a84d8f1e71af6251cd87d5` |
| `runtime.json` | Pre-inference environment provenance | `2172a3dd4d4139d6a72c2e7a0436aa8a8bb761acf59ebfa5bfd587441fb41c9c` |
| `runtime.sha256` | Receipt for `runtime.json` | `70f0cf3a5e369e9b0d335fee7314fb9a3e5e74f7061a98654e73b5272f405af0` |
| `run_observation.json` | Timing, exit code, GPU, sampled memory | `7a559bb5db6dff18818b4de4357ef9f94cb34b68fb01fcd8c0edc38aa407fdee` |
| `logs/boltz_predict.log` | Inference stdout/stderr | `2803ec5e644b39217453284b53c259852d34d1563b8cd67f2d10f46dcfcbb28b` |
| `outputs/.../predictions/4men/4men_model_0.cif` | Predicted complex | `2d06e475e8a7cd96eb17bd81c8e9b4ce2572d2a3c057fbd1a09318852a1fb37a` |
| `outputs/.../predictions/4men/confidence_4men_model_0.json` | Confidence summary | `8aad7010d2466b9da65f8a10ad11d8d27c7f62fbaec8d0e88b598026510f8483` |
| `outputs/.../predictions/4men/pae_4men_model_0.npz` | PAE array | `1f2e873fb922f9f82dc35772fc41c58f3fd672bd9c57ede7a46dafdb790978de` |
| `outputs/.../predictions/4men/pde_4men_model_0.npz` | PDE array | `d5b89e9ef390437b5e662e77b26291a359a00b1c64d1c14b639cdcc19ed61194` |
| `outputs/.../predictions/4men/plddt_4men_model_0.npz` | pLDDT array | `bb31d174e85561de4a5b97e4329cb12c8a4c411a8c8f2faf459a15a9f235d47d` |
| `outputs/.../processed/manifest.json` | Boltz processed manifest | `3fce071748f5d1a30732350164de67a0eafbc4bfde6490df1d028de40fee9784` |
| `outputs/.../processed/records/4men.json` | Processed record (no templates/constraints) | `04f978919551509b01a40fa814c61acdf6a5e75b13fee949e9f55ba399e21c9d` |
| `outputs/.../processed/mols/4men.pkl` | Prepared ligand as consumed | `822cc49a30d3aee7ac8978f41c9657cbfe290825f7ac2ddbc2ab7fbbbcf303eb` |
| `outputs/.../processed/msa/4men_0.npz` | Processed MSA as consumed | `f8bea2f180ba237e76cc8cb7cd69ca9b3e46a1539b980600efe4992515f763ad` |
| `outputs/.../processed/structures/4men.npz` | Processed input structure tokens | `5b749bbba663eae6c1d44c9e1117637f256f488f76ee79fe98f657266ca4d2b8` |
| `outputs/.../processed/constraints/4men.npz` | Processed constraints | `5bb6207de735409d85bb089dc1e6ec0d8359cd296c8f893bd4856fcb6eda94ce` |
| `outputs/.../lightning_logs/version_0/hparams.yaml` | Model hyperparameters as loaded | `11e4031595a3a6bedeaeb8e0402b3bb278280b072a65e9435e96fefb883863c6` |
| `analysis/metrics.json` | Scientific evaluation | `e0bf5efb83167df4c43ed9457374b545979deed892aaf3fd03c2cae10eb0ea0f` |

`outputs/...` stands for `outputs/boltz_results_4men`. The `processed/`
files include a pickle. Load it only from this trusted, hash-verified
archive.

**Archival gaps:**

- No transcript or output from either external audit is in the
  repository.
- No contemporaneous, signed execution-approval record exists (see the
  execution-provenance statement in section D).
- The ColabFold database release is unknown.
- No test in the repository reproduces the 12% expansion result. This is
  planned as a Phase 4 unit test.

## J. Conclusions and next steps

**1. What this experiment establishes:**

- The prepared pipeline runs end to end on a consumer GPU in about
  1.5 minutes: pinned inputs, verified model and cache, Boltz 2.2.1
  inference, and analysis in a separate environment.
- Its provenance chain links the frozen inputs, runtime, outputs and
  metrics by hash.
- On this one known complex, the single preregistered sample reproduces
  the deposited ligand placement after the protein fit (0.14 Å) and the
  4.0 Å residue contact set, with well-formed heavy-atom ligand geometry.
- The analysis code produces sensible values on a real Boltz output: atom
  naming, symmetry handling, partial-residue coverage and contacts all
  work.

**2. What remains unverified:**

- Accuracy on complexes outside the training data.
- Seed-to-seed and run-to-run variability.
- Atom-level and hydrogen-bond fidelity.
- Water-mediated recognition.
- The calibration of the confidence scores.
- Whether 4MEN was in this checkpoint's training data.
- Bitwise reproducibility of a rerun.
- Anything about affinity.

**3. Justified follow-up validations.** Each needs its own experiment ID
and its own preregistration:

- A multi-seed replicate of identical inputs, with the number of seeds
  and the summary statistic fixed in advance.
- Validator sensitivity tests: committed unit tests with synthetic
  perturbations, starting with a reproducible version of the audit's 12%
  centroid expansion and extending to torsion changes, translations and
  reflections, plus scale-sensitive or PoseBusters-style checks. These
  need no GPU. Any stricter gates apply to future experiments only, not
  to this one.
- Atom-level interaction analysis of the existing prediction, for example
  per-residue minimum distances and the Asn140 contact distance. This
  would be a new, separately recorded analysis that leaves
  `metrics.json` unchanged.
- Evaluation on complexes released after the cutoff, with a recorded
  audit of sequence, pocket and chemotype neighbours. 32SC/9FWX must be
  counted once, at its 2024 date.
- A properly preregistered MSA ablation. The v4 no-MSA smoke bundle is
  not an ablation.

**4. Premature claims:**

- "Boltz-2 accurately predicts BRD4 ligand poses."
- "This validates co-folding for drug discovery."
- "The confidence score shows the pose is reliable."
- "The model has learned the physics of KAc recognition."
- "The model memorised 4MEN."
- "Hydrogen bonds were reproduced."
- Any affinity statement, or any claim about 4MEO.

**5. Relationship to Phase 4.** Phase 4 is in offline development on a
separate branch. **No further GPU inference has been approved.** Offline
work, such as validator unit tests and leakage audits, can proceed
without it. Any experiment that needs inference (multiple seeds, held-out
or perturbed cases) needs its own preregistration and its own approval
first. Phase 4 aims to study robustness and possible memorisation. It will use new experiment IDs and
output directories. It must not edit, rerun into or reanalyse the
`4men_msa_s20261019_v1` bundle, or alter `metrics.json` or the
preregistered thresholds. Phase 4 results will be reported alongside
this result, not in place of it. The design of Phase 4 is not described
or endorsed here.

### Interview summary

> "I built a preregistered, hash-pinned pipeline around Boltz-2 and ran
> it on BRD4 BD1 with compound 3 from Vidler et al., PDB 4MEN, as a
> control. I fixed the construct, the MSA, the seed, the alignment and
> the acceptance gates before running anything, and scored the only
> sample without choosing it against the crystal structure. It
> reproduced the ligand position to 0.14 Å after a protein fit and
> recovered all 11 contacting residues at 4 Å. Even so, I treat it as
> evidence that the pipeline works, not that the model does. 4MEN was
> public ten years before Boltz-2's structural training cutoff. It is
> one case with one seed, and at 3.5 Å contact recall falls to 0.4,
> so residue-level agreement does not mean every interaction is exact.
> An independent audit also showed that my geometry screen would pass a
> ligand inflated by 12%, so the validator needs strengthening before I
> rely on it for harder cases. The informative next experiments are
> multi-seed runs and held-out complexes released after the cutoff,
> with a recorded audit of their neighbours."
