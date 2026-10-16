# Controlled Boltz-2 pilot: reconciled preregistration

Status: **design reviewed against all six supplied deliverables; no
prediction executed**. The primary scientific condition follows the
supplied protocol's Arm A: sequence + full cached MSA, no template.
Its first GPU run is one seed and one sample for timing and structural
evaluation. The no-MSA condition is an optional smoke test only.

The primary design is
[`configs/phase2/4men_primary.json`](../../configs/phase2/4men_primary.json).
Its real MSA and provenance have not been acquired. Null fields mean
**pending**, not a fabricated alignment or a valid checksum. `plan`
validates the design offline; strict `validate`, `prepare` and subsequent
runtime checks prevent executing an unpinned input.

The separate smoke configuration is
[`configs/phase2/4men_smoke.json`](../../configs/phase2/4men_smoke.json),
with prepared bundle `experiments/phase2/4men_smoke_s20261019_v4/`.
The earlier `4men_single_sequence_s20261019`, v2 and v3 smoke bundles are
unchanged, superseded preparation records. V4 freezes the corrected MSA,
CCD-cache and runtime-provenance gates. No bundle contains model predictions.
See [the reconciliation](PROTOCOL_RECONCILIATION.md) for every material
departure from the supplied protocol. The library has no model-download,
MSA-search, prediction-launch or seed-campaign command.

## Software and interface verification

PyPI and the upstream release identify **boltz 2.2.1**, published
8 September 2025, as the release checked on 15 October 2026. The small wheel
was downloaded into memory, its SHA-256 verified against PyPI, and selected
source files extracted for inspection. Boltz and Torch are not installed
in the Phase 1 environment. No `boltz` process has been invoked.

| Item | Pin |
|---|---|
| PyPI wheel SHA-256 | `b8c62bbdede1922931d9203118f62c858f11aa699bf91fd4c05a5ed6a6d8b4fc` |
| Upstream tag | `v2.2.1`, commit `cb04aeccdd480fd4db707f0bbafde538397fa2ac` |
| Model repository | `boltz-community/boltz-2` |
| Model revision | `6fdef46d763fee7fbb83ca5501ccceff43b85607` |
| Structure checkpoint SHA-256 | `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1` |
| RDKit for SMILES atom naming | `2026.3.6` |

CLI options are checked against the wheel's Click decorators, including the
mixed spelling `--preprocessing-threads` and `--recycling_steps`.
YAML fields and output paths were checked against the wheel parser/writer.
This is static interface verification; compatibility with an actual GPU
runtime and real emitted output remains an execution-stage check.
The pinned official GitHub prediction guide and parser/CLI source were
re-fetched during this follow-up and match the saved guide/wheel extracts.
The guide's options table still shows 1.638 for `step_scale`; the Boltz-2
branch in the actual v2.2.1 code uses 1.5. We specify 1.5 explicitly.

Boltz requires, among other dependencies, NumPy `<2.0`, Biopython `1.84`
and requests `2.32.3`. These conflict with the pinned Phase 1 environment.
Use a separate Linux/Python 3.12 environment; do not install Boltz into
`.venv`. Pin Boltz/RDKit as above, resolve Torch/CUDA against the selected
host, and freeze the complete resolved environment **before prediction**.
The host image, Torch/CUDA build and driver are still unselected.
Environment setup is part of the future approval, not completed validation.

The source unconditionally calls `download_boltz2`, which checks
`mols.tar`, its extracted `mols/` directory, `boltz2_conf.ckpt` **and**
`boltz2_aff.ckpt`, even when affinity is not requested. The immutable URLs,
expected sizes and hashes for all three files are recorded in the prepared
bundle. A cached file must be verified before invoking Boltz; an existing
filename alone is insufficient. No model asset was downloaded here.

## Construct and ligand

Use one protein chain A with the full deposited 127-residue sequence:
the `SM` expression remnant at author positions 42–43, followed by UniProt
**O60885 residues 44–168**. RCSB polymer metadata, the PDB `SEQRES`,
`DBREF`/`SEQADV` annotations and the observed residue sequence agree.
All 127 Cα positions are present. The 125 native positions are the fit set.
The generated residue mapping records input index, author chain/number,
insertion code, residue identity and UniProt position; the cloning remnant
has no UniProt position. Landmark checks include W81, D88, Y97 and N140.

Use ligand chain B with SMILES:

```text
Cc1ccc(CN(C)c2cc(C)nc3ncnn23)cc1
```

RDKit confirms the same full graph as the pinned 25K SDF and CCD mmCIF:
C15H17N5, 20 heavy atoms, formal charge zero, no stereocentres or specified
stereobonds, InChIKey `JAEVUMCKPBRDAJ-UHFFFAOYSA-N`. The alternate SMILES
in the supplied evidence brief gives the same graph. This is template
chemistry consistency, not a new experimental sample identification.
The neutral input is a fixed modelling assumption consistent with the
CCD; it is not an experimentally established pKa or proof that other
protomers cannot matter.

The generated YAML supplies only sequence, SMILES and an explicit MSA
policy. It supplies no receptor coordinates, templates, binding-site
constraints or affinity property. Reference coordinates are used only for
evaluation. SMILES input does not eliminate training leakage.

## One-run design and MSA policy

| Setting | Preregistered value |
|---|---|
| Purpose | First Arm A timing/structural evaluation on a known complex |
| Model | `boltz2` |
| Seed | `20261019` |
| Diffusion samples / parallel samples | 1 / 1 |
| Recycling / diffusion sampling steps | 3 / 200 |
| Step scale | 1.5 |
| MSA sequence limit | 8192, upstream default; no `--subsample_msa` depth ablation |
| Steering potentials | Off |
| Optimized kernels | Disabled with `--no_kernels` for a simpler initial runtime |
| Templates / constraints / affinity inference | None |
| Output | mmCIF plus the model's confidence JSON |
| Sample selection | The sole `model_0`; no reranking against the reference |
| MSA | Full cached A3M, supplied as `msa: inputs/protein.a3m`; actual bytes/hash/provenance required before preparation |

The supplied protocol's `--use_msa_server` syntax is valid when `msa` is
omitted. For controlled reuse, acquire that alignment once after approval,
archive it, and bind its checksum and provenance **before** freezing the
primary experiment. No MSA service was called in this preparation.
The proposed service is `https://api.colabfold.com`; record its actual
generator/version when available, database identity/release, retrieval
time, query, returned row count and file hash. Use an explicit `unknown`
when the service does not expose a database release; do not invent it.
Record acquisition duration separately from inference duration.

The primary inference command consumes the local file and does not use
`--use_msa_server`. “Full MSA” means the obtained alignment with the
upstream 8192-row limit and no intentional depth ablation, not every
homolog that might exist. Save the complete returned file even if the
model caps the rows it consumes. Reuse exactly these bytes in later
approved arms; a changed MSA requires a new experiment ID.

The validator checks the full 127-residue query, aligned row lengths and
lowercase-insertion handling. At least one non-query record in the first
8192 records must contain an aligned standard amino-acid symbol. Gap-only,
unknown-only and lowercase-insertion-only padding cannot establish a hit.
Padding rows may coexist with an informative row; their bytes are not
rewritten. This is a structural input check, not independently verified
homology: the real source and provenance must still be reviewed.

The no-MSA smoke uses explicit `msa: empty` and otherwise the same seed,
model and sampling settings. It is labelled `smoke_test` in validation,
preregistration and analysis. Its result does not substitute for Arm A or
estimate the effect of MSA removal. Upstream warns of reduced accuracy;
neither unchanged peak VRAM nor a fixed runtime ratio is assumed.

## Metrics fixed before looking at predictions

**Protein-aligned, symmetry-aware ligand RMSD.** Require the entire
reference ligand graph and all expected protein heavy atoms except
optional `OXT` at the final residue. Boltz 2.2.1's canonical standard
residue atom set omits this terminal oxygen; 4MEN models it on GLU168.
Any omission is explicitly reported in `atom_coverage`. All Cα, backbone
`O` and other expected heavy atoms remain required. The reference OXT is
outside the 8 Å pocket and all registered contact cutoffs; reference
coordinates and metrics are retained.

The reference also lacks three distal side-chain atoms at Lys91 and three
at Gln123. A complete prediction can legitimately include these atoms.
They are accepted only if their names and elements match a residue-type
template assembled from other occurrences in the pinned reference, and
are listed in `atom_coverage`. That template inventory was independently
checked against Boltz 2.2.1's canonical atom lists. No missing experimental
coordinates are invented. These atom-name/element templates are used for
evaluation only; the model input still contains no structural template or
reference coordinates. Align the
predicted protein onto the reference using an unweighted proper Kabsch
rotation/translation of Cα atoms at author residues 44–168. Reflections,
missing residues, guessed correspondences and fit-set selection are
forbidden. Apply this transform to the ligand.

For each full heavy-atom graph isomorphism preserving elements, bond
orders, charges, isotopes and specified stereochemistry, compute
`sqrt(mean_i ||x_reference[i] - x_aligned_prediction[mapping[i]]||²)`.
Report the minimum, the chosen mapping and the number of mappings
evaluated. **Do not fit the ligand separately.** Abort if enumeration
exceeds 4,096 mappings rather than reporting a truncated search.
Boltz's SMILES atom names are reproduced from its versioned parser and
checked in the runtime environment; atom-order coincidence and
distance-inferred ligand bonds are not used.

Report protein Cα RMSD and Cα RMSD for reference residues within 8 Å of
the ligand as diagnostics. The pocket diagnostic uses the same global
protein fit, without a second fit or selection using the prediction.
Reference protein alternate conformers follow Phase 1's Biopython
highest-occupancy selection, with file order breaking ties. Predicted
alternate conformers or partial occupancy are rejected.
This deliberately retains the existing 125-native-residue fit and
8 Å pocket diagnostic, rather than adding the draft's all-127 fit and a
second 4.2 Å pocket fit. The difference is declared before predictions;
there is no choice of the best-fitting metric after seeing an output.

**Contact precision and recall.** Contacts are sets of unique mapped
polymer residues, each having at least one heavy atom within the inclusive
cutoff of any ligand heavy atom. Evaluate the reference and prediction in
their own respective complex frames, at full coordinate precision.
Primary cutoff: 4.0 Å; 3.5 and 5.0 Å are prespecified sensitivity analyses.
Waters, additives, symmetry-mate crystal contacts and inferred interaction
types are excluded.

For reference set R and predicted set P:
`precision = |P ∩ R| / |P|`, `recall = |P ∩ R| / |R|`.
An empty predicted set against a nonempty reference scores zero for both.
An empty reference is unscorable; empty/empty is never a perfect score.
This is residue-level contact recovery, not atom-pair or hydrogen-bond
recovery. The 4MEN reference contains 5 / 11 / 13 residues at 3.5 / 4 / 5 Å.

The primary contact metric retains **all predicted protein heavy atoms**,
including valid atoms at partially resolved reference side chains; its
definition and thresholds have not changed. A separate
`contacts_reference_observed_atoms` diagnostic restricts the prediction
to atoms with reference coordinates. Interpret any difference in the
context of the unresolved reference side chains. A primary “false
positive” means absent from the deposited-model contact set, not proof
of a physically wrong interaction. Geometry/clash checks still use all
predicted protein heavy atoms. The diagnostic never selects or replaces
the primary acceptance score.

**Geometry screen.** Check finite coordinates, exact atom coverage and
graph/element consistency before metrics. Evaluate bond lengths and
heavy-atom bond angles against the pinned CCD ideal geometry, aromatic
ring planarity, intraligand nonbonded overlap, and protein–ligand severe
overlap. No ligand relaxation, minimization or post hoc repair is allowed.
25K is achiral; the utility rejects ligands needing stereochemical
geometry validation rather than pretending to support them.

| Primary acceptance gate | Threshold |
|---|---|
| Ligand RMSD after the protein fit | ≤ 2.0 Å |
| Protein Cα RMSD on the full fit set | ≤ 2.0 Å |
| Residue contact precision at 4.0 Å | ≥ 0.70 |
| Residue contact recall at 4.0 Å | ≥ 0.70 |
| Absolute bond-length deviation from CCD ideal | ≤ 0.25 Å for every heavy-atom bond |
| Absolute heavy-atom bond-angle deviation | ≤ 30° for every angle; degenerate angles fail |
| Aromatic ring distance from best-fit plane | ≤ 0.15 Å maximum |
| Intraligand nonbonded distance / Bondi radius sum | ≥ 0.75; exclude bonded and 1–3 neighbours |
| Protein–ligand Bondi overlap | No heavy-atom pair with overlap ≥ 0.80 Å |

All gates must pass for a **harness success**. These are operational pilot
thresholds, not calibrated probabilities or a validated screen for drug
activity. The geometry screen is deliberately limited: it is not a
force-field energy, complete PoseBusters validation, hydrogen-bond
analysis, or proof of physical stability. The deposited 4MEN reference
passes this screen before predictions.

Confidence fields (`confidence_score`, `ptm`, `iptm`, `ligand_iptm`,
`complex_plddt` and associated output) are saved as diagnostics. They do
not override failed geometry/pose metrics. One case cannot calibrate
confidence or quantify seed variance.

## Leakage register and later case

4MEN/4MEO predate the report's structural cutoff. The same study's 4MEP
and 4MEQ, docking/template structures 3MXF and 2OSS, other BRD4 complexes,
other BET bromodomains and close ligand analogues create additional risks.
Structural, distilled-prediction, CCD/chemical-prior, affinity-label and
MSA exposures must be considered separately. No exact checkpoint training
membership has been established here.

32SC is excluded as an independent 2026 case: 9FWX already exposed the
same compared coordinates in 2024. The
[assessment](../../data/phase2/review/32SC_assessment.json) records the
release dates, exact-graph search, coordinate comparison and a citation
metadata conflict. A future evaluation would deduplicate these entries,
use the earliest public complex date, audit sequence/pocket relatives and
ligand/scaffold neighbours across pre-cutoff structures and labels, and
report unresolved training coverage. It would remain a retrospective
evaluation. No 32SC or 4MEO prediction is included in this pilot.

## Provenance, directories and execution boundary

Preparation writes a new experiment directory and refuses to reuse it.
The frozen manifest hashes inputs, metrics, atom/residue maps and source
snapshots. New analysis code and reused Phase 1 utilities are hashed with
source line endings normalized to LF; scientific data and prepared
artefacts retain byte-exact hashes. Changing a frozen file or the analysis
code prevents scoring until a new reviewed experiment is prepared.

After user approval, run `capture-runtime` in the separate Linux GPU
environment before invoking Boltz. It checks Boltz/RDKit versions and
atom names and checks relevant package versions against the pinned Boltz
wheel requirements. It verifies cached model assets and compares every
extracted CCD filename and file byte with the verified `mols.tar`.
Missing, altered or unexpected files/directories, archive path traversal,
absolute/drive paths, links, special files and duplicate normalized archive
entries fail before any GPU check. The comparison reads the archive; it
does not extract it or load its pickle contents.

Capture records a versioned `runtime.json` plus a `runtime.sha256`
companion receipt, including the prepared-input hash, actual copied-MSA
hash, model revision/assets, verified CCD inventory digest and count,
installed distributions, Python, OS, Torch/CUDA, cuDNN, driver/GPU identity
and relevant environment settings. Both files are written exclusively.
It does not download weights or invoke Boltz.

From the **prepared experiment directory**, the planned model command is:

```bash
boltz predict inputs/4men.yaml \
  --out_dir outputs --cache /workspace/boltz-cache-v2.2.1 \
  --checkpoint /workspace/boltz-cache-v2.2.1/boltz2_conf.ckpt \
  --model boltz2 --accelerator gpu --devices 1 --output_format mmcif \
  --num_workers 0 --preprocessing-threads 1 --seed 20261019 \
  --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 \
  --max_parallel_samples 1 --step_scale 1.5 --max_msa_seqs 8192 --no_kernels
```

This command is a specification, **not permission to execute it**.
On an approved host, wrap it with the scheduler's 30-minute termination
limit (or GNU `timeout --signal=TERM --kill-after=60s 30m`) and capture
stdout/stderr, exit status, wall time and peak GPU memory. Cache verification
must succeed first; block inference-time network access so an unexpected
missing asset cannot trigger a download. Do not use `--override`.

Expected outputs are under
`outputs/boltz_results_4men/predictions/4men/`, including
`4men_model_0.cif` and `confidence_4men_model_0.json`. Retain processed
inputs and logs with that experiment. Record start/end UTC, actual
hardware, billed duration/cost and peak memory in `run_observation.json`;
these are presently unmeasured. A seed alone does not ensure identical
GPU results across software/hardware.

Then use the Phase 1 analysis environment to run:

```text
.venv/Scripts/python.exe -m brd4_analysis.phase2.cli analyse experiments/phase2/4men_msa_s20261019_v1
```

Analysis requires the frozen bundle, `runtime.json`, its `runtime.sha256`
receipt and exactly the sole planned model. Before scientific scoring it
checks the receipt, the complete runtime schema, and consistency with the
frozen config/input, command, MSA checksum/provenance, model revision and
checkpoint/archive metadata. Runtime package pins/requirements and CUDA
build metadata must also agree. Copy both runtime files with the outputs
when returning to the analysis environment; model weights need not be
copied there for manifest validation.

The receipt detects changes to a captured record; it is not cryptographic
attestation against deliberately forging the record and receipt together,
nor proof that a subsequent manual command used the recorded environment.
Execution logs and operator provenance remain necessary.

Analysis writes `analysis/metrics.json` once, including input,
prediction, confidence and runtime hashes. Stop after this one run and
review it, regardless of whether its acceptance gates pass.
