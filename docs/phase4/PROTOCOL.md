# Phase 4 experimental protocols — DRAFTED, NOT EXECUTED

These protocols are design documents for review. None of them has been
run; running any of them requires explicit approval, new preregistration,
and (for Boltz cases) the GPU runtime — all outside the bounds of this
Arm-A implementation task. Each protocol states its question, controls,
registered analysis and decision criteria *before* any execution.

Conventions carried forward: frozen acceptance thresholds live in the
experiment config (never in code); every run records seed, checkpoint
SHA-256, input SHA-256, environment provenance; no prediction is scored
against more than the single preregistered reference; "success" language
follows `preregistration.json` (pipeline control, not generalisation
evidence).

---

## Protocol P4-1 — multi-seed variability on the original complex

**Question.** Is the Phase 3 reconstruction (ligand RMSD 0.1405 A) a stable
property of the model+MSA+input, or a lucky draw from the diffusion
sampler?

**Design.** Repeat `boltz predict` on the *unchanged* primary inputs
(`inputs/4men.yaml`, pinned `4men_full.a3m`) at `K = 8` prespecified seeds
distinct from 20261019, all other flags frozen (recycling 3, sampling 200,
diffusion_samples 1, step_scale 1.5, max_msa_seqs 8192, no templates,
no potentials, boltz2_conf.ckpt pinned by SHA-256). Same GPU class;
CPU-feasible alternative: none (sampling is GPU-bound — flagged as a hard
dependency).

**Registered analysis.** For each seed: the five frozen gates on model_0,
plus *stability statistics* — distribution of ligand RMSD, contact-set
Jaccard across seeds, per-residue contact distance variance, and the
Phase-4 candidate checks (similarity scale, internal distance RMSD, bond
bias) as diagnostics. Multi-seed agreement measured by pairwise
protein-aligned ligand RMSD between seeds.

**Controls.** The committed prediction (seed 20261019) is the reference
replicate; `control_reference_replay` gates the scorer each run.
**Decision criteria (proposed).** Stable if ≥ 7/8 seeds pass all five
gates AND pairwise seed ligand RMSD median ≤ 1.0 A. Interpretation is
sampler variability only — not accuracy.

---

## Protocol P4-2 — chemical analogue perturbations

**Question.** Do the model+evaluator distinguish near-neighbour chemotypes of
25K at the structure-reconstruction level, and do the geometry gates reject
implausible analogue poses? This is not an affinity-ranking or prospective
predictive-accuracy protocol.

**Design.** A prespecified panel of ≤ 6 CCD-deposited BET-fragment
analogues sharing the triazolo-pyrimidine/tolyl scaffold (e.g. 25V-type
N-oxide, des-methyl, regioisomeric methyl shifts), each with: its own CCD
ideal geometry + a *registered* experimental complex reference if one
exists in the PDB, else it is an unlabelled probe. For each analogue, one
seed, one sample — same frozen command shape. The adaptor
(`read_prediction`) must be generalised per-ligand (component ID, atom
names, graph) and that generalisation is itself registered and hashed
before execution.

**Registered analysis.** The five gates where a deposited complex of the
same analogue exists (e.g. 4MEO/25V with the declared publication/MS
uncertainty); otherwise graph-consistency diagnostics only (does the
predicted ligand's own internal geometry pass the screen — not "is it
correct"). Ligand-identity checks (InChIKey, element multiset) are scored
as the *input-validation* layer they are.

**Controls.** 25K on 4MEN re-run in the same session; a deliberately
mismatched analogue (wrong scaffold) as a negative control expected to be
rejected by graph checks. That negative control tests ligand-identity
rejection, not binding discrimination. **Note:** chemically realistic
substitutions are a different question from coordinate-level adversaries —
a methyl→Cl variant changes the graph and is *supposed* to be caught by
identity checks if scored against the wrong reference.

---

## Protocol P4-3 — BRD4 pocket mutations

**Question.** Is pose quality robust to conservative pocket-proximal
mutations, and can the contact metrics detect a mutation-induced binding
mode change?

**Design.** Prespecified conservative substitutions on pocket residues
that preserve the BET fold — e.g. the ZA-channel Trp-shelf context and the
Asn-140 anchor neighbourhood — chosen by conservation analysis, NOT by
energetics prediction. Options: N140A (destroys the KAc-mimetic H-bond —
positive control for a *changed* contact set), W81F (conservative),
V87A (cavity), L94V. Each mutant: one seed, one sample, same MSA policy
(wild-type MSA; the mutation enters only the input sequence — declared
limitation, since the MSA won't reflect the variant).

**Registered analysis.** Gates evaluated *against the wild-type 4MEN
reference* — explicitly NOT as "did it still bind" but "does the model
keep or change the pose when the sequence changes". Mutant-vs-wildtype
contact-set difference and pose delta are the readouts. N140A is expected
to lose the ASN140 contact — a positive control that the metric can see
it.

**Controls.** Wild-type re-run in the same session; a far-surface mutation
(e.g. G108A) as a negative control expected to preserve the pose.
**Limitations declared up front:** single-point mutants at diffusion
sampling level are noisy; any apparent "mutation sensitivity" needs the
multi-seed protocol (P4-1) run on the mutant too.

---

## Protocol P4-4 — post-training-cutoff evaluation cases

**Question.** Does the reconstruction result hold for BRD4-BD1 ligand
complexes publicly deposited *after* the Boltz-2.2.1 training cutoff?
Post-cutoff deposition reduces direct memorisation risk, but checkpoint
membership remains unverifiable; it does not prove the case was unseen.

**Design.** Candidate complexes: BRD4 BD1 structures deposited after the
verified checkpoint training-data cutoff, screened by prespecified
overlap controls before inclusion:
- **sequence control:** BD1 sequence identity to the 4MEN construct ≥ 95 %
  (same domain) or explicitly reported otherwise;
- **pocket control:** pocket-residue identity vs the 33-residue registered
  pocket; report any mutations explicitly;
- **chemotype control:** ligand Morgan Tanimoto to 25K reported, PLUS max
  similarity to a pre-cutoff BET-fragment corpus — pairwise-to-25K alone
  is insufficient (the exact limitation recorded for the 32SC/9FWX case
  in `phase2/audit.py`, whose deposited coordinates were found identical
  between two entries and whose earliest public date — 2024-08-07 — postdates
  the report's own structural cutoff anyway);
- **model-knowledge control:** checkpoint membership remains unverifiable;
  the honest claim is "post-public-deposition", not "definitely unseen".

**Registered analysis.** Identical five gates on each case + candidate
checks; failure modes classified (pose error vs internal distortion vs
contact-set shift). At least one structural outlier (e.g. open/closed WPF
shelf state) included as a documented hard case.

**Gate:** candidate set frozen and hashed before any inference. Estimated
2–4 qualifying complexes at current deposition rates.

---

## Protocol P4-5 — independent comparison with classical docking

**Question.** Is the Boltz reconstruction better than a classical
structure-based method at the same task — the calibrated baseline every
ML claim needs?

**Design.** Dock 25K into the prepared 4MEN BD1 pocket (apo reference or
holoreference — declared which; self-docking into the holo pocket is the
minimal version) with an established CPU-only engine (e.g. AutoDock Vina /
Smina — exact tool + version pinned at preregistration, checked for
license). Exhaustiveness, box centre/size and seed prespecified. Up to 10
poses per run.

**Registered analysis.** Score each docked pose with the SAME five frozen
gates (poses converted into the evaluator's expected frame — the
conversion tooling is Phase-4 code, reviewed). Compare best-pose and
top-ranked-pose ligand RMSD and contact P/R against the Boltz results.
Energy scores reported as diagnostics only, never as gates.

**Controls.** Redocking the deposited 25K pose is the positive control;
a different-graph decoy ligand is a negative control for identity/pipeline
handling, not a validated test of binding discrimination.
**Interpretation guardrail:** single-target docking comparison supports
only "on this complex, at this preparation" claims.

---

## Cross-cutting requirements for all protocols

- Each protocol gets its own `preregistration.json` + manifest hashing,
  mirroring `phase2/workflow.prepare`.
- New analysis code lands in `brd4_analysis.phase4` (or a new isolated
  package), never inside `phase2` — the frozen `analysis_source_hashes`
  contract is what keeps executed bundles verifiable.
- GPU inference, checkpoint download and new dataset acquisition are
  out of scope for Arm A and require the existing user-approval gate.
- Reporting language: "pipeline control under declared experiment role",
  never clinical/predictive accuracy claims.
