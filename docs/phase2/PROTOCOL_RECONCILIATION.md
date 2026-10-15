# Reconciliation of the two Downloads deliverables

Both files were located and read in full on 15 October 2026. Their actual
Downloads filenames omit the `_claude_science` suffix requested in the
follow-up. The archived copies have that suffix and preserve their bytes.

| Located in `C:\Users\LLR User\Downloads\` | Repository copy | SHA-256 |
|---|---|---|
| `03_boltz2_protocol_4MEN.md` | `docs/protocols/03_boltz2_protocol_4MEN_claude_science.md` | `749e90b04a981602b28c21353f20b1fdeb43a3cf4ea33bcdf31c8f0924747269` |
| `04_interview_questions.md` | `docs/research/04_interview_questions_claude_science.md` | `ce9cead585ffe35540dc06c60b18c5732e84370fcd431f83d248316fe3f8ebd1` |

The files are source material, not automatically verified conclusions.
The first four deliverables and the prior review remain preserved.
The [integration memo](../research/REVIEWED_INTEGRATION.md) applies to all
six, including the reused map, sample and training claims in the questions.

## Protocol decisions

| Draft provision | Reconciled implementation or disposition |
|---|---|
| Arm A: no template, full MSA | **Primary scientific condition.** Its minimal first run is one seed and one sample, as requested in draft §11–12. Cache and hash the real alignment before freezing the run. |
| Arm B: `msa: empty` | Retained only as an optional **smoke test** at this stage. It cannot stand in for Arm A or establish MSA dependence by itself. |
| Live `--use_msa_server` | Valid upstream syntax, but inference reuses the acquired local A3M. No search or server request has been performed; its hash and actual provenance remain pending. |
| Five samples and five seeds | Deferred campaign, requiring separate approval. Five samples within one seeded invocation are not five independently seeded runs. |
| Include SM and use the given 127-residue sequence | Adopted; verified against the deposition. Keeping SM still requires a numbering table: Boltz input/output sequence positions begin at 1; reference author numbering begins at 42. |
| Neutral 25K by SMILES | Adopted and graph-checked. “No ionisable centre near physiological pH” is not independently established and is not an input-validation fact. |
| Ligand chain L / PDB output | Use B and mmCIF, both officially supported. The adapter explicitly checks the chain and residue namespaces; the choice is recorded. |
| Unpinned CLI / hosted `boltz-2.1` | Pin the open-source package to 2.2.1 and the model asset hashes. Hosted API settings are not assumed interchangeable and are outside this pilot. |
| Affinity property in the sample YAML | Omitted. Compound 3 has no quantitative IC50 comparator established here, and affinity introduces a different task and leakage risks. The CLI still checks/downloads the affinity checkpoint unless it is pre-cached. |
| All 127 Cα fit, plus separate 4.2 Å pocket fit | Retain the implemented 125-native-Cα primary fit, excluding expression tags. Report the 8 Å pocket Cα error in that same frame. A second fit is deferred; do not select the more flattering frame after prediction. |
| `GetBestRMS` with a warning about fitting | Use explicit complete graph isomorphisms and in-place RMSD after the protein transform. Pre-transforming coordinates does **not** stop `GetBestRMS` from performing a further ligand fit. |
| `<2 Å`, 2–5 Å, >5 Å site categories | Retain the preregistered inclusive ≤2 Å pose gate. Higher RMSD is a pose error, not by itself proof of a different binding site; inspect contacts and placement. |
| Contacts at <4.2 Å and the N2–Asn140 anchor | Retain Phase 1's inclusive 4.0 Å primary residue metric and 3.5/5.0 Å sensitivity metrics. The draft's 4.2 Å counts and per-atom anchor are separate geometric diagnostics, not silently substituted metrics or verified H bonds. |
| Validity alongside RMSD | Adopted as the documented limited heavy-atom geometry screen. It is not claimed to be the full PoseBusters implementation or an energy calculation. |
| Template arm C with 3MXF | Deferred. A co-folder with a protein template is not equivalent to Glide using a fixed holo receptor and explicit waters; the contrast does not isolate a universal “cost of no holo structure.” |
| Arm D: 32SC as a clean held-out case | Corrected. 32SC and 9FWX are **one duplicated public complex**, first established here as public in 2024, not two independent observations or a new 2026 holdout. |
| Confident training membership / guaranteed pocket recall | Replaced with date eligibility and unresolved direct/related-complex, distillation and affinity-label exposure. No exact checkpoint membership manifest was checked. |
| Confidence thresholds and seed correlation | Diagnostic only for this one run. No calibration or variance estimate is possible from one observation. |
| No-MSA has unchanged VRAM; no-kernels is 2× slower and numerically identical | Not accepted as verified. MSA depth can affect feature/memory requirements; kernel/runtime differences need measurement on the pinned host. |
| Self-hosted has no job cost; fixed MSA wait time | Corrected to an unmeasured resource/cost envelope. Infrastructure can be charged, and MSA service queue time varies. |
| DMSO in the 4MEN coordinate environment | Phase 1 identifies EDO and waters in 4MEN; DMS is deposited in 4MEO. This does not establish the complete crystallisation solution composition. |

The primary condition remains a retrospective reproduction of a known
complex. Calling it “primary” selects the main experimental condition; it
does not remove leakage or license a generalisation claim.

## Official syntax check

The immutable upstream source at commit
`cb04aeccdd480fd4db707f0bbafde538397fa2ac` was re-fetched and matches the
saved 2.2.1 wheel's `main.py`, `parse/schema.py` and `parse/yaml.py`.
The corresponding official `docs/prediction.md` also matches the saved
copy. These establish `boltz predict PATH`, local `msa: FILE`, explicit
`msa: empty`, optional `--use_msa_server`, `--seed`, `--max_msa_seqs`,
sampling controls and mmCIF output.

For Boltz-2, source code sets the default step scale to 1.5; the generic
options table's 1.638 also appears as the Boltz-1 default. The generated
command explicitly pins 1.5. The source checks model assets before full
input processing, which is another reason not to invoke prediction as an
MSA-only or validation command.

No Boltz package, model, MSA search or GPU process was launched to perform
this check. Runtime dependency resolution and actual emitted prediction
integration remain gates for a separately approved host.

## Interview document

The five questions were retained; corrected answer guidance is in
[`04_interview_questions_reviewed.md`](../research/04_interview_questions_reviewed.md).
The changes remove claims of proven training membership, a clean 32SC
holdout, model maps decisively settling identity, a universal ΔR-free
threshold and preferential impurity selection as the “best explanation.”
They also distinguish published Vidler results from the work performed
here, the smoke test from the full-MSA primary condition, and a limited
geometry screen from physical validation.
