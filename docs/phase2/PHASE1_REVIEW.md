# Phase 1 merge recommendation

**The Phase 1 correction is ready locally; do not merge without user
approval and native Linux CI.** The scientific framing is acceptable as a
baseline of deposited models. The explicit CSV fix preserves the existing
numeric results and file bytes.

Reviewed PR:
<https://github.com/pkaysantana/BRD4-structure-guided-discovery/pull/1>,
remote head `5c06db21811a29e9925ab008ca555129c71c0a35` at the recorded
review. The separately prepared local commit is **`cba7549`** on
`phase-1-baseline`, one commit ahead of its origin tracking branch.
It was created in `C:\Users\LLR User\BRD4\BRD4-phase1-csv-portability`.
No push, PR update or merge was performed.

## Evidence and remaining gate

| Review item | Finding |
|---|---|
| Earlier review: template-assigned protonation | Addressed: the report and code state graph consistency and its limits |
| Earlier review: stale ECFP4 similarity | Addressed: heavy-atom Morgan value is 6/72; regression and report consistency tests cover it |
| 25V atom roles | Corrected and tested at this head |
| 4MEO identity | Properly unresolved in the baseline; new research does not settle it |
| Corrected Windows baseline | All 45 existing tests pass against the worktree's actual corrected source |
| CSV line endings | Fixed locally: all five `to_csv` call sites explicitly use `lineterminator="\r\n"`; all eight committed CSVs remain unchanged |
| Actual writer under LF-default simulation | All eight CSVs and summary JSON reproduce their bytes; temporary report rendering matches the committed LF form and original workspace |
| Native Linux validation | Pending; WSL is not installed and Docker is unavailable |
| CI | Windows 2022 / Ubuntu 24.04 matrix configured, Python 3.12.10, official actions pinned by SHA; not run or dispatched |

Before this correction, pandas' default LF on Linux made all eight CSVs differ from the committed
bytes even with unchanged numerical results. `.gitattributes` correctly
preserves the files, but cannot fix a writer that varies by operating
system. The existing determinism test therefore exposes a portability
failure on Linux. The actual corrected writer was exercised under both
default-EOL policies, without the previous `DataFrame.to_csv` wrapper.

[`phase1_csv_fix.patch`](phase1_csv_fix.patch) is a reviewable fix against
the PR head: explicitly write CRLF at each CSV call site and correct the
`.gitattributes` comment. This preserves the committed tables and avoids
mixing a result reserialization into Phase 2. `git apply --check` validates
that the patch applies to the reviewed `5c06db2` base. The local Phase 1
commit applies this writer/comment fix and adds
`.github/workflows/phase1-reproducibility.yml`. Its diff is exactly those
three files: no Phase 2 code, research, config or experiment files.

The [actual-writer validation record](../../data/phase2/review/followup/phase1_validation.json)
records available tests and temporary `cli analyse` / `cli report` checks.
No existing input, result or report file was overwritten. Fresh Windows
checkouts can have a CRLF report; the renderer emits the committed LF form,
which also matches the original workspace's report bytes.

CI installs only the pinned Phase 1 requirements, runs the offline test
suite, regenerates numeric results/report, and fails on a tracked diff in
inputs, summary, tables or report. The varying run timestamp/platform is
excluded from that result-drift check. On a later Phase 2 branch, the same
workflow can also run the added offline tests without installing Boltz.
After an approved push/dispatch, require a green native Linux run before
merging. CI has not yet established that result.

A small PR-description edit is also recommended: replace the shorthand
“25V, compound 5” with “deposited N-oxide 25V; the paper attributes 4MEO
to compound 5, and sample identity remains unresolved.” Do not add the
unreproduced map/PDB-REDO results as new baseline findings.

This review neither merges the PR nor posts a message to GitHub. Phase 1's
fix remains a separate commit and worktree. Phase 2 work is on
`phase-2-controlled-pilot`; its review diff is taken against the corrected
local `phase-1-baseline`. The remote PR and Phase 1 results remain intact.
