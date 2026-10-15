# Phase 1 merge recommendation

**Do not merge PR #1 at its current head solely on the strength of the
supplied research.** The scientific framing is acceptable as a baseline of
deposited models, but a concrete CSV portability issue remains. A narrow
fix can preserve every existing result byte.

Reviewed PR:
<https://github.com/pkaysantana/BRD4-structure-guided-discovery/pull/1>,
head `5c06db21811a29e9925ab008ca555129c71c0a35`, checked 15 October 2026.
GitHub reports OPEN / MERGEABLE / CLEAN, with no status checks listed.
Mergeability is not a validation result.

## Evidence and remaining gate

| Review item | Finding |
|---|---|
| Earlier review: template-assigned protonation | Addressed: the report and code state graph consistency and its limits |
| Earlier review: stale ECFP4 similarity | Addressed: heavy-atom Morgan value is 6/72; regression and report consistency tests cover it |
| 25V atom roles | Corrected and tested at this head |
| 4MEO identity | Properly unresolved in the baseline; new research does not settle it |
| Native Windows baseline | All 45 existing tests pass, including exact numeric outputs, report rendering and determinism |
| Review on the current head: CSV line endings | Still open in code: five `to_csv` call sites use the OS default; all eight committed CSVs are CRLF |
| Native Linux validation | Not performed in this Windows workspace |

With pandas' default LF on Linux, all eight CSVs differ from the committed
bytes even with unchanged numerical results. `.gitattributes` correctly
preserves the files, but cannot fix a writer that varies by operating
system. The existing determinism test therefore exposes a portability
failure on Linux.

[`phase1_csv_fix.patch`](phase1_csv_fix.patch) is a reviewable fix against
the PR head: explicitly write CRLF at each CSV call site and correct the
`.gitattributes` comment. This preserves the committed tables and avoids
mixing a result reserialization into Phase 2. `git apply --check` validates
that the patch applies. It has **not** been applied to PR #1 or to this
workspace's Phase 1 source.

The separate validation record tests LF-default behaviour and the
explicit-CRLF policy using temporary outputs; it is an OS-default
simulation, not a claim of a native Linux run. After applying the fix to
the Phase 1 branch, run the existing tests on Windows and Linux. Regenerated
summary, tables and report should retain their current contents.
That closes the technical merge gate; retain the uncertainty language.

A small PR-description edit is also recommended: replace the shorthand
“25V, compound 5” with “deposited N-oxide 25V; the paper attributes 4MEO
to compound 5, and sample identity remains unresolved.” Do not add the
unreproduced map/PDB-REDO results as new baseline findings.

This review neither merges the PR nor posts a message to GitHub. Phase 2
work is isolated on `phase-2-controlled-pilot`; the PR head and Phase 1
results remain intact.
