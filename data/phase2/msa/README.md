# Primary MSA input gate

`4men_query.fasta` is the verified 127-residue input query, including the
SM expression remnant. It is **not an MSA**. No real alignment has been
acquired, and no MSA service has been called.

After approval, obtain the full alignment for that exact query from the
chosen source (proposed: MMseqs2 through the ColabFold service). Preserve
the original response as `4men_full.a3m`, record its SHA-256, generator,
database and available database-release identifier, retrieval time and
service URL in `configs/phase2/4men_primary.json`. Record unavailable
database-version information as `unknown`, not a guessed snapshot.

Review the source and query mapping before preparing the experiment.
A file containing only the query is rejected for the primary condition.
Synthetic test alignments live only in pytest temporary directories and
must never be copied here as scientific input.

No command in this repository performs the acquisition. Do not run
`boltz predict` merely to fetch an MSA: the upstream command can first
download model assets and subsequently start inference.
