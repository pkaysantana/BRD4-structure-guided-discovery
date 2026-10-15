# Pilot configurations

- `4men_primary.json`: primary Arm A timing/structural evaluation,
  one seed and one sample using a real cached full MSA. Null MSA metadata
  and SHA-256 fields are pending inputs, not valid provenance. `plan`
  can inspect this design offline; `validate` and `prepare` must fail until
  the actual file and provenance are pinned.
- `4men_smoke.json`: optional no-MSA smoke test only. The code records
  that role explicitly in the prepared and analysed artefacts. The current
  ID is `4men_smoke_s20261019_v3`; v2 is preserved as a superseded record
  after the terminal-OXT compatibility fix.

Both use the same verified construct, neutral 25K graph, Boltz version,
model asset pins, seed and sampling settings. Full MSA is the primary
scientific condition. No live MSA refresh, weights download, model launch
or multi-seed command is implemented.

The earlier `4men_pilot.json` was renamed to `4men_smoke.json` and given
a new experiment ID. Its original prepared bundle remains an unchanged
historical record; it is not the current primary protocol.
