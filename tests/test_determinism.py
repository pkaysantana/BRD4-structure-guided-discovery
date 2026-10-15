"""Reproducibility: re-running the analysis on the committed inputs must
reproduce the committed numerical outputs byte-for-byte, and two runs must
produce identical figures."""

from pathlib import Path

import pytest

from brd4_analysis.pipeline import run_analysis

PROJECT_ROOT = Path(__file__).resolve().parents[1]
COMMITTED = PROJECT_ROOT / "results"
NUMERIC_OUTPUTS = sorted(
    p.relative_to(COMMITTED).as_posix()
    for p in [COMMITTED / "summary.json", *(COMMITTED / "tables").glob("*.csv")]
)


@pytest.fixture(scope="module")
def two_runs(tmp_path_factory, data_dir):
    runs = []
    for i in range(2):
        out = tmp_path_factory.mktemp(f"run{i}")
        run_analysis(data_dir, out)
        runs.append(out)
    return runs


def test_numeric_output_inventory():
    assert "summary.json" in NUMERIC_OUTPUTS
    assert len(NUMERIC_OUTPUTS) == 9  # summary.json + 8 CSV tables


@pytest.mark.parametrize("rel", NUMERIC_OUTPUTS)
def test_numeric_outputs_match_committed(two_runs, rel):
    assert (two_runs[0] / rel).read_bytes() == (COMMITTED / rel).read_bytes()


def test_runs_identical_including_figures(two_runs):
    a, b = two_runs
    files = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file() and p.name != "run_metadata.json")
    assert files == sorted(p.relative_to(b) for p in b.rglob("*") if p.is_file() and p.name != "run_metadata.json")
    for rel in files:
        assert (a / rel).read_bytes() == (b / rel).read_bytes(), rel
