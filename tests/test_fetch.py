"""fetch must preserve original inputs: no network access or rewrite when
checksums match, and refusal (not silent overwrite) on mismatch."""

import shutil

import pytest

import brd4_analysis.fetch as fetch


@pytest.fixture
def data_copy(tmp_path, data_dir):
    dest = tmp_path / "data"
    shutil.copytree(data_dir, dest)
    return dest


def _no_network(*args, **kwargs):
    raise AssertionError("fetch attempted a download although inputs verified")


def test_verified_inputs_are_kept_without_download(data_copy, monkeypatch):
    monkeypatch.setattr(fetch.requests, "get", _no_network)
    manifest = data_copy / "provenance" / "provenance.json"
    before = manifest.read_bytes()
    fetch.fetch_all(data_copy)
    assert manifest.read_bytes() == before


def test_checksum_mismatch_refuses_to_overwrite(data_copy, monkeypatch):
    monkeypatch.setattr(fetch.requests, "get", _no_network)
    pdb = data_copy / "raw" / "4MEN.pdb"
    pdb.write_bytes(pdb.read_bytes() + b"REMARK tampered\n")
    with pytest.raises(RuntimeError, match="does not match its recorded SHA-256"):
        fetch.fetch_all(data_copy)
    assert pdb.read_bytes().endswith(b"REMARK tampered\n")
