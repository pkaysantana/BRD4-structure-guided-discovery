"""Runtime-gate regressions with tiny synthetic assets and mocked GPU services.

No real model files, hardware probes, Torch imports or inference are needed.
"""

from copy import deepcopy
from datetime import datetime, timedelta
import io
import json
import os
from pathlib import Path
import tarfile
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.phase2 import config as config_module, runtime, workflow
from brd4_analysis.phase2.config import read_json, sha256, source_path

CONFIG = PROJECT_ROOT / "configs/phase2/4men_smoke.json"
PACKAGES = {
    "boltz": "2.2.1", "rdkit": "2026.3.6", "torch": "2.7.1+cu128",
    "numpy": "1.26.4", "biopython": "1.84", "pytorch-lightning": "2.5.0",
    "scipy": "1.13.1", "requests": "2.32.3", "pyyaml": "6.0.2",
}


def archive_file(path, entries):
    with tarfile.open(path, "w") as archive:
        for entry in entries:
            if isinstance(entry, tarfile.TarInfo):
                archive.addfile(entry)
            else:
                name, data = entry
                info = tarfile.TarInfo(name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))


@pytest.fixture
def cache_tree(tmp_path):
    cache = tmp_path / "cache"
    mols = cache / "mols"
    (mols / "nested").mkdir(parents=True)
    contents = {"ALA.pkl": b"synthetic ALA definition", "nested/GLY.pkl": b"synthetic GLY definition"}
    for relative, data in contents.items():
        (mols / relative).write_bytes(data)
    archive = cache / "mols.tar"
    archive_file(archive, [("mols/" + name, data) for name, data in contents.items()])
    return cache, archive, mols


def test_verified_archive_and_exact_extracted_tree_pass(cache_tree):
    _, archive, mols = cache_tree
    result = runtime.verify_ccd_cache(archive, mols, sha256(archive))
    assert result["archive_sha256"] == sha256(archive)
    assert result["file_count"] == 2
    assert result["verification"] == "exact_archive_inventory_and_contents"
    assert len(result["tree_sha256"]) == 64


@pytest.mark.parametrize("change", ["missing", "same_size_tampering", "unexpected_file", "unexpected_directory"])
def test_incomplete_or_changed_extraction_is_rejected(cache_tree, change):
    _, archive, mols = cache_tree
    if change == "missing":
        (mols / "ALA.pkl").unlink()
    elif change == "same_size_tampering":
        path = mols / "ALA.pkl"
        path.write_bytes(path.read_bytes()[::-1])
    elif change == "unexpected_file":
        (mols / "unexpected.pkl").write_bytes(b"not in the archive")
    else:
        (mols / "unexpected").mkdir()
    with pytest.raises(ValueError, match="CCD"):
        runtime.verify_ccd_cache(archive, mols, sha256(archive))


def test_unpinned_archive_is_rejected(cache_tree):
    _, archive, mols = cache_tree
    with pytest.raises(ValueError, match="checksum"):
        runtime.verify_ccd_cache(archive, mols, "0" * 64)


def test_archive_file_hardlinks_are_rejected(cache_tree, tmp_path):
    _, archive, mols = cache_tree
    linked = tmp_path / "linked.tar"
    os.link(archive, linked)
    with pytest.raises(ValueError, match="non-link"):
        runtime.verify_ccd_cache(archive, mols, sha256(archive))


@pytest.mark.parametrize("name", [
    "../escape.pkl", "/mols/ALA.pkl", "mols/../escape.pkl",
    "mols\\ALA.pkl", "C:/mols/ALA.pkl", "mols/ALA.pkl:stream", "mols/ALA\n.pkl",
])
def test_unsafe_archive_paths_are_rejected(cache_tree, name):
    _, archive, mols = cache_tree
    archive_file(archive, [(name, b"synthetic")])
    with pytest.raises(ValueError, match="Unsafe"):
        runtime.verify_ccd_cache(archive, mols, sha256(archive))


@pytest.mark.parametrize("kind", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE, tarfile.CHRTYPE])
def test_archive_links_and_special_entries_are_rejected(cache_tree, kind):
    _, archive, mols = cache_tree
    member = tarfile.TarInfo("mols/not_a_regular_file")
    member.type = kind
    member.linkname = "ALA.pkl" if kind in {tarfile.SYMTYPE, tarfile.LNKTYPE} else ""
    archive_file(archive, [member])
    with pytest.raises(ValueError, match="Links and special"):
        runtime.verify_ccd_cache(archive, mols, sha256(archive))


@pytest.mark.parametrize("alias", ["mols/ALA.pkl", "./mols/ALA.pkl", "mols/./ALA.pkl"])
def test_duplicate_archive_entries_are_rejected_after_normalization(cache_tree, alias):
    _, archive, mols = cache_tree
    archive_file(archive, [("mols/ALA.pkl", b"first"), (alias, b"second")])
    with pytest.raises(ValueError, match="Duplicate"):
        runtime.verify_ccd_cache(archive, mols, sha256(archive))


def test_extracted_hardlinks_are_rejected(cache_tree, tmp_path):
    _, archive, mols = cache_tree
    path = mols / "ALA.pkl"
    other = tmp_path / "external.pkl"
    other.write_bytes(path.read_bytes())
    path.unlink()
    os.link(other, path)
    with pytest.raises(ValueError, match="linked"):
        runtime.verify_ccd_cache(archive, mols, sha256(archive))


def test_extracted_symlinks_are_rejected(cache_tree, tmp_path):
    _, archive, mols = cache_tree
    path = mols / "ALA.pkl"
    other = tmp_path / "external.pkl"
    other.write_bytes(path.read_bytes())
    path.unlink()
    try:
        path.symlink_to(other)
    except OSError as exc:
        pytest.skip(f"Symlink creation is unavailable on this host: {exc}")
    with pytest.raises(ValueError, match="Link"):
        runtime.verify_ccd_cache(archive, mols, sha256(archive))


@pytest.fixture
def prepared_runtime(tmp_path, cache_tree, monkeypatch, request):
    """Use real preparation/validation, replacing only assets and GPU services."""
    cache, archive, mols = cache_tree
    for filename in ["boltz2_conf.ckpt", "boltz2_aff.ckpt"]:
        (cache / filename).write_bytes(b"synthetic fixture only; not model weights")
    assets = [
        {"filename": name, "size_bytes": (cache / name).stat().st_size, "sha256": sha256(cache / name),
         "url": "https://example.invalid/synthetic/" + name}
        for name in ["mols.tar", "boltz2_conf.ckpt", "boltz2_aff.ckpt"]
    ]
    model = {
        "repository": "synthetic-fixture-only", "revision": "a" * 40,
        "files": assets, "downloaded": False, "note": "synthetic fixture only",
    }
    monkeypatch.setattr(workflow, "model_provenance", lambda: deepcopy(model))
    configuration = read_json(CONFIG)
    if getattr(request, "param", "smoke") == "local_a3m":
        relative = "data/phase2/msa/runtime_synthetic.a3m"
        msa = tmp_path / relative
        msa.parent.mkdir(parents=True)
        query = configuration["protein"]["sequence"]
        msa.write_bytes(f">query\n{query}\n>synthetic_test_only\n{query[:2]}Q{query[3:]}\n".encode())
        configuration["msa"] = {
            "mode": "local_a3m",
            "file": {"path": relative, "sha256": sha256(msa)},
            "provenance": {
                "generator": "synthetic unit test", "database": "not a biological MSA",
                "database_release": "not applicable", "retrieved_utc": "2026-10-15T00:00:00Z",
                "server": "none",
            },
        }

        def lookup(spec, root=PROJECT_ROOT, **kwargs):
            return source_path(spec, tmp_path if spec["path"] == relative else root, **kwargs)

        monkeypatch.setattr(config_module, "source_path", lookup)
        monkeypatch.setattr(workflow, "source_path", lookup)
    config_path = tmp_path / "fixture_config.json"
    config_path.write_text(json.dumps(configuration), encoding="utf-8")
    run = workflow.prepare(config_path, tmp_path / "runs")
    config = read_json(run / "config.json")
    actual_path = Path

    def mapped_path(value):
        return cache if value == config["boltz"]["cache_dir"] else actual_path(value)

    monkeypatch.setattr(runtime, "Path", mapped_path)
    monkeypatch.setattr(runtime, "runtime_packages", lambda: dict(PACKAGES))
    monkeypatch.setattr(runtime.platform, "system", lambda: "Linux")
    monkeypatch.setattr(runtime.platform, "platform", lambda: "Linux-synthetic-fixture")
    monkeypatch.setattr(runtime, "sys", SimpleNamespace(
        version_info=(3, 12, 10), version="3.12.10 (synthetic fixture)",
    ))
    gpu = Mock(return_value=SimpleNamespace(stdout="Synthetic GPU, GPU-fixture, 600.00, 12288\n"))
    monkeypatch.setattr(runtime.subprocess, "run", gpu)
    monkeypatch.setattr(runtime.shutil, "disk_usage", lambda _: SimpleNamespace(free=64 * 10**9))
    # Deliberately no real Torch import or CUDA call.
    monkeypatch.setitem(__import__("sys").modules, "torch", SimpleNamespace(
        __version__=PACKAGES["torch"], version=SimpleNamespace(cuda="12.8"),
        cuda=SimpleNamespace(is_available=lambda: True, device_count=lambda: 1),
        backends=SimpleNamespace(cudnn=SimpleNamespace(version=lambda: 90701)),
    ))
    return run, config, cache, archive, mols, gpu


@pytest.mark.parametrize("change", ["missing", "tampered", "unexpected"])
def test_bad_ccd_cache_stops_capture_before_gpu_probes(prepared_runtime, change):
    run, _, _, _, mols, gpu = prepared_runtime
    if change == "missing":
        (mols / "ALA.pkl").unlink()
    elif change == "tampered":
        path = mols / "ALA.pkl"
        path.write_bytes(path.read_bytes()[::-1])
    else:
        (mols / "unexpected.pkl").write_bytes(b"unexpected")
    with pytest.raises(ValueError, match="CCD"):
        runtime.capture_runtime(run)
    gpu.assert_not_called()
    assert not (run / "runtime.json").exists()
    assert not (run / "runtime.sha256").exists()


def test_capture_emits_a_complete_validated_record_and_checksum(prepared_runtime):
    run, config, _, archive, _, gpu = prepared_runtime
    runtime.capture_runtime(run)
    record = runtime.load_runtime_manifest(run, config)
    assert gpu.call_count == 1  # Mock only; no actual subprocess ran.
    assert record["ccd_cache"]["archive_sha256"] == sha256(archive)
    assert record["ccd_cache"]["file_count"] == 2
    assert record["schema_version"] == runtime.RUNTIME_SCHEMA_VERSION
    assert record["boltz_executed_by_this_command"] is False
    gpu.reset_mock()
    with pytest.raises(ValueError, match="already exists"):
        runtime.capture_runtime(run)
    gpu.assert_not_called()


def rewrite_runtime(run, record, *, resign=True):
    (run / "runtime.json").write_bytes((json.dumps(record, indent=2) + "\n").encode())
    if resign:
        (run / "runtime.sha256").write_text(sha256(run / "runtime.json") + "\n", encoding="ascii")


@pytest.fixture
def scoreable_run(prepared_runtime, monkeypatch):
    run, config, _, _, _, _ = prepared_runtime
    runtime.capture_runtime(run)
    folder = run / "outputs/boltz_results_4men/predictions/4men"
    folder.mkdir(parents=True)
    (folder / "4men_model_0.cif").write_bytes(b"synthetic placeholder; scoring is mocked")
    (folder / "confidence_4men_model_0.json").write_text(json.dumps({
        key: 0.9 for key in ["confidence_score", "ptm", "iptm", "ligand_iptm", "complex_plddt"]
    }), encoding="utf-8")
    score = Mock(return_value={"synthetic_score_only": True})
    monkeypatch.setattr(workflow, "score_prediction", score)
    return run, config, read_json(run / "runtime.json"), score


def test_complete_runtime_record_allows_scoring(scoreable_run):
    run, _, _, score = scoreable_run
    output = workflow.analyse_run(run)
    score.assert_called_once()
    assert read_json(output)["synthetic_score_only"] is True


def test_single_hash_runtime_record_is_rejected_before_scoring(scoreable_run):
    run, _, record, score = scoreable_run
    rewrite_runtime(run, {"prepared_manifest_sha256": record["prepared_manifest_sha256"]})
    with pytest.raises(ValueError, match="runtime manifest"):
        workflow.analyse_run(run)
    score.assert_not_called()
    assert not (run / "analysis").exists()


@pytest.mark.parametrize("field,value", [
    ("schema_version", 0), ("experiment_id", "different-experiment"),
    ("prepared_manifest_sha256", "0" * 64), ("input_sha256", "0" * 64),
    ("model_repository_revision", "wrong-revision"),
    ("command", ["boltz", "predict", "other.yaml"]),
    ("msa", {"mode": "local_a3m"}), ("msa_input_sha256", "0" * 64),
    ("python", "3.13.0"), ("platform", "Windows"),
    ("torch_version", "2.8.0+cu128"), ("torch_cuda_build", "12.6"),
    ("gpu_nvidia_smi_csv", ""), ("free_disk_bytes", 1),
    ("reproducibility_environment", {}),
    ("boltz_executed_by_this_command", True), ("user_approval", "approved"),
])
def test_inconsistent_runtime_fields_are_rejected_even_with_a_matching_checksum(scoreable_run, field, value):
    run, _, record, score = scoreable_run
    record[field] = value
    rewrite_runtime(run, record)
    with pytest.raises(ValueError):
        workflow.analyse_run(run)
    score.assert_not_called()
    assert not (run / "analysis").exists()


@pytest.mark.parametrize("change", ["missing_field", "extra_field", "package_version", "package_missing", "package_name",
                                     "asset_hash", "asset_missing", "asset_duplicate", "ccd_archive", "ccd_method",
                                     "naive_timestamp", "old_timestamp"])
def test_malformed_or_inconsistent_provenance_cannot_reach_scoring(scoreable_run, change):
    run, _, record, score = scoreable_run
    if change == "missing_field":
        del record["gpu_nvidia_smi_csv"]
    elif change == "extra_field":
        record["unreviewed"] = True
    elif change == "package_version":
        record["packages"]["numpy"] = "2.5.3"
    elif change == "package_missing":
        del record["packages"]["pytorch-lightning"]
    elif change == "package_name":
        record["packages"]["invalid package name"] = "1.0"
    elif change == "asset_hash":
        record["model_assets"][0]["sha256"] = "0" * 64
    elif change == "asset_missing":
        record["model_assets"].pop()
    elif change == "asset_duplicate":
        record["model_assets"].append(deepcopy(record["model_assets"][0]))
    elif change == "ccd_archive":
        record["ccd_cache"]["archive_sha256"] = "0" * 64
    elif change == "ccd_method":
        record["ccd_cache"]["verification"] = "tree_hash_only"
    elif change == "naive_timestamp":
        record["recorded_utc"] = "2026-10-15T00:00:00"
    else:
        prepared = datetime.fromisoformat(read_json(run / "manifest.json")["prepared_utc"])
        record["recorded_utc"] = (prepared - timedelta(seconds=1)).isoformat()
    rewrite_runtime(run, record)
    with pytest.raises(ValueError):
        workflow.analyse_run(run)
    score.assert_not_called()


def test_altered_runtime_manifest_is_rejected_before_scoring(scoreable_run):
    run, _, record, score = scoreable_run
    record["gpu_nvidia_smi_csv"] = "Different GPU, GPU-other, 600.00, 24576"
    rewrite_runtime(run, record, resign=False)
    with pytest.raises(ValueError, match="checksum"):
        workflow.analyse_run(run)
    score.assert_not_called()


def test_missing_runtime_receipt_is_rejected_before_scoring(scoreable_run):
    run, _, _, score = scoreable_run
    (run / "runtime.sha256").unlink()
    with pytest.raises(OSError):
        workflow.analyse_run(run)
    score.assert_not_called()


@pytest.mark.parametrize("prepared_runtime", ["local_a3m"], indirect=True)
def test_primary_runtime_binds_the_actual_frozen_msa(scoreable_run):
    run, config, record, score = scoreable_run
    assert record["msa_input_sha256"] == sha256(run / "inputs/protein.a3m")
    assert record["msa_input_sha256"] == config["msa"]["file"]["sha256"]
    workflow.analyse_run(run)
    score.assert_called_once()


@pytest.mark.parametrize("prepared_runtime", ["local_a3m"], indirect=True)
@pytest.mark.parametrize("change", ["runtime_hash", "runtime_provenance", "frozen_msa"])
def test_primary_msa_inconsistency_blocks_scoring(scoreable_run, change):
    run, _, record, score = scoreable_run
    if change == "runtime_hash":
        record["msa_input_sha256"] = "0" * 64
    elif change == "runtime_provenance":
        record["msa"]["provenance"]["generator"] = "different generator"
    else:
        # Even resealing the prepared file list must not break the binding to
        # the checksum in the frozen configuration and runtime record.
        path = run / "inputs/protein.a3m"
        path.write_bytes(path.read_bytes() + b"\n")
        manifest = read_json(run / "manifest.json")
        manifest["files"]["inputs/protein.a3m"] = sha256(path)
        (run / "manifest.json").write_bytes(json.dumps(manifest).encode())
        record["prepared_manifest_sha256"] = sha256(run / "manifest.json")
    rewrite_runtime(run, record)
    with pytest.raises(ValueError, match="MSA"):
        workflow.analyse_run(run)
    score.assert_not_called()
