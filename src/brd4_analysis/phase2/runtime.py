"""Read-only environment checks for use on the GPU host AFTER user approval.

Captures provenance; never downloads a file or invokes the Boltz CLI.
"""

from __future__ import annotations

import hashlib
import csv
import io
import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import tarfile
from datetime import datetime, timezone
from email.parser import Parser
from importlib.metadata import distributions
from pathlib import Path, PurePosixPath

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name
from packaging.version import InvalidVersion, Version

from brd4_analysis.phase2.config import SOURCE_DIR, keys, read_json, sha256, validate_a3m, write_json
from brd4_analysis.phase2.inputs import atom_names, boltz_ligand

RUNTIME_SCHEMA_VERSION = 1
APPROVAL_NOTE = "must be recorded separately; this manifest does not grant approval"
ENVIRONMENT_KEYS = {
    "CUDA_VISIBLE_DEVICES", "CUBLAS_WORKSPACE_CONFIG", "CUEQ_DEFAULT_CONFIG",
    "CUEQ_DISABLE_AOT_TUNING", "PYTHONHASHSEED",
}
RUNTIME_KEYS = {
    "schema_version", "experiment_id", "recorded_utc", "prepared_manifest_sha256",
    "input_sha256", "msa_input_sha256", "python", "platform", "packages",
    "torch_version", "torch_cuda_build", "cudnn_version", "gpu_nvidia_smi_csv", "free_disk_bytes",
    "reproducibility_environment", "model_assets", "model_repository_revision",
    "ccd_cache", "msa", "command", "boltz_executed_by_this_command", "user_approval",
}
REQUIRED_METADATA_PACKAGES = {
    "torch", "numpy", "rdkit", "biopython", "pytorch-lightning",
    "scipy", "requests", "pyyaml",
}


def _stream_sha256(handle) -> str:
    digest = hashlib.sha256()
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
    return digest.hexdigest()


def _is_link(path: Path) -> bool:
    return path.is_symlink() or getattr(path, "is_junction", lambda: False)()


def verify_ccd_cache(archive_path: Path, mols_dir: Path, expected_sha256: str) -> dict:
    """Compare a link-free extracted tree with the pinned tar, without extracting.

    The upstream archive is rooted at mols/. Names are normalized as POSIX
    paths; traversal, absolute/drive paths, links, special files and duplicate
    entries are refused before any GPU operation can be reached.
    """
    archive_path, mols_dir = Path(archive_path), Path(mols_dir)
    if _is_link(archive_path) or not archive_path.is_file() or archive_path.stat().st_nlink != 1:
        raise ValueError("CCD archive must be a regular, non-link file")
    if _is_link(mols_dir) or not mols_dir.is_dir():
        raise ValueError("CCD cache must be an extracted, non-link directory")
    with archive_path.open("rb") as raw:
        if _stream_sha256(raw) != expected_sha256:
            raise ValueError("CCD archive checksum mismatch")
        raw.seek(0)
        with tarfile.open(fileobj=raw, mode="r:*") as archive:
            expected_files, expected_dirs, seen = {}, set(), set()
            for member in archive:
                name = member.name
                path = PurePosixPath(name)
                if (
                    not name or path.is_absolute() or ".." in path.parts
                    or "\\" in name or ":" in name
                    or any(ord(c) < 32 or ord(c) == 127 for c in name)
                ):
                    raise ValueError(f"Unsafe CCD archive path: {name!r}")
                normalized = path.as_posix()
                if normalized in seen:
                    raise ValueError(f"Duplicate CCD archive entry: {normalized}")
                seen.add(normalized)
                if not member.isdir() and member.type not in {tarfile.REGTYPE, tarfile.AREGTYPE}:
                    raise ValueError(f"Links and special CCD archive entries are forbidden: {normalized}")
                if normalized == "." and member.isdir():
                    continue
                if not path.parts or path.parts[0] != "mols":
                    raise ValueError(f"CCD archive entry is outside mols/: {normalized}")
                if len(path.parts) == 1:
                    if not member.isdir():
                        raise ValueError("CCD archive root mols must be a directory")
                    continue
                relative = PurePosixPath(*path.parts[1:])
                for parent in relative.parents:
                    if parent.as_posix() != ".":
                        expected_dirs.add(parent.as_posix())
                if member.isdir():
                    expected_dirs.add(relative.as_posix())
                else:
                    if name.endswith("/"):
                        raise ValueError("Ambiguous CCD archive filename")
                    expected_files[relative.as_posix()] = member
            if not expected_files or set(expected_files) & expected_dirs:
                raise ValueError("Empty or conflicting CCD archive inventory")

            actual_files, actual_dirs = {}, set()
            for path in mols_dir.rglob("*"):
                if _is_link(path):
                    raise ValueError(f"Link in extracted CCD cache: {path.name}")
                info = path.lstat()
                relative = path.relative_to(mols_dir).as_posix()
                if stat.S_ISDIR(info.st_mode):
                    actual_dirs.add(relative)
                elif stat.S_ISREG(info.st_mode) and info.st_nlink == 1:
                    actual_files[relative] = path
                else:
                    raise ValueError(f"Non-regular or linked CCD file: {relative}")
            if set(actual_files) != set(expected_files) or actual_dirs != expected_dirs:
                raise ValueError("Extracted CCD inventory mismatch (missing or unexpected files/directories)")
            file_hashes = {}
            for relative, member in expected_files.items():
                path = actual_files[relative]
                if path.stat().st_size != member.size:
                    raise ValueError(f"Extracted CCD file size mismatch: {relative}")
                handle = archive.extractfile(member)
                if handle is None:
                    raise ValueError(f"Unreadable CCD archive member: {relative}")
                with handle:
                    expected_hash = _stream_sha256(handle)
                if sha256(path) != expected_hash:
                    raise ValueError(f"Extracted CCD file content mismatch: {relative}")
                file_hashes[relative] = expected_hash
    tree_digest = hashlib.sha256()
    for relative, digest in sorted(file_hashes.items()):
        tree_digest.update(f"{relative}\t{digest}\n".encode())
    return {
        "archive_sha256": expected_sha256,
        "tree_sha256": tree_digest.hexdigest(),
        "file_count": len(file_hashes),
        "verification": "exact_archive_inventory_and_contents",
    }


def runtime_packages() -> dict:
    packages = {}
    for distribution in distributions():
        name = distribution.metadata.get("Name")
        if name:
            name = canonicalize_name(name, validate=True)
            if name in packages and packages[name] != distribution.version:
                raise ValueError(f"Conflicting installed package versions: {name}")
            packages[name] = distribution.version
    return dict(sorted(packages.items()))


def validate_runtime_packages(packages: dict, config: dict) -> None:
    if not isinstance(packages, dict) or not packages:
        raise ValueError("Runtime packages must be recorded")
    for name, value in packages.items():
        if not isinstance(name, str) or canonicalize_name(name, validate=True) != name or not isinstance(value, str):
            raise ValueError("Runtime package names/versions must be canonical strings")
        try:
            Version(value)
        except InvalidVersion as exc:
            raise ValueError(f"Invalid runtime package version: {name}") from exc
    for name, expected in [("boltz", config["boltz"]["version"]), ("rdkit", config["boltz"]["rdkit_version"])]:
        if packages.get(name) != expected:
            raise ValueError(f"Runtime {name} version differs from the frozen configuration")
    # Reuse the already pinned wheel metadata; never consult a live package index.
    metadata_text = (SOURCE_DIR / "wheel_boltz-2.2.1.dist-info_METADATA.txt").read_text(encoding="utf-8")
    checked = set()
    for value in Parser().parsestr(metadata_text, headersonly=True).get_all("Requires-Dist", []):
        requirement = Requirement(value)
        name = canonicalize_name(requirement.name)
        if requirement.marker is None:
            if name not in packages or not requirement.specifier.contains(packages[name], prereleases=False):
                raise ValueError(f"Runtime {name} does not satisfy the pinned Boltz requirements")
            checked.add(name)
    if not REQUIRED_METADATA_PACKAGES <= checked:
        raise ValueError("Pinned Boltz runtime requirements are incomplete")


def _digest(value, label: str) -> None:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{label} must be a SHA-256 digest")


def _utc_timestamp(value: str):
    if not isinstance(value, str):
        raise ValueError("Runtime timestamp must be a string")
    try:
        stamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError("Invalid runtime timestamp") from exc
    if stamp.tzinfo is None or stamp.utcoffset() is None:
        raise ValueError("Runtime timestamp must include a timezone")
    return stamp


def validate_runtime_manifest(run: Path, config: dict, record: dict) -> None:
    """Check complete provenance against frozen inputs, without GPU/model access."""
    from brd4_analysis.phase2.workflow import input_document, model_provenance, planned_command

    run = Path(run)
    keys(record, RUNTIME_KEYS, "runtime manifest")
    if type(record["schema_version"]) is not int or record["schema_version"] != RUNTIME_SCHEMA_VERSION:
        raise ValueError("Unsupported runtime manifest schema")
    if record["experiment_id"] != config["experiment_id"]:
        raise ValueError("Runtime experiment identity mismatch")
    if record["prepared_manifest_sha256"] != sha256(run / "manifest.json"):
        raise ValueError("Runtime provenance belongs to a different prepared experiment")
    frozen = read_json(run / "manifest.json")
    if _utc_timestamp(record["recorded_utc"]) < _utc_timestamp(frozen["prepared_utc"]):
        raise ValueError("Runtime timestamp predates preparation")
    if record["input_sha256"] != sha256(run / "inputs/4men.yaml"):
        raise ValueError("Runtime model input checksum mismatch")
    if read_json(run / "inputs/4men.yaml") != input_document(config):
        raise ValueError("Frozen model input disagrees with the configuration")
    command = read_json(run / "command.json")["argv"]
    if record["command"] != command or command != planned_command(config):
        raise ValueError("Runtime command differs from the frozen plan")
    if record["msa"] != config["msa"]:
        raise ValueError("Runtime MSA provenance differs from the frozen configuration")
    if config["msa"]["mode"] == "local_a3m":
        msa = run / "inputs/protein.a3m"
        expected = config["msa"]["file"]["sha256"]
        if record["msa_input_sha256"] != expected or sha256(msa) != expected:
            raise ValueError("Runtime MSA checksum mismatch")
        validate_a3m(msa, config["protein"]["sequence"])
    elif record["msa_input_sha256"] is not None:
        raise ValueError("Single-sequence runtime must not claim an MSA checksum")

    expected_model = read_json(run / "preregistration.json")["model"]
    if expected_model != model_provenance():
        raise ValueError("Frozen model metadata disagrees with the pinned sources")
    if record["model_repository_revision"] != expected_model["revision"]:
        raise ValueError("Runtime model repository revision mismatch")
    assets = record["model_assets"]
    if not isinstance(assets, list) or any(not isinstance(item, dict) for item in assets):
        raise ValueError("Runtime model assets must be a complete list")
    for item in assets:
        keys(item, {"filename", "sha256", "size_bytes", "url"}, "runtime model asset")
        if not isinstance(item["filename"], str) or type(item["size_bytes"]) is not int:
            raise ValueError("Invalid runtime model asset metadata")
    by_name = {item["filename"]: item for item in assets}
    if len(by_name) != len(assets) or by_name != {a["filename"]: a for a in expected_model["files"]}:
        raise ValueError("Runtime checkpoint/archive metadata differs from the frozen model")
    ccd = record["ccd_cache"]
    keys(ccd, {"archive_sha256", "tree_sha256", "file_count", "verification"}, "CCD verification")
    if ccd["archive_sha256"] != by_name["mols.tar"]["sha256"]:
        raise ValueError("CCD verification is not bound to the pinned archive")
    _digest(ccd["tree_sha256"], "CCD tree checksum")
    if type(ccd["file_count"]) is not int or ccd["file_count"] < 1:
        raise ValueError("Invalid verified CCD file count")
    if ccd["verification"] != "exact_archive_inventory_and_contents":
        raise ValueError("CCD cache was not verified against the archive")

    validate_runtime_packages(record["packages"], config)
    if record["torch_version"] != record["packages"]["torch"]:
        raise ValueError("Imported Torch version differs from the recorded distribution")
    if not isinstance(record["python"], str) or not re.match(r"^3\.12\.\d+(?:\s|$)", record["python"]):
        raise ValueError("Runtime must record Python 3.12")
    if not isinstance(record["platform"], str) or not record["platform"].startswith("Linux"):
        raise ValueError("Runtime must record the planned Linux environment")
    cuda = record["torch_cuda_build"]
    if not isinstance(cuda, str) or not re.fullmatch(r"\d+\.\d+(?:\.\d+)?", cuda):
        raise ValueError("Runtime CUDA build is missing or invalid")
    local = Version(record["packages"]["torch"]).local or ""
    if local.startswith(("cpu", "rocm")):
        raise ValueError("Torch build is inconsistent with a CUDA runtime")
    cuda_tag = re.fullmatch(r"cu(\d{2})(\d{1,2})", local)
    if cuda_tag and tuple(map(int, cuda.split(".")[:2])) != tuple(map(int, cuda_tag.groups())):
        raise ValueError("Torch wheel tag and CUDA build disagree")
    if record["cudnn_version"] is not None and (
        type(record["cudnn_version"]) is not int or record["cudnn_version"] <= 0
    ):
        raise ValueError("Invalid recorded cuDNN version")
    gpu = record["gpu_nvidia_smi_csv"]
    if not isinstance(gpu, str) or not gpu.strip():
        raise ValueError("Runtime GPU/driver provenance is missing")
    for row in csv.reader(io.StringIO(gpu)):
        if len(row) != 4 or any(not field.strip() for field in row):
            raise ValueError("Incomplete GPU/driver provenance")
        if not re.fullmatch(r"\d+(?:\.\d+)+", row[2].strip()) or not re.fullmatch(r"\d+", row[3].strip()) or int(row[3]) <= 0:
            raise ValueError("Invalid GPU driver or memory record")
    if type(record["free_disk_bytes"]) is not int or record["free_disk_bytes"] < config["resources"]["min_free_disk_GB"] * 10**9:
        raise ValueError("Runtime disk allocation is below the frozen requirement")
    keys(record["reproducibility_environment"], ENVIRONMENT_KEYS, "runtime environment")
    if any(value is not None and not isinstance(value, str) for value in record["reproducibility_environment"].values()):
        raise ValueError("Invalid runtime environment values")
    if record["boltz_executed_by_this_command"] is not False or record["user_approval"] != APPROVAL_NOTE:
        raise ValueError("Runtime provenance is not an execution approval")


def load_runtime_manifest(run: Path, config: dict) -> dict:
    """Reject incomplete/altered records before scientific scoring.

    The companion hash detects changes to the captured record. It is not
    cryptographic attestation against someone deliberately forging both files.
    """
    run = Path(run)
    digest = (run / "runtime.sha256").read_text(encoding="ascii").strip()
    _digest(digest, "Runtime receipt")
    if digest != sha256(run / "runtime.json"):
        raise ValueError("Runtime manifest checksum mismatch")
    record = read_json(run / "runtime.json")
    validate_runtime_manifest(run, config, record)
    return record


def capture_runtime(run: Path) -> Path:
    from brd4_analysis.phase2.workflow import model_provenance, planned_command, verify_run

    run = Path(run).resolve()
    config = verify_run(run)
    if (run / "runtime.json").exists() or (run / "runtime.sha256").exists():
        raise ValueError("Runtime provenance already exists; do not overwrite it")
    packages = runtime_packages()
    validate_runtime_packages(packages, config)
    if platform.system() != "Linux" or sys.version_info[:2] != (3, 12):
        raise ValueError("Pilot specification requires a separate Linux / Python 3.12 environment")
    names = atom_names(boltz_ligand(config["ligand"]["smiles"]))
    if names != read_json(run / "ligand_atom_names.json")["smiles_input_order"]:
        raise ValueError("Runtime ligand atom naming changed")
    cache = Path(config["boltz"]["cache_dir"])
    weights = []
    for item in model_provenance()["files"]:
        path = cache / item["filename"]
        if not path.is_file() or path.stat().st_size != item["size_bytes"]:
            raise ValueError(f"Missing or changed model asset: {path.name}; no download was attempted")
        if item["filename"] != "mols.tar" and sha256(path) != item["sha256"]:
            raise ValueError(f"Changed model asset: {path.name}; no download was attempted")
        weights.append(item)
    archive = next(item for item in weights if item["filename"] == "mols.tar")
    ccd = verify_ccd_cache(cache / "mols.tar", cache / "mols", archive["sha256"])
    gpu = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,uuid,driver_version,memory.total", "--format=csv,noheader,nounits"],
        check=True, capture_output=True, text=True, timeout=15,
    ).stdout.strip()
    import torch  # only in this explicit future preflight, never during offline preparation

    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise ValueError("Expose exactly one CUDA GPU for the pilot")
    free_disk = shutil.disk_usage(run).free
    if free_disk < config["resources"]["min_free_disk_GB"] * 10**9:
        raise ValueError("Insufficient free disk for the planned run")
    output = run / "runtime.json"
    record = {
        "schema_version": RUNTIME_SCHEMA_VERSION,
        "experiment_id": config["experiment_id"],
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "prepared_manifest_sha256": sha256(run / "manifest.json"),
        "input_sha256": sha256(run / "inputs/4men.yaml"),
        "msa_input_sha256": sha256(run / "inputs/protein.a3m") if config["msa"]["mode"] == "local_a3m" else None,
        "python": sys.version, "platform": platform.platform(),
        "packages": packages,
        "torch_version": str(torch.__version__), "torch_cuda_build": torch.version.cuda,
        "cudnn_version": torch.backends.cudnn.version(),
        "gpu_nvidia_smi_csv": gpu, "free_disk_bytes": free_disk,
        "reproducibility_environment": {name: os.environ.get(name) for name in sorted(ENVIRONMENT_KEYS)},
        "model_assets": weights, "model_repository_revision": model_provenance()["revision"],
        "ccd_cache": ccd,
        "msa": config["msa"], "command": planned_command(config),
        "boltz_executed_by_this_command": False,
        "user_approval": APPROVAL_NOTE,
    }
    validate_runtime_manifest(run, config, record)
    write_json(output, record)
    with (run / "runtime.sha256").open("xb") as handle:
        handle.write((sha256(output) + "\n").encode("ascii"))
    return output
