"""Read-only environment checks for use on the GPU host AFTER user approval.

Captures provenance; never downloads a file or invokes the Boltz CLI.
"""

from __future__ import annotations

import hashlib
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, distributions, version
from pathlib import Path

from brd4_analysis.phase2.config import read_json, sha256, write_json
from brd4_analysis.phase2.inputs import atom_names, boltz_ligand
from brd4_analysis.phase2.workflow import model_provenance, planned_command, verify_run


def capture_runtime(run: Path) -> Path:
    run = Path(run).resolve()
    config = verify_run(run)
    try:
        installed_boltz = version("boltz")
    except PackageNotFoundError as exc:
        raise ValueError("Boltz is not installed; this command belongs in the separately approved GPU environment") from exc
    if installed_boltz != config["boltz"]["version"] or version("rdkit") != config["boltz"]["rdkit_version"]:
        raise ValueError("Unreviewed Boltz/RDKit runtime; stop before prediction")
    if platform.system() != "Linux" or sys.version_info[:2] != (3, 12):
        raise ValueError("Pilot specification requires a separate Linux / Python 3.12 environment")
    names = atom_names(boltz_ligand(config["ligand"]["smiles"]))
    if names != read_json(run / "ligand_atom_names.json")["smiles_input_order"]:
        raise ValueError("Runtime ligand atom naming changed")
    cache = Path(config["boltz"]["cache_dir"])
    weights = []
    for item in model_provenance()["files"]:
        path = cache / item["filename"]
        if not path.is_file() or path.stat().st_size != item["size_bytes"] or sha256(path) != item["sha256"]:
            raise ValueError(f"Missing or changed model asset: {path.name}; no download was attempted")
        weights.append(item)
    mols = cache / "mols"
    if not mols.is_dir():
        raise ValueError("Verified mols.tar must be extracted in the approved cache before inference")
    tree_digest = hashlib.sha256()
    count = 0
    for path in sorted(mols.rglob("*")):
        if path.is_symlink():
            raise ValueError("Unexpected link in the CCD cache")
        if path.is_file():
            tree_digest.update(f"{path.relative_to(mols).as_posix()}\t{sha256(path)}\n".encode())
            count += 1
    if not count:
        raise ValueError("Empty CCD cache")
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
    write_json(output, {
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "prepared_manifest_sha256": sha256(run / "manifest.json"),
        "python": sys.version, "platform": platform.platform(),
        "packages": dict(sorted((d.metadata["Name"], d.version) for d in distributions() if d.metadata.get("Name"))),
        "torch_cuda_build": torch.version.cuda,
        "cudnn_version": torch.backends.cudnn.version(),
        "gpu_nvidia_smi_csv": gpu, "free_disk_bytes": free_disk,
        "reproducibility_environment": {name: os.environ.get(name) for name in [
            "CUDA_VISIBLE_DEVICES", "CUBLAS_WORKSPACE_CONFIG", "CUEQ_DEFAULT_CONFIG",
            "CUEQ_DISABLE_AOT_TUNING", "PYTHONHASHSEED",
        ]},
        "model_assets": weights, "model_repository_revision": model_provenance()["revision"],
        "ccd_tree_sha256": tree_digest.hexdigest(), "ccd_file_count": count,
        "msa": config["msa"], "command": planned_command(config),
        "boltz_executed_by_this_command": False,
        "user_approval": "must be recorded separately; this manifest does not grant approval",
    })
    return output
