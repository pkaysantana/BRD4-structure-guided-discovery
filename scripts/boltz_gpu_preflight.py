"""Legion-only CUDA preflight. Run manually in the dedicated environment.

No installation, network request, model download, Boltz import or inference.
Importing this file performs no diagnostics. main() is the only entry point
that collects hardware facts and runs the tiny CUDA operations.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

ENV_NAME = "brd4-boltz221-cu128"
TORCH_PIN = "2.7.1+cu128"
CUDA_PIN = "12.8"
BOLTZ_PIN = "2.2.1"
GIB = 1024**3


def version_parts(value: str | None) -> tuple[int, ...]:
    match = re.match(r"^\d+(?:\.\d+)*", value or "")
    return tuple(map(int, match[0].split("."))) if match else ()


def is_target_gpu(name: str) -> bool:
    normalized = re.sub(r"[^a-z0-9]", "", name.lower())
    return "rtx5070ti" in normalized and "laptop" in normalized


def parse_smi(text: str, *, with_capability: bool = True) -> list[dict]:
    """Parse the explicitly requested, unit-free nvidia-smi CSV columns."""
    fields = ["index", "name", "driver_version"]
    if with_capability:
        fields.append("compute_cap")
    fields += ["memory_total_mib", "memory_used_mib", "memory_free_mib", "power_limit_w"]
    rows = []
    for row in csv.reader(io.StringIO(text)):
        if not row:
            continue
        if len(row) != len(fields):
            raise ValueError("Unexpected nvidia-smi CSV column count")
        item = dict(zip(fields, (value.strip() for value in row)))
        for key in fields:
            if key.endswith(("_mib", "_w")):
                try:
                    value = float(item[key])
                    item[key] = value if math.isfinite(value) else None
                except ValueError:
                    item[key] = None  # WSL/laptop power telemetry may be unavailable.
        rows.append(item)
    return rows


def assess(facts: dict, min_free_gib: float, require_boltz: bool = False) -> dict:
    """Pure decision logic; tests use synthetic facts, never local hardware."""
    gpu = facts.get("gpu", {})
    packages = facts.get("packages", {})
    memory = facts.get("torch_memory_before_probe", {})
    observed_free = [
        snapshot["free_bytes"] for snapshot in [
            memory, facts.get("torch_memory_after_probe", {}),
        ] if snapshot.get("free_bytes") is not None
    ]
    free = min(observed_free) if observed_free else None
    total = memory.get("total_bytes")
    # Toolkit baselines, not a guarantee that an old driver supports this laptop.
    floor = (570, 65) if facts.get("platform") == "Windows" or facts.get("is_wsl") else (570, 26)
    checks = []

    def check(name, passed, group="compatibility"):
        checks.append({"name": name, "passed": bool(passed), "group": group})

    check("dedicated_virtual_environment", facts.get("isolated_env") is True)
    check("python_3_12", facts.get("python_major_minor") == [3, 12])
    check("expected_laptop_gpu", is_target_gpu(gpu.get("name", "")))
    check("cuda_12_8_driver_baseline", version_parts(gpu.get("driver_version")) >= floor)
    check("torch_pin", facts.get("torch_version") == TORCH_PIN)
    check("torch_cuda_build", facts.get("torch_cuda_build") == CUDA_PIN)
    check("numpy_pin", packages.get("numpy") == "1.26.4")
    check("cuda_available", facts.get("cuda_available") is True)
    check("one_visible_cuda_device", facts.get("cuda_device_count") == 1)
    check("compute_capability_sm120", facts.get("compute_capability") == [12, 0])
    check("native_sm120_in_torch_build", "sm_120" in facts.get("torch_arch_list", []))
    check("bf16_supported", facts.get("bf16_supported") is True)
    check("tiny_fp32_and_bf16_operations", facts.get("tensor_test_passed") is True)
    check("expected_vram_class", total is not None and total >= 11 * GIB)
    check("free_vram_policy", free is not None and free >= min_free_gib * GIB, "headroom")
    check("no_collection_error", not facts.get("error"))
    if require_boltz:
        check("boltz_pin", packages.get("boltz") == BOLTZ_PIN, "boltz_environment")
        check("rdkit_atom_mapping_pin", packages.get("rdkit") == "2026.3.6", "boltz_environment")
        check("biopython_pin", packages.get("biopython") == "1.84", "boltz_environment")
        check("current_pilot_linux_requirement", facts.get("platform") == "Linux", "boltz_environment")
    return {
        "status": "PASS" if all(item["passed"] for item in checks) else "FAIL_OR_INCOMPLETE",
        "compatibility_checks_pass": all(c["passed"] for c in checks if c["group"] == "compatibility"),
        "headroom_policy_pass": all(c["passed"] for c in checks if c["group"] == "headroom"),
        "minimum_observed_free_bytes": free,
        "min_free_gib_policy": min_free_gib,
        "boltz_environment_required": require_boltz,
        "checks": checks,
        "failed_checks": [c["name"] for c in checks if not c["passed"]],
        "scope": "CUDA compatibility and tiny tensor operations only; not a Boltz memory-fit guarantee",
        "model_weights_loaded": False,
        "boltz_inference_executed": False,
    }


def find_smi() -> str:
    candidates = [
        shutil.which("nvidia-smi"),
        str(Path(os.environ.get("SystemRoot", "C:/Windows")) / "System32/nvidia-smi.exe"),
        "C:/Program Files/NVIDIA Corporation/NVSMI/nvidia-smi.exe",
        "/usr/lib/wsl/lib/nvidia-smi",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return candidate
    raise RuntimeError("nvidia-smi was not found; no CUDA operation was attempted")


def collect_facts() -> dict:
    """Called ONLY by a deliberate invocation on the powered-on Legion."""
    facts = {
        "platform": platform.system(),
        "platform_release": platform.release(),
        "is_wsl": "microsoft" in platform.release().lower(),
        "python": sys.version,
        "python_major_minor": list(sys.version_info[:2]),
        "python_executable": sys.executable,
        "environment_prefix": sys.prefix,
        "isolated_env": (
            sys.prefix != sys.base_prefix
            and Path(sys.prefix).name.casefold() == ENV_NAME
        ),
        "packages": {},
    }
    # Refuse the analysis environment or an arbitrary interpreter before hardware queries.
    if not facts["isolated_env"] or facts["python_major_minor"] != [3, 12]:
        facts["error"] = f"Use a fresh Python 3.12 venv named {ENV_NAME}; do not use Phase 1 .venv"
        return facts
    try:
        for package in ["torch", "numpy", "boltz", "rdkit", "biopython", "pytorch-lightning"]:
            try:
                facts["packages"][package] = metadata.version(package)
            except metadata.PackageNotFoundError:
                facts["packages"][package] = None
        smi = find_smi()
        facts["nvidia_smi_path"] = smi
        version = subprocess.run([smi, "--version"], capture_output=True, text=True, timeout=20)
        facts["nvidia_smi_version_output"] = version.stdout.strip()
        query = "index,name,driver_version,compute_cap,memory.total,memory.used,memory.free,power.limit"
        result = subprocess.run(
            [smi, f"--query-gpu={query}", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=20,
        )
        with_capability = result.returncode == 0
        if result.returncode:
            # Some nvidia-smi releases omit compute_cap; Torch still checks it.
            query = query.replace("compute_cap,", "")
            result = subprocess.run(
                [smi, f"--query-gpu={query}", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=20,
            )
        if result.returncode:
            detail = (result.stderr or result.stdout).strip()[:1000]
            raise RuntimeError("nvidia-smi failed: " + detail)
        rows = parse_smi(result.stdout, with_capability=with_capability)
        facts["nvidia_gpus"] = rows
        if len(rows) != 1 or not is_target_gpu(rows[0]["name"]):
            raise RuntimeError("Expected one RTX 5070 Ti Laptop GPU; no CUDA operation was attempted")
        facts["gpu"] = rows[0]

        import torch  # Deferred: importing this script never imports Torch.

        facts.update({
            "torch_version": str(torch.__version__),
            "torch_cuda_build": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "torch_arch_list": torch.cuda.get_arch_list(),
        })
        if not facts["cuda_available"]:
            return facts
        facts["cuda_device_count"] = torch.cuda.device_count()
        if facts["cuda_device_count"] != 1:
            return facts
        facts["torch_gpu_name"] = torch.cuda.get_device_name(0)
        facts["compute_capability"] = list(torch.cuda.get_device_capability(0))
        facts["bf16_supported"] = torch.cuda.is_bf16_supported()
        free, total = torch.cuda.mem_get_info(0)
        facts["torch_memory_before_probe"] = {
            "free_bytes": free, "total_bytes": total,
            "allocated_bytes": torch.cuda.memory_allocated(0),
            "reserved_bytes": torch.cuda.memory_reserved(0),
        }
        if (
            facts["torch_version"] != TORCH_PIN or facts["torch_cuda_build"] != CUDA_PIN
            or facts["compute_capability"] != [12, 0]
            or "sm_120" not in facts["torch_arch_list"] or not facts["bf16_supported"]
        ):
            return facts  # Do not launch kernels on an unreviewed build/architecture.
        if free < GIB:
            raise RuntimeError("Less than 1 GiB free after CUDA initialization; tensor probe skipped")
        torch.cuda.reset_peak_memory_stats(0)
        with torch.no_grad():
            vector = torch.arange(4, dtype=torch.float32, device="cuda:0")
            doubled = vector * 2
            matrix = torch.ones((16, 16), dtype=torch.bfloat16, device="cuda:0")
            product = matrix @ matrix
            torch.cuda.synchronize()
            facts["tensor_test_passed"] = (
                doubled.cpu().tolist() == [0.0, 2.0, 4.0, 6.0]
                and bool(torch.all(product == 16).item())
            )
        free_after, total_after = torch.cuda.mem_get_info(0)
        facts["torch_memory_after_probe"] = {
            "free_bytes": free_after, "total_bytes": total_after,
            "allocated_bytes": torch.cuda.memory_allocated(0),
            "reserved_bytes": torch.cuda.memory_reserved(0),
            "probe_peak_allocated_bytes": torch.cuda.max_memory_allocated(0),
            "probe_peak_reserved_bytes": torch.cuda.max_memory_reserved(0),
        }
        # These peaks describe ONLY this tiny preflight process, not a model.
    except Exception as exc:
        facts["error"] = f"{type(exc).__name__}: {exc}"
    return facts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New JSON report; existing files are refused")
    parser.add_argument("--min-free-gib", type=float, default=9.0, help="Pilot headroom policy, NOT measured model demand")
    parser.add_argument("--require-boltz", action="store_true", help="Also check installed Boltz pins and the Linux pilot requirement")
    args = parser.parse_args(argv)
    if not math.isfinite(args.min_free_gib) or not 1 <= args.min_free_gib <= 11:
        parser.error("--min-free-gib must be between 1 and 11")
    if args.output.exists() or not args.output.parent.is_dir():
        parser.error("Choose a new report filename in an existing directory")
    facts = collect_facts()
    report = assess(facts, args.min_free_gib, args.require_boltz)
    report.update({
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "observations": facts,
        "target_hardware": "Legion Pro 7 Gen 10 / RTX 5070 Ti Laptop / 12 GB VRAM / 64 GB RAM (user-specified)",
        "notes": [
            "No model/checkpoint, MSA search, Boltz import or inference is performed.",
            "nvidia-smi CUDA version is driver capability; torch.version.cuda identifies the wheel runtime.",
            "The 9 GiB default is a conservative starting policy, not an experimentally measured Boltz requirement.",
            "A passing tiny operation does not prove all Boltz operations or the full-MSA input fit in memory.",
        ],
    })
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(f"{report['status']}: {args.output}")
    print("This is a CUDA preflight only; no Boltz inference was run.")
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
