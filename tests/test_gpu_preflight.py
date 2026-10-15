"""Synthetic decision/parser tests only; never execute hardware collection."""

import builtins
import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/boltz_gpu_preflight.py"


@pytest.fixture
def preflight(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Synthetic tests must not query hardware or execute commands")

    monkeypatch.setattr("subprocess.run", forbidden)
    original_import = builtins.__import__

    def no_torch(name, *args, **kwargs):
        if name == "torch" or name.startswith("torch."):
            raise AssertionError("Synthetic tests must not import Torch")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_torch)
    spec = importlib.util.spec_from_file_location("synthetic_preflight", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_facts():
    return {
        "platform": "Windows", "is_wsl": False, "isolated_env": True,
        "python_major_minor": [3, 12],
        "gpu": {"name": "NVIDIA GeForce RTX 5070 Ti Laptop GPU", "driver_version": "600.00"},
        "torch_version": "2.7.1+cu128", "torch_cuda_build": "12.8",
        "packages": {"numpy": "1.26.4"},
        "cuda_available": True, "cuda_device_count": 1,
        "compute_capability": [12, 0], "torch_arch_list": ["sm_90", "sm_120"],
        "bf16_supported": True, "tensor_test_passed": True,
        "torch_memory_before_probe": {"total_bytes": 12 * 1024**3, "free_bytes": 10 * 1024**3},
    }


def test_smi_csv_accepts_unavailable_laptop_power_telemetry(preflight):
    rows = preflight.parse_smi(
        "0, NVIDIA GeForce RTX 5070 Ti Laptop GPU, 600.00, 12.0, 12227, 1000, 11227, [N/A]\n"
    )
    assert rows[0]["memory_free_mib"] == 11227
    assert rows[0]["power_limit_w"] is None
    with pytest.raises(ValueError, match="column"):
        preflight.parse_smi("0, incomplete\n")


def test_expected_synthetic_stack_passes_without_claiming_model_fit(preflight):
    report = preflight.assess(synthetic_facts(), 9)
    assert report["status"] == "PASS"
    assert report["model_weights_loaded"] is False
    assert report["boltz_inference_executed"] is False


@pytest.mark.parametrize("key,value,failed_check", [
    ("torch_cuda_build", "12.6", "torch_cuda_build"),
    ("torch_arch_list", ["sm_90"], "native_sm120_in_torch_build"),
    ("bf16_supported", False, "bf16_supported"),
    ("isolated_env", False, "dedicated_virtual_environment"),
])
def test_wrong_build_or_environment_is_not_accepted(preflight, key, value, failed_check):
    facts = synthetic_facts()
    facts[key] = value
    report = preflight.assess(facts, 9)
    assert report["status"] == "FAIL_OR_INCOMPLETE"
    assert failed_check in report["failed_checks"]


def test_headroom_policy_is_separate_from_cuda_compatibility(preflight):
    facts = synthetic_facts()
    facts["torch_memory_before_probe"]["free_bytes"] = 8 * 1024**3
    report = preflight.assess(facts, 9)
    assert report["compatibility_checks_pass"]
    assert not report["headroom_policy_pass"]
    assert report["failed_checks"] == ["free_vram_policy"]
    facts["torch_memory_before_probe"]["free_bytes"] = 10 * 1024**3
    facts["torch_memory_after_probe"] = {"free_bytes": 8 * 1024**3}
    assert not preflight.assess(facts, 9)["headroom_policy_pass"]


def test_boltz_stage_requires_linux_and_matching_package_metadata(preflight):
    facts = synthetic_facts()
    facts["packages"].update({"boltz": "2.2.1", "rdkit": "2026.3.6", "biopython": "1.84"})
    assert "current_pilot_linux_requirement" in preflight.assess(facts, 9, True)["failed_checks"]
    facts["platform"] = "Linux"
    facts["is_wsl"] = True
    assert preflight.assess(facts, 9, True)["status"] == "PASS"
