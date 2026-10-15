"""Strict, offline validation. No network or model imports."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

from brd4_analysis.constants import PROJECT_ROOT

BOLTZ_VERSION = "2.2.1"
RDKIT_VERSION = "2026.3.6"
SOURCE_DIR = PROJECT_ROOT / "data" / "phase2" / "sources"
SEQUENCE = (
    "SMNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDM"
    "GTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEE"
)
GEOMETRY_KEYS = {
    "max_bond_length_error_A", "max_bond_angle_error_deg",
    "max_aromatic_plane_deviation_A", "min_nonbonded_vdw_ratio",
    "max_protein_ligand_overlap_A",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path):
    def no_duplicates(pairs):
        data = {}
        for k, v in pairs:
            if k in data:
                raise ValueError(f"Duplicate JSON key: {k}")
            data[k] = v
        return data

    return json.loads(
        Path(path).read_text(encoding="utf-8"), object_pairs_hook=no_duplicates,
        parse_constant=lambda s: (_ for _ in ()).throw(ValueError(f"Non-finite JSON: {s}")),
    )


def write_json(path: Path, value):
    # Exclusive creation is intentional. An experiment is never silently resumed.
    with Path(path).open("xb") as handle:
        handle.write((json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode())


def keys(value, expected, where):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(f"{where}: expected exactly {sorted(expected)}")


def positive(value, name, maximum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} exceeds pilot limit {maximum}")


def source_path(spec: dict, root: Path = PROJECT_ROOT, *, allow_pending_hash: bool = False) -> Path:
    keys(spec, {"path", "sha256"}, "source")
    if not isinstance(spec["path"], str) or "\\" in spec["path"] or ":" in spec["path"]:
        raise ValueError("Source paths must use repository-relative POSIX spelling")
    rel = Path(spec["path"])
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError("Source path escapes the repository")
    resolved = (root / rel).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError("Source path escapes the repository through a symlink")
    if spec["sha256"] is None and allow_pending_hash:
        return resolved
    if not isinstance(spec["sha256"], str) or not re.fullmatch(r"[0-9a-f]{64}", spec["sha256"]):
        raise ValueError("Source requires a SHA-256 digest")
    if not resolved.is_file() or sha256(resolved) != spec["sha256"]:
        raise ValueError(f"Missing or changed source: {rel}")
    return resolved


def validate_a3m(path: Path, sequence: str) -> int:
    records, current = [], None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current is not None:
                records.append("".join(current))
            current = []
        elif line.strip():
            if current is None:
                raise ValueError("A3M sequence before first header")
            current.append(line.strip())
    if current is not None:
        records.append("".join(current))
    if not records or records[0] != sequence:
        raise ValueError("A3M first sequence must exactly match the full 127-residue query")
    for record in records:
        if not re.fullmatch(r"[A-Za-z-]+", record):
            raise ValueError("Invalid A3M characters")
        aligned = "".join(c for c in record if not c.islower())
        if len(aligned) != len(sequence):
            raise ValueError("A3M aligned length differs from query")
    return len(records)


def experiment_role(config: dict) -> str:
    """No-MSA outputs can never be labelled as the primary scientific condition."""
    return "smoke_test" if config["msa"]["mode"] == "single_sequence" else "primary_structural_evaluation"


def msa_is_pinned(config: dict) -> bool:
    msa = config["msa"]
    return msa["mode"] == "single_sequence" or (
        msa["file"]["sha256"] is not None
        and all(value is not None for value in msa["provenance"].values())
    )


def validate_config(
    config: dict, root: Path = PROJECT_ROOT, *, allow_pending_msa: bool = False,
) -> dict:
    keys(config, {
        "schema_version", "experiment_id", "requires_user_approval", "reference",
        "protein", "ligand", "msa", "boltz", "metrics", "resources",
    }, "configuration")
    if type(config["schema_version"]) is not int or config["schema_version"] != 1 or config["requires_user_approval"] is not True:
        raise ValueError("This preparation schema requires explicit user approval before execution")
    if not isinstance(config["experiment_id"], str) or not re.fullmatch(
        r"4men_[a-z0-9][a-z0-9_-]{2,79}", config["experiment_id"]
    ):
        raise ValueError("Invalid experiment_id")
    ref = config["reference"]
    keys(ref, {"pdb", "ccd_sdf", "ccd_cif", "polymer_metadata"}, "reference")
    for spec in ref.values():
        source_path(spec, root)
    protein = config["protein"]
    keys(protein, {"chain_id", "reference_chain", "sequence", "fit_auth_residue_range"}, "protein")
    if protein != {
        "chain_id": "A", "reference_chain": "A", "sequence": SEQUENCE,
        "fit_auth_residue_range": [44, 168],
    }:
        raise ValueError("Pilot construct/chain/fit is fixed to 4MEN SM + UniProt O60885 residues 44-168")
    ligand = config["ligand"]
    keys(ligand, {"chain_id", "ccd_id", "smiles", "expected_inchikey"}, "ligand")
    if ligand["chain_id"] != "B" or ligand["ccd_id"] != "25K" or ligand["expected_inchikey"] != "JAEVUMCKPBRDAJ-UHFFFAOYSA-N":
        raise ValueError("This pilot supports the deposited 25K reference graph only")
    if not isinstance(ligand["smiles"], str) or not ligand["smiles"]:
        raise ValueError("Ligand SMILES is required")
    msa = config["msa"]
    if not isinstance(msa, dict):
        raise ValueError("msa must be an explicit policy object")
    if msa.get("mode") == "single_sequence":
        keys(msa, {"mode"}, "msa")
    elif msa.get("mode") == "local_a3m":
        keys(msa, {"mode", "file", "provenance"}, "msa")
        keys(msa["provenance"], {"generator", "database", "database_release", "retrieved_utc", "server"}, "msa provenance")
        for value in msa["provenance"].values():
            if value is None and allow_pending_msa:
                continue
            if not isinstance(value, str) or not value.strip():
                raise ValueError("MSA provenance is pending; record its origin before preparation")
        msa_path = source_path(msa["file"], root, allow_pending_hash=allow_pending_msa)
        if msa["file"]["sha256"] is not None:
            if validate_a3m(msa_path, protein["sequence"]) < 2:
                raise ValueError("The primary MSA requires the query and at least one alignment hit")
    else:
        raise ValueError("Only explicit single_sequence or pinned local_a3m is supported; no live server")
    boltz = config["boltz"]
    keys(boltz, {
        "version", "rdkit_version", "model", "seed", "recycling_steps", "sampling_steps",
        "diffusion_samples", "max_parallel_samples", "step_scale", "use_potentials",
        "no_kernels", "affinity", "templates", "constraints", "cache_dir", "max_msa_seqs",
    }, "boltz")
    if boltz["version"] != BOLTZ_VERSION or boltz["rdkit_version"] != RDKIT_VERSION or boltz["model"] != "boltz2":
        raise ValueError("Unreviewed Boltz or RDKit version")
    for key in ["seed", "recycling_steps", "sampling_steps", "diffusion_samples", "max_parallel_samples", "max_msa_seqs"]:
        if type(boltz[key]) is not int:
            raise ValueError(f"{key} must be an integer")
    if not 0 <= boltz["seed"] <= 2**32 - 1:
        raise ValueError("Seed must be a uint32")
    if boltz["diffusion_samples"] != 1 or boltz["max_parallel_samples"] != 1:
        raise ValueError("This pilot specification supports only one sample in one seed")
    if boltz["max_msa_seqs"] != 8192:
        raise ValueError("Keep the upstream 8192-sequence MSA limit; a depth ablation is a separate experiment")
    positive(boltz["recycling_steps"], "recycling_steps", 3)
    positive(boltz["sampling_steps"], "sampling_steps", 200)
    positive(boltz["step_scale"], "step_scale", 2)
    if boltz["step_scale"] < 1:
        raise ValueError("step_scale must be in [1, 2]")
    if boltz["use_potentials"] is not False or boltz["affinity"] is not False:
        raise ValueError("Steering and affinity are outside this pilot")
    if boltz["no_kernels"] is not True or boltz["templates"] != [] or boltz["constraints"] != []:
        raise ValueError("Pilot requires no kernels, templates or structural constraints")
    if not isinstance(boltz["cache_dir"], str) or not re.fullmatch(r"/[A-Za-z0-9_./-]+", boltz["cache_dir"]) or ".." in boltz["cache_dir"].split("/"):
        raise ValueError("cache_dir must be an explicit absolute Linux path")
    metrics = config["metrics"]
    keys(metrics, {
        "contact_cutoffs_A", "primary_contact_cutoff_A", "min_contact_precision",
        "min_contact_recall", "max_ligand_rmsd_A", "max_protein_ca_rmsd_A",
        "pocket_cutoff_A", "geometry",
    }, "metrics")
    if metrics["contact_cutoffs_A"] != [3.5, 4.0, 5.0] or metrics["primary_contact_cutoff_A"] != 4.0:
        raise ValueError("Contact thresholds must retain the Phase 1 conventions")
    for name in ["min_contact_precision", "min_contact_recall"]:
        positive(metrics[name], name, 1)
    for name in ["max_ligand_rmsd_A", "max_protein_ca_rmsd_A", "pocket_cutoff_A"]:
        positive(metrics[name], name)
    keys(metrics["geometry"], GEOMETRY_KEYS, "geometry")
    for name, value in metrics["geometry"].items():
        positive(value, name)
    positive(metrics["geometry"]["min_nonbonded_vdw_ratio"], "vdw_ratio", 1)
    resources = config["resources"]
    keys(resources, {"gpu_count", "recommended_vram_GB", "min_ram_GB", "min_free_disk_GB",
                     "wall_time_limit_minutes", "max_cost_USD"}, "resources")
    if type(resources["gpu_count"]) is not int or resources["gpu_count"] != 1:
        raise ValueError("One GPU only")
    for key, value in resources.items():
        positive(value, key)
    positive(resources["wall_time_limit_minutes"], "wall_time_limit_minutes", 30)
    positive(resources["max_cost_USD"], "max_cost_USD", 5)
    return config
