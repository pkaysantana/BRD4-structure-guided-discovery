"""Project constants for the Phase 1 experimental baseline.

The ligand component IDs below were verified against RCSB non-polymer
entity metadata (data.rcsb.org/rest/v1/core/nonpolymer_entity) on
2026-10-14. All other HET groups present (EDO, DMS, HOH) are treated as
solvent/crystallisation additives, not as the bound ligand.
"""

from pathlib import Path

PDB_IDS: tuple[str, ...] = ("4MEN", "4MEO")

# Actual bound inhibitor per entry (Vidler et al. 2013, doi:10.1021/jm4011302).
LIGAND_COMPONENTS: dict[str, str] = {
    "4MEN": "25K",
    "4MEO": "25V",
}

# Heavy-atom distance thresholds, in Angstrom.
CONTACT_THRESHOLDS: tuple[float, ...] = (3.5, 4.0, 5.0)

# Residues with these component IDs are crystallisation additives / solvent,
# never the bound ligand. Water is handled separately via hetfield "W".
SOLVENT_COMPONENT_IDS: frozenset[str] = frozenset({"EDO", "DMS", "HOH", "SO4", "PO4", "CL", "NA", "MG", "GOL"})

PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
CCD_IDEAL_SDF_URL = "https://files.rcsb.org/ligands/download/{comp_id}_ideal.sdf"
ENTRY_METADATA_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "results"
