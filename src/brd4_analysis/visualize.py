"""Structural visualisations for the Phase 1 baseline (matplotlib only).

Produces:
  * per-residue minimum ligand heavy-atom distance bar charts with the
    3.5 / 4.0 / 5.0 Angstrom thresholds drawn as reference lines;
  * 2D depictions of each ligand from its CCD reference structure.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from brd4_analysis.constants import CONTACT_THRESHOLDS


def plot_residue_contacts(contacts_df: pd.DataFrame, pdb_id: str, out_path: Path) -> Path:
    df = contacts_df[contacts_df["within_5.0A"]].sort_values("min_distance_A")
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(df)), 4.5))
    colors = []
    for d in df["min_distance_A"]:
        if d <= 3.5:
            colors.append("#b2182b")
        elif d <= 4.0:
            colors.append("#ef8a62")
        else:
            colors.append("#67a9cf")
    ax.bar(df["residue"], df["min_distance_A"], color=colors, edgecolor="black", lw=0.4)
    for t, style in zip(CONTACT_THRESHOLDS, ["-", "--", ":"]):
        ax.axhline(t, color="grey", ls=style, lw=0.9)
        ax.text(len(df) - 0.4, t + 0.05, f"{t} Å", fontsize=8, color="grey", ha="right")
    ax.set_ylabel("min heavy-atom distance to ligand (Å)")
    ax.set_xlabel("protein residue")
    ax.set_title(f"{pdb_id}: residues within 5.0 Å of bound ligand")
    ax.set_ylim(0, 5.5)
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=7)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def save_ligand_2d(mol: Chem.Mol, comp_id: str, pdb_id: str, out_path: Path) -> Path:
    img = Draw.MolToImage(mol, size=(360, 300), legend=f"{comp_id} (CCD) — {pdb_id}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(out_path)
    return out_path
