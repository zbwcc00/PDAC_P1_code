"""Create publication-style EPAS1 PAS-B docking-pose visualizations from final Vina outputs.

The dashed lines denote geometry-based candidate polar contacts (heavy-atom distance
2.1-3.5 Å), not experimentally verified hydrogen bonds.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


BASE = Path(r"D:/PDAC_P1")
DOCK = BASE / "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer"
FINAL = DOCK / "final_exhaustiveness16"
OUT = BASE / "analysis/15_priority_figures"
OUT.mkdir(parents=True, exist_ok=True)

COMPLEXES = [
    {
        "compound": "Y-39983",
        "pdb_id": "6D0B",
        "affinity": -7.729,
        "ligand_file": FINAL / "6D0B_BRD-K56751279_out.pdbqt",
        "receptor_file": DOCK / "6D0B_receptor.pdb",
        "ligand_colour": "#00A6D6",
        "panel": "A",
    },
    {
        "compound": "triclabendazole",
        "pdb_id": "6CZW",
        "affinity": -6.967,
        "ligand_file": FINAL / "6CZW_BRD-K81916719_out.pdbqt",
        "receptor_file": DOCK / "6CZW_receptor.pdb",
        "ligand_colour": "#C43C8C",
        "panel": "C",
    },
]

ELEMENT_COLOURS = {"C": "#414141", "N": "#2E75B6", "O": "#D62728", "S": "#E6B800", "CL": "#41AB5D", "F": "#41AB5D", "H": "#EEEEEE"}
COVALENT_RADII = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "CL": 1.02, "F": 0.57, "H": 0.31}
POLAR = {"N", "O", "S"}
HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TYR", "TRP", "PRO"}


def infer_element(atom_name: str, element_field: str) -> str:
    element = element_field.strip().upper()
    if element:
        return element
    atom_name = atom_name.strip().upper()
    if atom_name.startswith("CL"):
        return "CL"
    return atom_name[0]


def parse_receptor(path: Path) -> list[dict]:
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atom = line[12:16].strip()
            atoms.append({
                "atom": atom,
                "resname": line[17:20].strip(),
                "chain": line[21].strip() or "_",
                "resseq": int(line[22:26]),
                "coord": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                "element": infer_element(atom, line[76:78]),
            })
        except (ValueError, IndexError):
            continue
    return atoms


def parse_first_pose(path: Path) -> list[dict]:
    atoms, started = [], False
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MODEL"):
            if started:
                break
            started = True
            continue
        if not started or not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atom = line[12:16].strip()
            atoms.append({
                "atom": atom,
                "serial": int(line[6:11]),
                "coord": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                "element": infer_element(atom, line[77:79]),
            })
        except (ValueError, IndexError):
            continue
    return atoms


def residue_key(atom: dict) -> tuple[str, str, int]:
    return atom["chain"], atom["resname"], atom["resseq"]


def residue_label(key: tuple[str, str, int]) -> str:
    return f"{key[1]}{key[2]}"


def select_pocket(receptor: list[dict], ligand: list[dict], cutoff: float = 4.5) -> tuple[list[dict], list[tuple[str, str, int]]]:
    ligand_coords = np.array([a["coord"] for a in ligand if a["element"] != "H"])
    residue_atoms = defaultdict(list)
    for atom in receptor:
        residue_atoms[residue_key(atom)].append(atom)
    distance_by_residue = {}
    for key, atoms in residue_atoms.items():
        coords = np.array([a["coord"] for a in atoms if a["element"] != "H"])
        if len(coords):
            distance_by_residue[key] = float(np.min(np.linalg.norm(coords[:, None, :] - ligand_coords[None, :, :], axis=2)))
    chosen = [key for key, value in sorted(distance_by_residue.items(), key=lambda item: item[1]) if value <= cutoff][:12]
    pocket = [atom for atom in receptor if residue_key(atom) in chosen]
    return pocket, chosen


def polar_contacts(receptor: list[dict], ligand: list[dict], max_distance: float = 3.5) -> list[dict]:
    contacts = []
    for ligand_atom in ligand:
        if ligand_atom["element"] not in POLAR:
            continue
        for receptor_atom in receptor:
            if receptor_atom["element"] not in POLAR:
                continue
            distance = float(np.linalg.norm(ligand_atom["coord"] - receptor_atom["coord"]))
            if 2.1 <= distance <= max_distance:
                contacts.append({
                    "ligand_atom": ligand_atom["atom"],
                    "residue": residue_label(residue_key(receptor_atom)),
                    "receptor_atom": receptor_atom["atom"],
                    "distance_A": distance,
                    "ligand_coord": ligand_atom["coord"],
                    "receptor_coord": receptor_atom["coord"],
                })
    contacts.sort(key=lambda item: item["distance_A"])
    kept, seen = [], set()
    for contact in contacts:
        key = (contact["ligand_atom"], contact["residue"])
        if key not in seen:
            kept.append(contact)
            seen.add(key)
        if len(kept) == 4:
            break
    return kept


def bonds(atoms: list[dict]) -> list[tuple[int, int]]:
    pairs = []
    for first in range(len(atoms)):
        if atoms[first]["element"] == "H":
            continue
        for second in range(first + 1, len(atoms)):
            if atoms[second]["element"] == "H":
                continue
            distance = float(np.linalg.norm(atoms[first]["coord"] - atoms[second]["coord"]))
            cutoff = COVALENT_RADII.get(atoms[first]["element"], 0.77) + COVALENT_RADII.get(atoms[second]["element"], 0.77) + 0.35
            if 0.7 < distance <= cutoff:
                pairs.append((first, second))
    return pairs


def set_equal_3d(ax, coords: np.ndarray, padding: float = 2.0) -> None:
    center = coords.mean(axis=0)
    radius = max(np.ptp(coords, axis=0).max() / 2, 1.0) + padding
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1))


def style_axis(ax) -> None:
    ax.set_axis_off()
    ax.grid(False)
    ax.set_facecolor("white")


def protein_cartoon(ax, receptor: list[dict], ligand: list[dict], ligand_colour: str, title: str, label: str) -> None:
    by_chain = defaultdict(list)
    for atom in receptor:
        if atom["atom"] == "CA":
            by_chain[atom["chain"]].append(atom)
    for chain_atoms in by_chain.values():
        chain_atoms = sorted(chain_atoms, key=lambda item: item["resseq"])
        coordinates = np.array([atom["coord"] for atom in chain_atoms])
        ax.plot(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], color="#152D62", linewidth=3.3, alpha=0.92, solid_capstyle="round")
    ligand_coords = np.array([atom["coord"] for atom in ligand if atom["element"] != "H"])
    for first, second in bonds(ligand):
        start, end = ligand[first]["coord"], ligand[second]["coord"]
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=ligand_colour, linewidth=3.0, zorder=5)
    ax.scatter(ligand_coords[:, 0], ligand_coords[:, 1], ligand_coords[:, 2], s=34, color=ligand_colour, edgecolor="white", linewidth=0.35, depthshade=False, zorder=6)
    protein_coords = np.array([atom["coord"] for atom in receptor if atom["atom"] == "CA"])
    set_equal_3d(ax, np.vstack((protein_coords, ligand_coords)), padding=1.5)
    ax.view_init(elev=19, azim=-54)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=-3)
    ax.text2D(0.01, 0.96, label, transform=ax.transAxes, fontsize=14, fontweight="bold")
    style_axis(ax)


def pocket_view(ax, pocket: list[dict], residues: list[tuple[str, str, int]], ligand: list[dict], contacts: list[dict], ligand_colour: str, title: str, label: str) -> None:
    for first, second in bonds(pocket):
        start, end = pocket[first]["coord"], pocket[second]["coord"]
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color="#7D8FA5", linewidth=1.1, alpha=0.85)
    pocket_coords = np.array([atom["coord"] for atom in pocket if atom["element"] != "H"])
    pocket_colours = ["#B7824A" if atom["resname"] in HYDROPHOBIC_RESIDUES else "#97AABC" for atom in pocket if atom["element"] != "H"]
    ax.scatter(pocket_coords[:, 0], pocket_coords[:, 1], pocket_coords[:, 2], s=12, color=pocket_colours, alpha=0.85, depthshade=False)
    for first, second in bonds(ligand):
        start, end = ligand[first]["coord"], ligand[second]["coord"]
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=ligand_colour, linewidth=3.2, zorder=5)
    for atom in ligand:
        if atom["element"] == "H":
            continue
        ax.scatter(*atom["coord"], s=42, color=ELEMENT_COLOURS.get(atom["element"], ligand_colour), edgecolor="black", linewidth=0.25, depthshade=False, zorder=6)
    for contact in contacts:
        start, end = contact["ligand_coord"], contact["receptor_coord"]
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], linestyle=(0, (3, 2)), color="#E27D1D", linewidth=1.6, zorder=8)
        midpoint = (start + end) / 2
        ax.text(*midpoint, f"{contact['distance_A']:.1f}", fontsize=6.5, color="#9A4F00")
    residue_atoms = defaultdict(list)
    for atom in pocket:
        residue_atoms[residue_key(atom)].append(atom)
    for key in residues[:7]:
        atoms = residue_atoms[key]
        anchor = next((atom["coord"] for atom in atoms if atom["atom"] == "CA"), atoms[0]["coord"])
        ax.text(*(anchor + np.array([0.35, 0.35, 0.35])), residue_label(key), fontsize=7, color="#263746")
    ligand_coords = np.array([atom["coord"] for atom in ligand if atom["element"] != "H"])
    set_equal_3d(ax, np.vstack((pocket_coords, ligand_coords)), padding=1.0)
    ax.view_init(elev=18, azim=-56)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=-3)
    ax.text2D(0.01, 0.96, label, transform=ax.transAxes, fontsize=14, fontweight="bold")
    polar_text = "polar contacts: " + ", ".join(f"{item['residue']} ({item['distance_A']:.2f} A)" for item in contacts) if contacts else "no candidate polar contact <=3.5 A"
    ax.text2D(0.04, 0.04, polar_text, transform=ax.transAxes, fontsize=7.2, color="#9A4F00" if contacts else "#555555")
    style_axis(ax)


def main() -> None:
    figure = plt.figure(figsize=(10.0, 8.4), dpi=300)
    contact_rows = []
    axes = [figure.add_subplot(2, 2, index, projection="3d") for index in range(1, 5)]
    for offset, complex_info in enumerate(COMPLEXES):
        receptor = parse_receptor(complex_info["receptor_file"])
        ligand = parse_first_pose(complex_info["ligand_file"])
        pocket, residues = select_pocket(receptor, ligand)
        contacts = polar_contacts(pocket, ligand)
        protein_cartoon(axes[offset * 2], receptor, ligand, complex_info["ligand_colour"], f"{complex_info['compound']} in EPAS1 PAS-B ({complex_info['pdb_id']})", complex_info["panel"])
        pocket_view(axes[offset * 2 + 1], pocket, residues, ligand, contacts, complex_info["ligand_colour"], f"Pocket view; Vina = {complex_info['affinity']:.3f} kcal/mol", chr(ord(complex_info["panel"]) + 1))
        for contact in contacts:
            contact_rows.append({
                "compound": complex_info["compound"], "representative_pdb": complex_info["pdb_id"], "vina_affinity_kcal_mol": complex_info["affinity"],
                "ligand_atom": contact["ligand_atom"], "residue": contact["residue"], "receptor_atom": contact["receptor_atom"],
                "distance_A": round(contact["distance_A"], 3), "annotation": "geometry_based_candidate_polar_contact",
            })
    figure.suptitle("Supplementary Figure S23. Representative EPAS1 PAS-B docking binding modes", fontsize=15, fontweight="bold", y=0.985)
    figure.text(0.5, 0.945, "Best-scoring pose per compound from the final exhaustiveness=16 five-conformer screen. Orange dashed lines: candidate polar contacts (heavy-atom distance 2.1-3.5 A).", ha="center", fontsize=8.8, color="#444444")
    legend = [Line2D([0], [0], color="#152D62", lw=3, label="EPAS1 backbone"), Line2D([0], [0], color="#00A6D6", lw=3, label="Y-39983 pose"), Line2D([0], [0], color="#C43C8C", lw=3, label="triclabendazole pose"), Line2D([0], [0], color="#E27D1D", lw=1.5, linestyle=(0, (3, 2)), label="candidate polar contact")]
    figure.legend(handles=legend, loc="lower center", ncol=4, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.012))
    figure.text(0.5, 0.04, "Visualization is structural hypothesis support only; dashed contacts are not experimentally validated hydrogen bonds or direct-binding proof.", ha="center", fontsize=8, color="#555555")
    figure.subplots_adjust(left=0.02, right=0.98, top=0.91, bottom=0.09, wspace=0.02, hspace=0.04)
    stem = OUT / "Figure_S23_EPAS1_complex_binding_modes"
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    pd.DataFrame(contact_rows).to_csv(OUT / "S23_EPAS1_representative_pose_polar_contacts.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
