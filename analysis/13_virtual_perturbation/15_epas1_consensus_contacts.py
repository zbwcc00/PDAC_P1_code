from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from meeko import PDBQTMolecule


BASE = Path(r"D:\PDAC_P1")
D = BASE / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer"
PDB_IDS = ["5TBM", "6CZW", "6D09", "6D0B", "6X21"]
NATIVE = {"5TBM": "79A", "6CZW": "FO7", "6D09": "FOJ", "6D0B": "FOV", "6X21": "UKJ"}
CANDIDATES = {"Y-39983": "BRD-K56751279", "triclabendazole": "BRD-K81916719", "negative_control": "BRD-K49448285_negative"}


def receptor_atoms(pdb_id):
    atoms = []
    for line in (D / f"{pdb_id}_receptor.pdb").read_text(errors="ignore").splitlines():
        if line.startswith("ATOM"):
            atoms.append({"coord": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), "residue": f"{line[21:22].strip()}:{line[17:20].strip()}{int(line[22:26])}"})
    return atoms


def pose_data(path):
    pose = PDBQTMolecule.from_file(str(path))
    idx = [i for i, atom in enumerate(pose.atoms()) if atom[2] != "H"]
    coords = np.asarray(pose.positions())[idx]
    affinity = float(pose.score) if pose.score is not None else float("nan")
    return coords, affinity


def contacts_for(pdb_id, pose_path, cutoff=4.0):
    receptor = receptor_atoms(pdb_id)
    ligand, affinity = pose_data(pose_path)
    contact_residues = set()
    for atom in receptor:
        if np.min(np.linalg.norm(ligand - atom["coord"], axis=1)) <= cutoff:
            contact_residues.add(atom["residue"])
    return contact_residues, affinity


def main():
    records = []
    all_contacts = defaultdict(Counter)
    pose_groups = {}
    for label, prefix in [("positive_control", "native")]:
        for pdb_id in PDB_IDS:
            pose = D / f"{pdb_id}_79A_redock_out.pdbqt"
            # Native redock output is the positive control for each EPAS1 construct.
            contacts, affinity = contacts_for(pdb_id, pose)
            pose_groups[(label, pdb_id)] = contacts
            for residue in contacts:
                all_contacts[label][residue] += 1
            records.append({"pdb_id": pdb_id, "group": label, "compound": NATIVE[pdb_id], "affinity_kcal_mol": affinity, "contact_count": len(contacts), "contacts": ";".join(sorted(contacts))})
    for label, compound in CANDIDATES.items():
        for pdb_id in PDB_IDS:
            pose = D / f"{pdb_id}_{compound}_out.pdbqt"
            if label == "negative_control":
                pose = D / f"{pdb_id}_{compound}_out.pdbqt"
            contacts, affinity = contacts_for(pdb_id, pose)
            pose_groups[(label, pdb_id)] = contacts
            for residue in contacts:
                all_contacts[label][residue] += 1
            records.append({"pdb_id": pdb_id, "group": label, "compound": compound, "affinity_kcal_mol": affinity, "contact_count": len(contacts), "contacts": ";".join(sorted(contacts))})

    with (D / "epas1_contact_by_conformer.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(records)

    consensus = []
    for group, counter in all_contacts.items():
        for residue, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
            consensus.append({"group": group, "residue": residue, "n_conformers": count, "frequency": f"{count / len(PDB_IDS):.2f}"})
    with (D / "epas1_consensus_binding_residues.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(consensus[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(consensus)

    groups = ["positive_control", "Y-39983", "triclabendazole", "negative_control"]
    residue_counts = sorted({r["residue"] for r in consensus}, key=lambda res: (-max((x["n_conformers"] for x in consensus if x["residue"] == res), default=0), res))[:30]
    matrix = np.zeros((len(groups), len(residue_counts)))
    for i, group in enumerate(groups):
        for j, residue in enumerate(residue_counts):
            matrix[i, j] = all_contacts[group][residue] / len(PDB_IDS)
    fig, ax = plt.subplots(figsize=(14, 3.8))
    im = ax.imshow(matrix, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(range(len(groups)), groups)
    ax.set_xticks(range(len(residue_counts)), residue_counts, rotation=75, ha="right", fontsize=8)
    ax.set_xlabel("EPAS1 residue (chain:residue)"); ax.set_ylabel("Docking group")
    ax.set_title("EPAS1 PAS-B contact frequency across five crystal conformers")
    fig.colorbar(im, ax=ax, label="Frequency")
    fig.tight_layout(); fig.savefig(D / "epas1_contact_frequency_heatmap.png", dpi=300); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    affinities = {group: [r["affinity_kcal_mol"] for r in records if r["group"] == group] for group in groups}
    ax.boxplot([affinities[g] for g in groups], labels=groups, showmeans=True)
    ax.set_ylabel("Vina affinity (kcal/mol)"); ax.set_title("EPAS1 docking controls and candidates")
    fig.tight_layout(); fig.savefig(D / "epas1_affinity_controls.png", dpi=300); plt.close(fig)

    # 3D contact schematic in the 5TBM coordinate frame.
    target_res = [res for res in residue_counts if all_contacts["Y-39983"][res] >= 3]
    atoms = receptor_atoms("5TBM")
    xyz = []; labels = []
    for atom in atoms:
        if atom["residue"] in target_res and atom["residue"] not in labels:
            xyz.append(atom["coord"]); labels.append(atom["residue"])
    ligand, _ = pose_data(D / "5TBM_BRD-K56751279_out.pdbqt")
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    fig = plt.figure(figsize=(7, 6)); ax = fig.add_subplot(111, projection="3d")
    if xyz: ax.scatter(np.array(xyz)[:, 0], np.array(xyz)[:, 1], np.array(xyz)[:, 2], s=55, c="crimson", label="consensus residues")
    ax.scatter(ligand[:, 0], ligand[:, 1], ligand[:, 2], s=12, c="royalblue", label="Y-39983 pose")
    ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)"); ax.set_title("5TBM EPAS1 PAS-B: Y-39983 consensus contacts")
    ax.legend(loc="upper left"); fig.tight_layout(); fig.savefig(D / "epas1_Y39983_3d_contact_schematic.png", dpi=300); plt.close(fig)
    print("wrote contact tables and figures")


if __name__ == "__main__":
    main()
