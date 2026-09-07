from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from meeko import PDBQTMolecule


BASE = Path(r"D:\PDAC_P1")
D = BASE / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer"
F = D / "final_exhaustiveness16"
PDB_IDS = ["5TBM", "6CZW", "6D09", "6D0B", "6X21"]
COMPOUNDS = {"Y-39983": "BRD-K56751279", "triclabendazole": "BRD-K81916719"}
ALL_GROUPS = {**COMPOUNDS, "negative_control": "BRD-K49448285_negative"}


def receptor_atoms(pdb_id):
    atoms = []
    for line in (D / f"{pdb_id}_receptor.pdb").read_text(errors="ignore").splitlines():
        if line.startswith("ATOM"):
            atoms.append((np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), f"{line[21:22].strip()}:{line[17:20].strip()}{int(line[22:26])}"))
    return atoms


def load_pose(path):
    molecule = PDBQTMolecule.from_file(str(path))
    indices = [i for i, atom in enumerate(molecule.atoms()) if atom[2] != "H"]
    return np.asarray(molecule.positions())[indices], float(molecule.score)


def contacts(pdb_id, path):
    ligand, score = load_pose(path)
    residues = set()
    for coordinate, residue in receptor_atoms(pdb_id):
        if np.min(np.linalg.norm(ligand - coordinate, axis=1)) <= 4.0:
            residues.add(residue)
    return residues, score


def main():
    records = []
    frequencies = defaultdict(Counter)
    for group, compound in ALL_GROUPS.items():
        for pdb_id in PDB_IDS:
            if group == "negative_control":
                pose = D / f"{pdb_id}_BRD-K49448285_negative_out.pdbqt"
            else:
                pose = F / f"{pdb_id}_{compound}_out.pdbqt"
            residues, score = contacts(pdb_id, pose)
            frequencies[group].update(residues)
            records.append({"group": group, "pdb_id": pdb_id, "compound": compound, "affinity_kcal_mol": score, "contact_count": len(residues), "contacts": ";".join(sorted(residues))})
    with (F / "final_contact_by_conformer.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t"); writer.writeheader(); writer.writerows(records)
    consensus = []
    for group, counter in frequencies.items():
        for residue, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
            consensus.append({"group": group, "residue": residue, "n_conformers": count, "frequency": f"{count / len(PDB_IDS):.2f}"})
    with (F / "final_consensus_binding_residues.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(consensus[0]), delimiter="\t"); writer.writeheader(); writer.writerows(consensus)
    print("final contact tables written")


if __name__ == "__main__":
    main()
