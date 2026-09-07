from pathlib import Path
import csv
import numpy as np
from meeko import PDBQTMolecule


D = Path(r"D:\PDAC_P1\analysis\13_virtual_perturbation\docking\EPAS1_multiconformer")
rows = list(csv.DictReader((D / "multiconformer_summary.tsv").open(encoding="utf-8-sig"), delimiter="\t"))
for row in rows:
    if row["kind"] != "cocrystal_redock":
        continue
    pdb_id = row["pdb_id"]
    ligand = row["ligand"]
    reference = PDBQTMolecule.from_file(str(D / f"{pdb_id}_{ligand}_crystal.pdbqt"))
    pose = PDBQTMolecule.from_file(str(D / f"{pdb_id}_79A_redock_out.pdbqt"))
    reference_idx = [i for i, atom in enumerate(reference.atoms()) if atom[2] != "H"]
    pose_idx = [i for i, atom in enumerate(pose.atoms()) if atom[2] != "H"]
    reference_xyz = np.asarray(reference.positions())[reference_idx]
    pose_xyz = np.asarray(pose.positions())[pose_idx]
    if reference_xyz.shape != pose_xyz.shape:
        row["rmsd_angstrom"] = f"shape_mismatch:{reference_xyz.shape}/{pose_xyz.shape}"
        row["rmsd_method"] = "Meeko-PDBQT heavy-atom order"
        row["qc"] = "fail"
    else:
        value = np.sqrt(np.mean(np.sum((pose_xyz - reference_xyz) ** 2, axis=1)))
        row["rmsd_angstrom"] = f"{value:.4f}"
        row["rmsd_method"] = "Meeko-PDBQT heavy-atom order"
        row["qc"] = "pass" if value <= 2.0 else "fail"

fields = list(rows[0])
if "rmsd_method" not in fields:
    fields.append("rmsd_method")
if "qc" not in fields:
    fields.append("qc")
for row in rows:
    row.setdefault("rmsd_method", "")
    row.setdefault("qc", "")
with (D / "multiconformer_summary_corrected.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    if row["kind"] == "cocrystal_redock":
        print(row["pdb_id"], row["rmsd_angstrom"], row["qc"])
