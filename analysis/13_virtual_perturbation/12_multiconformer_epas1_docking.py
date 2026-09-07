from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem


BASE = Path(r"D:\PDAC_P1")
STRUCT = BASE / "data_external" / "structure" / "pdb"
DOCK = BASE / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer"
DOCK.mkdir(parents=True, exist_ok=True)
VINA = BASE / "software" / "autodock_vina" / "vina_1.2.7_win.exe"
PREP_REC = Path(r"C:\Users\HUAWEI\Documents\ChatGPT\Try\python311\Scripts\mk_prepare_receptor.exe")
PREP_LIG = Path(r"C:\Users\HUAWEI\Documents\ChatGPT\Try\python311\Scripts\mk_prepare_ligand.exe")
PDB_IDS = ["5TBM", "6CZW", "6D09", "6D0B", "6X21"]
LIG_IDS = {"5TBM": "79A", "6CZW": "FO7", "6D09": "FOJ", "6D0B": "FOV", "6X21": "UKJ"}
CENTER = (23.974, -0.251, -10.571)
SIZE = (20.0, 20.0, 20.0)
# The high-exhaustiveness multi-conformer stage focuses on the two ligands
# that produced favorable preliminary affinities; the other three remain in
# the low-cost 5TBM screen and are not silently treated as validated hits.
CANDIDATES = ["BRD-K56751279", "BRD-K81916719"]


def atom_lines(path: Path, resname: str):
    lines = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == resname:
            lines.append(line)
    return lines


def make_receptor(pdb_id: str):
    src = STRUCT / f"{pdb_id}.pdb"
    receptor = DOCK / f"{pdb_id}_receptor.pdb"
    lines = [
        line
        for line in src.read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.startswith("ATOM  ") and line[21:22] in {"A", "B"}
    ]
    receptor.write_text("\n".join(lines) + "\nEND\n", encoding="ascii")
    base = DOCK / pdb_id
    cmd = [str(PREP_REC), "--read_pdb", str(receptor), "-o", str(base), "-p", "--delete_bad_res", "--default_altloc", "A", "--write_vina_box", "--box_center", *map(str, CENTER), "--box_size", *map(str, SIZE)]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return base.with_suffix(".pdbqt")


def make_reference_ligand(pdb_id: str):
    resname = LIG_IDS[pdb_id]
    ideal = DOCK / f"{resname}_ideal.sdf"
    if not ideal.exists():
        import urllib.request

        urllib.request.urlretrieve(f"https://files.rcsb.org/ligands/download/{resname}_ideal.sdf", ideal)
    mol = Chem.SDMolSupplier(str(ideal), removeHs=False)[0]
    mol = Chem.RemoveHs(mol)
    coords = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in atom_lines(STRUCT / f"{pdb_id}.pdb", resname)])
    if mol is None or mol.GetNumAtoms() != len(coords):
        raise ValueError(f"atom count mismatch for {pdb_id} {resname}")
    conf = mol.GetConformer()
    for idx, xyz in enumerate(coords):
        from rdkit.Geometry import Point3D

        conf.SetAtomPosition(idx, Point3D(*xyz))
    mol = Chem.AddHs(mol, addCoords=True)
    sdf = DOCK / f"{pdb_id}_{resname}_crystal.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    out = DOCK / f"{pdb_id}_{resname}_crystal.pdbqt"
    subprocess.run([str(PREP_LIG), "-i", str(sdf), "-o", str(out), "--add_index_map"], check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return out, coords


def run_vina(receptor, ligand, out, log, exhaustiveness=8, modes=3):
    cmd = [str(VINA), "--receptor", str(receptor), "--ligand", str(ligand), "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]), "--size_x", str(SIZE[0]), "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]), "--exhaustiveness", str(exhaustiveness), "--num_modes", str(modes), "--cpu", "1", "--out", str(out), "--verbosity", "1"]
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.STDOUT, text=True)
    text = out.read_text(encoding="utf-8", errors="ignore") if out.exists() else ""
    match = re.search(r"REMARK VINA RESULT:\s+([-0-9.]+)", text)
    return float(match.group(1)) if match else float("nan")


def rmsd(reference, pose_pdbqt):
    from meeko import PDBQTMolecule

    pose = PDBQTMolecule.from_file(str(pose_pdbqt))
    atoms = pose.atoms()
    heavy_idx = [i for i, atom in enumerate(atoms) if atom[2] != "H"]
    coords = np.asarray(pose.positions())[heavy_idx]
    if coords.shape != reference.shape:
        return float("nan")
    return float(np.sqrt(np.mean(np.sum((coords - reference) ** 2, axis=1))))


def main():
    receptors = {pdb_id: make_receptor(pdb_id) for pdb_id in PDB_IDS}
    references = {pdb_id: make_reference_ligand(pdb_id) for pdb_id in PDB_IDS}
    rows = []
    for pdb_id in PDB_IDS:
        receptor = receptors[pdb_id]
        ref_ligand, ref_coords = references[pdb_id]
        redock_out = DOCK / f"{pdb_id}_79A_redock_out.pdbqt"
        redock_log = DOCK / f"{pdb_id}_79A_redock.log"
        redock_aff = run_vina(receptor, ref_ligand, redock_out, redock_log, exhaustiveness=4, modes=1)
        rows.append({"pdb_id": pdb_id, "ligand": LIG_IDS[pdb_id], "kind": "cocrystal_redock", "compound": LIG_IDS[pdb_id], "affinity_kcal_mol": redock_aff, "rmsd_angstrom": rmsd(ref_coords, redock_out), "status": "completed"})
        for compound in CANDIDATES:
            ligand = DOCK.parent / "EPAS1" / "ligands_pdbqt" / f"{compound}.pdbqt"
            out = DOCK / f"{pdb_id}_{compound}_out.pdbqt"
            log = DOCK / f"{pdb_id}_{compound}.log"
            affinity = run_vina(receptor, ligand, out, log, exhaustiveness=2, modes=3)
            rows.append({"pdb_id": pdb_id, "ligand": LIG_IDS[pdb_id], "kind": "candidate", "compound": compound, "affinity_kcal_mol": affinity, "rmsd_angstrom": "", "status": "completed"})
    with (DOCK / "multiconformer_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"completed rows: {len(rows)}")


if __name__ == "__main__":
    main()
