from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path


BASE = Path(r"D:\PDAC_P1")
D = BASE / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer"
VINA = BASE / "software" / "autodock_vina" / "vina_1.2.7_win.exe"
PDB_IDS = ["5TBM", "6CZW", "6D09", "6D0B", "6X21"]
CANDIDATES = ["BRD-K56751279", "BRD-K81916719"]
CENTER = (23.974, -0.251, -10.571)
SIZE = (20.0, 20.0, 20.0)
OUT = D / "final_exhaustiveness16"
OUT.mkdir(parents=True, exist_ok=True)


def run_one(pdb_id: str, compound: str):
    receptor = D / f"{pdb_id}.pdbqt"
    ligand = D.parent / "EPAS1" / "ligands_pdbqt" / f"{compound}.pdbqt"
    pose = OUT / f"{pdb_id}_{compound}_out.pdbqt"
    log = OUT / f"{pdb_id}_{compound}.log"
    command = [
        str(VINA), "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]), "--size_x", str(SIZE[0]),
        "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
        "--exhaustiveness", "16", "--num_modes", "9", "--energy_range", "10", "--cpu", "8",
        "--out", str(pose), "--verbosity", "1",
    ]
    with log.open("w", encoding="utf-8") as handle:
        process = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)
    text = pose.read_text(encoding="utf-8", errors="ignore") if pose.exists() else ""
    scores = [float(value) for value in re.findall(r"REMARK VINA RESULT:\s+([-0-9.]+)", text)]
    return {
        "pdb_id": pdb_id,
        "compound": compound,
        "affinity_best_kcal_mol": min(scores) if scores else "",
        "affinity_all_modes_kcal_mol": ";".join(f"{value:.3f}" for value in scores),
        "n_modes": len(scores),
        "exhaustiveness": 16,
        "num_modes_requested": 9,
        "energy_range": 10,
        "cpu": 8,
        "pose_file": str(pose),
        "log_file": str(log),
        "return_code": process.returncode,
        "status": "completed" if scores else "failed",
    }


def main():
    rows = []
    total = len(PDB_IDS) * len(CANDIDATES)
    count = 0
    for pdb_id in PDB_IDS:
        for compound in CANDIDATES:
            count += 1
            print(f"[{count}/{total}] {pdb_id} {compound}", flush=True)
            row = run_one(pdb_id, compound)
            rows.append(row)
            print(f"  {row['status']} best={row['affinity_best_kcal_mol']}", flush=True)
    out = OUT / "final_exhaustiveness16_summary.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
