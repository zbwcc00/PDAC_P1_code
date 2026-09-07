from pathlib import Path
import csv
from rdkit import Chem
from rdkit.Chem import AllChem


BASE = Path(r"D:\PDAC_P1")
DRUG = BASE / "analysis" / "13_virtual_perturbation" / "drugreflector"
OUT = BASE / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1"
OUT.mkdir(parents=True, exist_ok=True)


def write_sdf(row):
    smiles = row["canonical_smiles"]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {row['compound']}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260903
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=20260903, useRandomCoords=True)
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    path = OUT / f"{row['compound']}.sdf"
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return path


def main():
    rows = []
    with (DRUG / "drugreflector_docking_shortlist.tsv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["axis"] == "EPAS1" and row["canonical_smiles"]:
                rows.append(row)
    seen = set()
    prepared = []
    for row in rows:
        if row["compound"] in seen:
            continue
        seen.add(row["compound"])
        path = write_sdf(row)
        prepared.append({"compound": row["compound"], "name": row["cmap_name"], "cid": row["pubchem_cid"], "sdf": str(path)})
    with (OUT / "ligand_preparation.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(prepared[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(prepared)
    print(f"prepared ligands: {len(prepared)}")


if __name__ == "__main__":
    main()
