from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem


BASE = Path(r"D:\PDAC_P1")
out = BASE / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer"
smiles = "O=C1NC(=O)C(=C1c2c[nH]c3ccccc23)c4c[nH]c5ccccc45"
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
AllChem.EmbedMolecule(mol, randomSeed=20260903)
AllChem.UFFOptimizeMolecule(mol, maxIters=500)
sdf = out / "BRD-K49448285_negative.sdf"
writer = Chem.SDWriter(str(sdf)); writer.write(mol); writer.close()
print(sdf)
