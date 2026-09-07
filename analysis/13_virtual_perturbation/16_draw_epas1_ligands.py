from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw


D = Path(r"D:\PDAC_P1\analysis\13_virtual_perturbation\docking\EPAS1_multiconformer")
items = {
    "Y-39983": D.parent / "EPAS1" / "BRD-K56751279.sdf",
    "triclabendazole": D.parent / "EPAS1" / "BRD-K81916719.sdf",
    "negative_control": D / "BRD-K49448285_negative.sdf",
}
mols = []
legends = []
for label, filename in items.items():
    mol = Chem.SDMolSupplier(str(filename), removeHs=True)[0]
    mols.append(mol)
    legends.append(label)
    Draw.MolToFile(mol, str(D / f"{label}_2d.png"), size=(900, 700), legend=label)
img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(800, 650), legends=legends, useSVG=False)
img.save(str(D / "epas1_ligand_2d_panel.png"))
print("wrote 2D ligand panel")
