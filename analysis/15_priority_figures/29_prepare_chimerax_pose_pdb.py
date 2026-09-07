from pathlib import Path


BASE = Path(r"D:/PDAC_P1/analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16")


def convert(source: Path, target: Path) -> None:
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    pose = []
    in_first_model = False
    for line in lines:
        if line.startswith("MODEL"):
            if in_first_model:
                break
            in_first_model = True
            continue
        if in_first_model and line.startswith(("ATOM  ", "HETATM")):
            pose.append(line)
    if not pose:
        raise RuntimeError(f"No first pose found in {source}")
    output = ["REMARK First (best-scoring) Vina pose extracted for ChimeraX visualization", "MODEL        1"]
    for index, line in enumerate(pose, start=1):
        atom_name = line[12:16]
        element = line[77:79].strip() or atom_name.strip()[0]
        output.append(f"HETATM{index:5d} {atom_name:<4s} LIG L   1    {float(line[30:38]):8.3f}{float(line[38:46]):8.3f}{float(line[46:54]):8.3f}  1.00  0.00          {element:>2s}")
    output.extend(["ENDMDL", "END"])
    target.write_text("\n".join(output) + "\n", encoding="ascii")


convert(BASE / "6D0B_BRD-K56751279_out.pdbqt", BASE / "6D0B_Y39983_best_pose.pdb")
convert(BASE / "6CZW_BRD-K81916719_out.pdbqt", BASE / "6CZW_triclabendazole_best_pose.pdb")
