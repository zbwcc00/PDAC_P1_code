from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


BASE = Path(r"D:\第二篇大论文")
DRUG = BASE / "analysis" / "13_virtual_perturbation" / "drugreflector"
META = BASE / "data_external" / "drugreflector" / "compoundinfo_beta.txt"


def read_table(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    priority = read_table(DRUG / "drugreflector_priority.tsv")
    overlap = read_table(DRUG / "drugreflector_cross_axis_overlap.tsv")
    metadata = read_table(META)
    by_id = defaultdict(list)
    for row in metadata:
        by_id[row.get("pert_id", "")].append(row)

    selected = {row["compound"] for row in priority}
    selected.update(row["compound"] for row in overlap)
    output = []
    for compound in sorted(selected):
        records = by_id.get(compound, [])
        if records:
            names = sorted({r.get("cmap_name", "") for r in records if r.get("cmap_name")})
            targets = sorted({r.get("target", "") for r in records if r.get("target")})
            moas = sorted({r.get("moa", "") for r in records if r.get("moa")})
            smiles = sorted({r.get("canonical_smiles", "") for r in records if r.get("canonical_smiles")})
            aliases = sorted({r.get("compound_aliases", "") for r in records if r.get("compound_aliases")})
            inchis = sorted({r.get("inchi_key", "") for r in records if r.get("inchi_key")})
            resolved = any(name and name != compound for name in names)
            output.append(
                {
                    "compound": compound,
                    "cmap_name": ";".join(names),
                    "target": ";".join(targets),
                    "moa": ";".join(moas),
                    "canonical_smiles": smiles[0] if smiles else "",
                    "inchi_key": inchis[0] if inchis else "",
                    "compound_aliases": ";".join(aliases),
                    "mapping_status": "name_resolved" if resolved else "structure_only",
                    "mapping_source": "LINCS/CLUE compoundinfo_beta.txt",
                }
            )
        else:
            output.append(
                {
                    "compound": compound,
                    "cmap_name": "",
                    "target": "",
                    "moa": "",
                    "canonical_smiles": "",
                    "inchi_key": "",
                    "compound_aliases": "",
                    "mapping_status": "not_in_metadata",
                    "mapping_source": "LINCS/CLUE compoundinfo_beta.txt",
                }
            )
    fields = [
        "compound",
        "cmap_name",
        "target",
        "moa",
        "canonical_smiles",
        "inchi_key",
        "compound_aliases",
        "mapping_status",
        "mapping_source",
    ]
    out = DRUG / "drugreflector_compound_mapping.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(f"Selected compounds: {len(selected)}")
    for status in sorted({row["mapping_status"] for row in output}):
        print(status, sum(row["mapping_status"] == status for row in output))


if __name__ == "__main__":
    main()
