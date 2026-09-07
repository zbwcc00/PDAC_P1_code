from __future__ import annotations

import csv
from pathlib import Path


DRUG = Path(r"D:\第二篇大论文\analysis\13_virtual_perturbation\drugreflector")


def read(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    priority = read(DRUG / "drugreflector_priority.tsv")
    mapping = {row["compound"]: row for row in read(DRUG / "drugreflector_compound_mapping_pubchem.tsv")}
    cross = read(DRUG / "drugreflector_cross_axis_overlap.tsv")
    selected = {}
    for target in ("EPAS1", "MARCKS", "HERPUD1"):
        target_rows = [row for row in priority if row["target"] == target]
        target_rows.sort(key=lambda row: float(row["priority_score"]))
        for row in target_rows[:5]:
            selected[(target, row["compound"])] = "top5 reverse priority for axis"
    for row in cross:
        if int(row["n_axes"]) == 3:
            for target in row["axes"].split(";"):
                selected[(target, row["compound"])] = "cross-axis recurrence in all three axes"
    output = []
    for (target, compound), reason in selected.items():
        m = mapping.get(compound, {})
        p = next((row for row in priority if row["target"] == target and row["compound"] == compound), {})
        output.append(
            {
                "axis": target,
                "compound": compound,
                "cmap_name": m.get("cmap_name", ""),
                "pubchem_cid": m.get("pubchem_cid", ""),
                "canonical_smiles": m.get("canonical_smiles", ""),
                "inchi_key": m.get("inchi_key", ""),
                "n_signatures": p.get("n_signatures", ""),
                "n_strata": p.get("n_strata", ""),
                "mean_rank": p.get("mean_rank", ""),
                "priority_score": p.get("priority_score", ""),
                "selection_reason": reason,
                "structure_status": m.get("pubchem_status", ""),
            }
        )
    output.sort(key=lambda row: (("EPAS1", "MARCKS", "HERPUD1").index(row["axis"]), float(row["priority_score"] or 999)))
    fields = list(output[0])
    out = DRUG / "drugreflector_docking_shortlist.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(f"Docking shortlist rows: {len(output)}")
    print(f"Unique compounds: {len({row['compound'] for row in output})}")


if __name__ == "__main__":
    main()
