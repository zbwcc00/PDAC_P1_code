from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\第二篇大论文\analysis\13_virtual_perturbation\drugreflector")
INPUT = ROOT / "drugreflector_top50.tsv"


def read_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    rows = read_rows(INPUT)
    for row in rows:
        row["rank_num"] = float(row["rank"])
        row["prob_num"] = float(row["prob"])
        row["target"] = row["target"].strip()
        row["direction"] = row["direction"].strip()

    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["target"], row["direction"], row["compound"])].append(row)

    summary = []
    for (target, direction, compound), items in grouped.items():
        ranks = [item["rank_num"] for item in items]
        probs = [item["prob_num"] for item in items]
        strata = sorted({item["signature"] for item in items})
        summary.append(
            {
                "target": target,
                "direction": direction,
                "compound": compound,
                "n_signatures": len(items),
                "n_strata": len(strata),
                "mean_rank": f"{sum(ranks) / len(ranks):.4f}",
                "median_rank": f"{sorted(ranks)[len(ranks) // 2]:.4f}",
                "best_rank": f"{min(ranks):.4f}",
                "worst_rank": f"{max(ranks):.4f}",
                "mean_probability": f"{sum(probs) / len(probs):.8g}",
                "signature_names": ";".join(strata),
            }
        )

    summary.sort(
        key=lambda row: (
            row["target"],
            row["direction"],
            -int(row["n_signatures"]),
            float(row["mean_rank"]),
        )
    )
    write_rows(
        ROOT / "drugreflector_stability_summary.tsv",
        summary,
        [
            "target",
            "direction",
            "compound",
            "n_signatures",
            "n_strata",
            "mean_rank",
            "median_rank",
            "best_rank",
            "worst_rank",
            "mean_probability",
            "signature_names",
        ],
    )

    # Therapeutic prioritization is based on reverse signatures only. State
    # signatures remain available as mechanistic/phenocopy sensitivity analyses.
    priority = []
    for row in summary:
        if row["direction"] != "reverse":
            continue
        n_signatures = int(row["n_signatures"])
        n_strata = int(row["n_strata"])
        mean_rank = float(row["mean_rank"])
        # Smaller is better; the score rewards recurrence across independent strata.
        score = (mean_rank + 1.0) / (50.0 * max(n_strata, 1))
        if n_signatures >= 2 or n_strata >= 2:
            priority.append(
                {
                    "target": row["target"],
                    "compound": row["compound"],
                    "n_signatures": row["n_signatures"],
                    "n_strata": row["n_strata"],
                    "mean_rank": row["mean_rank"],
                    "worst_rank": row["worst_rank"],
                    "mean_probability": row["mean_probability"],
                    "priority_score": f"{score:.8f}",
                    "evidence_rule": "reverse top50; recurrent in >=2 signatures or strata",
                }
            )
    priority.sort(key=lambda row: (row["target"], float(row["priority_score"])))
    write_rows(
        ROOT / "drugreflector_priority.tsv",
        priority,
        [
            "target",
            "compound",
            "n_signatures",
            "n_strata",
            "mean_rank",
            "worst_rank",
            "mean_probability",
            "priority_score",
            "evidence_rule",
        ],
    )

    # Compounds recurring across the three candidate axes (reverse only).
    by_compound = defaultdict(dict)
    for row in priority:
        by_compound[row["compound"]][row["target"]] = row
    overlap = []
    for compound, targets in by_compound.items():
        if len(targets) < 2:
            continue
        mean_score = sum(float(item["priority_score"]) for item in targets.values()) / len(targets)
        overlap.append(
            {
                "compound": compound,
                "n_axes": len(targets),
                "axes": ";".join(sorted(targets)),
                "mean_priority_score": f"{mean_score:.8f}",
                "axis_details": ";".join(
                    f"{target}:n={targets[target]['n_signatures']},mean_rank={targets[target]['mean_rank']}"
                    for target in sorted(targets)
                ),
            }
        )
    overlap.sort(key=lambda row: (-int(row["n_axes"]), float(row["mean_priority_score"])))
    write_rows(
        ROOT / "drugreflector_cross_axis_overlap.tsv",
        overlap,
        ["compound", "n_axes", "axes", "mean_priority_score", "axis_details"],
    )

    print(f"Input rows: {len(rows)}")
    print(f"Unique target-direction-compound summaries: {len(summary)}")
    print(f"Reverse priority rows: {len(priority)}")
    print(f"Cross-axis recurrent compounds: {len(overlap)}")
    for row in overlap[:10]:
        print(row)


if __name__ == "__main__":
    main()
