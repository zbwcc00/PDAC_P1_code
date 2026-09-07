"""Targeted GTEx v8 vascular-tissue significant cis-eQTL pre-screen."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import requests


BASE_URL = "https://gtexportal.org/api/v2"
GENES = ["EPAS1", "TACC1", "MARCKS", "HERPUD1", "AQP7"]
TISSUES = ["Artery_Aorta", "Artery_Coronary", "Artery_Tibial"]
OUTPUT_DIR = Path("D:/PDAC_P1/analysis/11_methodology_review/GTEx_v8_targeted")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def get_json(path: str, params: dict[str, str | int]) -> dict:
    response = requests.get(f"{BASE_URL}{path}", params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def main() -> None:
    gencode_ids = {}
    for gene in GENES:
        data = get_json("/reference/gene", {"geneId": gene, "datasetId": "gtex_v8"})["data"]
        gencode_ids[gene] = data[0]["gencodeId"] if data else ""
    all_rows = []
    summary = []
    for gene, gencode_id in gencode_ids.items():
        for tissue in TISSUES:
            data = get_json(
                "/association/singleTissueEqtl",
                {
                    "gencodeId": gencode_id,
                    "tissueSiteDetailId": tissue,
                    "datasetId": "gtex_v8",
                    "itemsPerPage": 100000,
                },
            )["data"]
            summary.append({"gene": gene, "gencode_id": gencode_id, "tissue": tissue, "n_significant_eqtls": len(data)})
            all_rows.extend(data)
    with (OUTPUT_DIR / "GTEx_v8_vascular_candidate_eQTL_prescreen.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "gencode_id", "tissue", "n_significant_eqtls"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    with (OUTPUT_DIR / "GTEx_v8_vascular_candidate_significant_eQTLs.json").open("w", encoding="utf-8") as handle:
        json.dump(all_rows, handle, indent=2)
    with (OUTPUT_DIR / "GTEx_v8_vascular_candidate_significant_eQTLs.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["geneSymbol", "gencodeId", "variantId", "snpId", "chromosome", "pos", "nes", "pValue", "tissueSiteDetailId", "datasetId"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)
    positive = [row for row in summary if row["n_significant_eqtls"] > 0]
    lines = [
        "# GTEx v8 vascular-tissue targeted cis-eQTL pre-screen",
        "",
        "Pre-specified tissues: aorta, coronary artery, and tibial artery. This is a tissue-relevant instrument availability screen, not endothelial cell-specific eQTL evidence.",
        "",
    ]
    lines.extend([f"- {r['gene']} / {r['tissue']}: {r['n_significant_eqtls']} significant eQTLs." for r in summary])
    lines.extend([
        "",
        f"Positive gene-tissue combinations: {len(positive)} of {len(summary)}.",
        "Full beta/SE/EAF summary statistics and a colocalization analysis remain mandatory before causal interpretation.",
    ])
    (OUTPUT_DIR / "GTEx_v8_vascular_candidate_eQTL_prescreen_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

