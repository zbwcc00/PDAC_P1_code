"""Targeted GTEx v8 Pancreas significant cis-eQTL pre-screen.

The GTEx Portal API returns significant associations only. It is used here to
assess instrument availability, not as final MR exposure summary statistics.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import requests


BASE_URL = "https://gtexportal.org/api/v2"
GENES = ["EPAS1", "TACC1", "MARCKS", "HERPUD1", "AQP7"]
OUTPUT_DIR = Path("D:/PDAC_P1/analysis/11_methodology_review/GTEx_v8_targeted")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def get_json(path: str, params: dict[str, str | int]) -> dict:
    response = requests.get(f"{BASE_URL}{path}", params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def main() -> None:
    all_rows: list[dict] = []
    summary: list[dict] = []
    for gene in GENES:
        gene_data = get_json("/reference/gene", {"geneId": gene, "datasetId": "gtex_v8"})["data"]
        if not gene_data:
            summary.append({"gene": gene, "gencode_id": "", "n_significant_eqtls": 0, "status": "gene_not_found"})
            continue
        gencode_id = gene_data[0]["gencodeId"]
        eqtl_data = get_json(
            "/association/singleTissueEqtl",
            {
                "gencodeId": gencode_id,
                "tissueSiteDetailId": "Pancreas",
                "datasetId": "gtex_v8",
                "itemsPerPage": 100000,
            },
        )["data"]
        summary.append(
            {
                "gene": gene,
                "gencode_id": gencode_id,
                "n_significant_eqtls": len(eqtl_data),
                "status": "significant_eqtls_found" if eqtl_data else "no_significant_pancreas_eqtl",
            }
        )
        all_rows.extend(eqtl_data)

    with (OUTPUT_DIR / "GTEx_v8_Pancreas_candidate_significant_eQTLs.json").open("w", encoding="utf-8") as handle:
        json.dump(all_rows, handle, indent=2)
    with (OUTPUT_DIR / "GTEx_v8_Pancreas_candidate_eQTL_prescreen.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["gene", "gencode_id", "n_significant_eqtls", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    with (OUTPUT_DIR / "GTEx_v8_Pancreas_candidate_significant_eQTLs.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["geneSymbol", "gencodeId", "variantId", "snpId", "chromosome", "pos", "nes", "pValue", "tissueSiteDetailId", "datasetId"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    report = [
        "# GTEx v8 Pancreas targeted cis-eQTL pre-screen",
        "",
        "This API query returns significant eQTLs only and lacks the full beta/SE/EAF fields required for final harmonised MR.",
        "It is used solely to determine whether a candidate has any significant Pancreas cis-eQTL signal.",
        "",
    ]
    for row in summary:
        report.append(f"- {row['gene']}: {row['n_significant_eqtls']} significant Pancreas eQTLs ({row['status']}).")
    report.extend(
        [
            "",
            "A candidate with no significant Pancreas eQTL should not proceed to Pancreas-tissue MR. It may still be evaluated in a pre-specified vascular tissue if biologically justified.",
            "Candidates passing this pre-screen require full exposure summary statistics with beta, SE and EAF plus PDAC outcome harmonisation before MR/colocalization.",
        ]
    )
    (OUTPUT_DIR / "GTEx_v8_Pancreas_candidate_eQTL_prescreen_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

