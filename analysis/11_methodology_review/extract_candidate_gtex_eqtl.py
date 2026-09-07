from __future__ import annotations

import csv
import gzip
from pathlib import Path

BASE = Path("D:/PDAC_P1")
EQTL_DIR = BASE / "data/04_genetics/cis_eQTL/GTEx_v8/extracted/GTEx_Analysis_v8_eQTL"
OUT = BASE / "analysis/11_methodology_review/GTEx_v8_targeted"
OUT.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "EPAS1": "ENSG00000116016",
    "TACC1": "ENSG00000147526",
    "MARCKS": "ENSG00000277443",
    "HERPUD1": "ENSG00000051108",
    "AQP7": "ENSG00000165269",
}
FILES = {
    "Artery_Aorta": EQTL_DIR / "Artery_Aorta.v8.signif_variant_gene_pairs.txt.gz",
    "Artery_Tibial": EQTL_DIR / "Artery_Tibial.v8.signif_variant_gene_pairs.txt.gz",
    "Pancreas": EQTL_DIR / "Pancreas.v8.signif_variant_gene_pairs.txt.gz",
}


def main() -> None:
    rows = []
    summary = []
    for tissue, path in FILES.items():
        with gzip.open(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                gene_id = row["gene_id"].split(".")[0]
                gene = next((name for name, ensg in TARGETS.items() if ensg == gene_id), None)
                if gene is not None:
                    row["gene"] = gene
                    row["tissue"] = tissue
                    rows.append(row)
        for gene, ensg in TARGETS.items():
            gene_rows = [r for r in rows if r["gene"] == gene and r["tissue"] == tissue]
            summary.append({"gene": gene, "ensembl_id": ensg, "tissue": tissue, "n_significant_pairs": len(gene_rows)})

    fields = ["gene", "tissue", "variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta"]
    with (OUT / "GTEx_v8_targeted_candidate_eQTLs.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    with (OUT / "GTEx_v8_targeted_candidate_eQTL_counts.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "ensembl_id", "tissue", "n_significant_pairs"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    report = [
        "# GTEx v8 candidate eQTL extraction",
        "",
        "The extracted records are significant variant–gene pairs from GTEx v8. They contain nominal slope and slope SE, but not an EAF for every variant; final MR still requires allele-frequency and outcome harmonisation.",
        "",
    ]
    report.extend([f"- {r['gene']} / {r['tissue']}: {r['n_significant_pairs']} significant pairs." for r in summary])
    report.extend(["", "Only EPAS1–Artery_Aorta, TACC1–Artery_Tibial and HERPUD1–Artery_Tibial pass the targeted tissue availability gate."])
    (OUT / "GTEx_v8_targeted_candidate_eQTL_extraction_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
