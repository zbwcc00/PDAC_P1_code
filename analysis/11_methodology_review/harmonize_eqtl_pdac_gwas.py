from __future__ import annotations

import csv
import gzip
import json
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests


BASE = Path("D:/PDAC_P1")
EQTL = BASE / "data/04_genetics/cis_eQTL/GTEx_v8/extracted/GTEx_Analysis_v8_eQTL"
GWAS = BASE / "data/04_genetics/gwas"
OUT = BASE / "analysis/11_methodology_review/GTEx_v8_targeted"
OUT.mkdir(parents=True, exist_ok=True)
TARGETS = {"EPAS1", "TACC1", "HERPUD1"}


def read_eqtl() -> list[dict]:
    rows: list[dict] = []
    for tissue in ["Artery_Aorta", "Artery_Tibial"]:
        path = EQTL / f"{tissue}.v8.signif_variant_gene_pairs.txt.gz"
        with gzip.open(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                gene_id = row["gene_id"].split(".")[0]
                if (tissue == "Artery_Aorta" and gene_id == "ENSG00000116016") or (
                    tissue == "Artery_Tibial" and gene_id in {"ENSG00000147526", "ENSG00000051108"}
                ):
                    row["tissue"] = tissue
                    row["gene"] = {"ENSG00000116016": "EPAS1", "ENSG00000147526": "TACC1", "ENSG00000051108": "HERPUD1"}[gene_id]
                    rows.append(row)
    return rows


def map_b38_to_b37(variant_id: str) -> dict | None:
    chrom, pos, ref, alt, _ = variant_id.split("_")
    url = f"https://rest.ensembl.org/map/human/GRCh38/{chrom.replace('chr', '')}:{pos}..{pos}/GRCh37"
    response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
    response.raise_for_status()
    mappings = response.json().get("mappings", [])
    if not mappings:
        return None
    mapped = mappings[0]["mapped"]
    return {"chrom37": str(mapped["seq_region_name"]), "pos37": int(mapped["start"]), "ref38": ref, "alt38": alt}


def index_eqtl(rows: list[dict]) -> tuple[dict[tuple[str, int], list[dict]], dict]:
    by_coord: dict[tuple[str, int], list[dict]] = {}
    mappings = {}
    tasks = {i: row for i, row in enumerate(rows) if row["variant_id"].endswith("_b38")}
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(map_b38_to_b37, row["variant_id"]): i for i, row in tasks.items()}
        for future in as_completed(futures):
            i = futures[future]
            row = tasks[i]
            try:
                mapping = future.result()
            except Exception as exc:
                row["mapping_error"] = str(exc)
                continue
            if mapping:
                row.update(mapping)
                by_coord.setdefault((mapping["chrom37"], mapping["pos37"]), []).append(row)
                mappings[row["variant_id"]] = mapping
    return by_coord, mappings


def scan_gwas(path: Path, study: str, by_coord: dict[tuple[str, int], list[dict]], gz: bool) -> list[dict]:
    opener = gzip.open if gz else open
    matches: list[dict] = []
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if study == "GCST90011858":
                key = (row.get("chromosome", "").replace("chr", ""), int(float(row["base_pair_location"])))
                ea, oa = row.get("effect_allele", ""), row.get("other_allele", "")
                beta, se, eaf, pval = row.get("beta", ""), row.get("standard_error", ""), row.get("effect_allele_frequency", ""), row.get("p_value", "")
            else:
                variant = row.get("SNP", "")
                parts = variant.split("_")
                if len(parts) < 5:
                    continue
                key = (parts[0].replace("chr", ""), int(parts[1]))
                ea, oa = row.get("EA", ""), row.get("NEA", "")
                beta, se, eaf, pval = row.get("BETA", ""), row.get("SE", ""), "", row.get("P", "")
            for eqtl in by_coord.get(key, []):
                matches.append({"study": study, "gene": eqtl["gene"], "tissue": eqtl["tissue"], "eqtl_variant_id": eqtl["variant_id"], "gwas_variant": row.get("variant_id", row.get("SNP", "")), "chrom37": key[0], "pos37": key[1], "eqtl_ref38": eqtl.get("ref38", ""), "eqtl_alt38": eqtl.get("alt38", ""), "gwas_effect_allele": ea, "gwas_other_allele": oa, "gwas_beta": beta, "gwas_se": se, "gwas_eaf": eaf, "gwas_p": pval, "eqtl_slope": eqtl.get("slope", ""), "eqtl_slope_se": eqtl.get("slope_se", ""), "eqtl_p": eqtl.get("pval_nominal", ""), "eqtl_maf": eqtl.get("maf", "")})
    return matches


def main() -> None:
    rows = read_eqtl()
    by_coord, mappings = index_eqtl(rows)
    (OUT / "GTEx_v8_candidate_b38_to_b37_coordinate_map.json").write_text(json.dumps(mappings, indent=2), encoding="utf-8")
    results = []
    gwas_files = [("GCST90011858", GWAS / "GCST90011858/GCST90011858_buildGRCh37.tsv", False), ("GCST1010616", GWAS / "GCST1010616/Japan_PC.tsv.gz", True)]
    for study, path, gz in gwas_files:
        if path.exists():
            results.extend(scan_gwas(path, study, by_coord, gz))
    fields = ["study", "gene", "tissue", "eqtl_variant_id", "gwas_variant", "chrom37", "pos37", "eqtl_ref38", "eqtl_alt38", "gwas_effect_allele", "gwas_other_allele", "gwas_beta", "gwas_se", "gwas_eaf", "gwas_p", "eqtl_slope", "eqtl_slope_se", "eqtl_p", "eqtl_maf"]
    with (OUT / "GTEx_PD​​AC_GWAS_coordinate_overlap.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    # Write a clean ASCII filename as well; the first path is retained only if a filesystem previously created it.
    with (OUT / "GTEx_PDAC_GWAS_coordinate_overlap.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    counts = {(study, gene): sum(1 for r in results if r["study"] == study and r["gene"] == gene) for study, _, _ in gwas_files for gene in sorted(TARGETS)}
    report = ["# GTEx eQTL–PDAC GWAS coordinate overlap", "", "GTEx GRCh38 coordinates were mapped to GRCh37 using Ensembl REST. PDAC GWAS files were scanned by chromosome and position; allele orientation remains to be harmonised after reference/alternate allele liftover.", ""]
    for (study, gene), count in counts.items():
        report.append(f"- {study} / {gene}: {count} coordinate overlaps.")
    report.extend(["", "The GCST90011858 file provides beta, SE and EAF. Japan_PC provides beta and SE but its EAF is not in the main table and must be recovered or treated as unavailable for palindromic-SNP checks.", "Coordinate overlap is not evidence of a valid instrument until allele matching, strand resolution, LD clumping and F-statistic filtering are completed."])
    (OUT / "GTEx_PDAC_GWAS_coordinate_overlap_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
