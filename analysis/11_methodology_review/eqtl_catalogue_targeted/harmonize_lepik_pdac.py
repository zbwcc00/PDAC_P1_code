import csv
import gzip
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests


BASE = Path("D:/PDAC_P1")
OUT = BASE / "analysis/11_methodology_review/eqtl_catalogue_targeted"
EQTL = OUT / "LEPIK_targeted_cis_eqtl.tsv"
OUT.mkdir(parents=True, exist_ok=True)


def load_variants():
    rows = []
    seen = set()
    with EQTL.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if float(row["pvalue"]) < 5e-8 and row["rsid"] != "NA" and row["variant"] not in seen:
                rows.append(row)
                seen.add(row["variant"])
    return rows


def map_variant(row):
    chrom, pos, ref, alt = row["variant"].split("_")[:4]
    url = f"https://rest.ensembl.org/map/human/GRCh38/{chrom[3:]}:{pos}..{pos}/GRCh37"
    for attempt in range(6):
        try:
            response = requests.get(url, headers={"Content-Type": "application/json", "Accept": "application/json"}, timeout=45)
            if response.status_code == 503:
                time.sleep(3 + attempt * 2)
                continue
            response.raise_for_status()
            mappings = response.json().get("mappings", [])
            if not mappings:
                return None
            mapped = mappings[0]["mapped"]
            return {**row, "chrom37": str(mapped["seq_region_name"]), "pos37": str(mapped["start"]), "eqtl_ref38": ref, "eqtl_alt38": alt}
        except Exception:
            if attempt == 5:
                return None
            time.sleep(3 + attempt * 2)
    return None


def scan_gwas(path, study, mapped, gz):
    by_coord = {(row["chrom37"], int(row["pos37"])): row for row in mapped}
    opener = gzip.open if gz else open
    matches = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if study == "GCST90011858":
                key = (row["chromosome"].replace("chr", ""), int(float(row["base_pair_location"])))
                fields = (row.get("variant_id", ""), row.get("effect_allele", ""), row.get("other_allele", ""), row.get("beta", ""), row.get("standard_error", ""), row.get("effect_allele_frequency", ""), row.get("p_value", ""))
            else:
                parts = row.get("SNP", "").split("_")
                if len(parts) < 5:
                    continue
                key = (parts[0].replace("chr", ""), int(parts[1]))
                fields = (row.get("SNP", ""), row.get("EA", ""), row.get("NEA", ""), row.get("BETA", ""), row.get("SE", ""), "", row.get("P", ""))
            if key not in by_coord:
                continue
            eqtl = by_coord[key]
            matches.append({"study": study, "gene": eqtl["target_gene"], "eqtl_variant": eqtl["variant"], "rsid": eqtl["rsid"], "chrom37": key[0], "pos37": key[1], "eqtl_p": eqtl["pvalue"], "eqtl_beta": eqtl["beta"], "eqtl_se": eqtl["se"], "eqtl_F": eqtl["F"], "eqtl_ref38": eqtl["eqtl_ref38"], "eqtl_alt38": eqtl["eqtl_alt38"], "gwas_variant": fields[0], "gwas_ea": fields[1], "gwas_oa": fields[2], "gwas_beta": fields[3], "gwas_se": fields[4], "gwas_eaf": fields[5], "gwas_p": fields[6]})
    return matches


def main():
    variants = load_variants()
    mapped = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(map_variant, row) for row in variants]
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                mapped.append(result)
    (OUT / "LEPIK_genomewide_b38_to_b37_map.json").write_text(json.dumps(mapped, indent=2), encoding="utf-8")
    gwas = [("GCST90011858", BASE / "data/04_genetics/gwas/GCST90011858/GCST90011858_buildGRCh37.tsv", False), ("GCST1010616", BASE / "data/04_genetics/gwas/GCST1010616/Japan_PC.tsv.gz", True)]
    matches = []
    for study, path, gz in gwas:
        matches.extend(scan_gwas(path, study, mapped, gz))
    fields = ["study", "gene", "eqtl_variant", "rsid", "chrom37", "pos37", "eqtl_p", "eqtl_beta", "eqtl_se", "eqtl_F", "eqtl_ref38", "eqtl_alt38", "gwas_variant", "gwas_ea", "gwas_oa", "gwas_beta", "gwas_se", "gwas_eaf", "gwas_p"]
    with (OUT / "LEPIK_PDAC_overlap.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(matches)
    print(f"eQTL genome-wide variants: {len(variants)}")
    print(f"Mapped to GRCh37: {len(mapped)}")
    print(f"PDAC coordinate overlaps: {len(matches)}")
    for row in matches:
        print(row)


if __name__ == "__main__":
    main()
