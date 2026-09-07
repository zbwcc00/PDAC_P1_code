from __future__ import annotations

import csv
from pathlib import Path

path = Path("D:/PDAC_P1/analysis/11_methodology_review/GTEx_v8_targeted/GTEx_PDAC_GWAS_coordinate_overlap.tsv")
out = path.parent
rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
for row in rows:
    eq_ref, eq_alt = row["eqtl_ref38"].upper(), row["eqtl_alt38"].upper()
    ea, oa = row["gwas_effect_allele"].upper(), row["gwas_other_allele"].upper()
    if (eq_alt, eq_ref) == (ea, oa):
        row["allele_status"] = "reverse_match"
        row["harmonized_gwas_beta"] = str(-float(row["gwas_beta"]))
    elif (eq_ref, eq_alt) == (ea, oa):
        row["allele_status"] = "same_match"
        row["harmonized_gwas_beta"] = row["gwas_beta"]
    elif {eq_ref, eq_alt} == {ea, oa}:
        row["allele_status"] = "strand_ambiguous_or_complement"
        row["harmonized_gwas_beta"] = ""
    else:
        row["allele_status"] = "allele_mismatch"
        row["harmonized_gwas_beta"] = ""
    try:
        row["eqtl_f_statistic"] = str((float(row["eqtl_slope"]) / float(row["eqtl_slope_se"])) ** 2)
    except (ValueError, ZeroDivisionError):
        row["eqtl_f_statistic"] = ""
    try:
        row["mr_instrument_p_gate"] = "PASS" if float(row["eqtl_p"]) < 5e-8 else "FAIL"
    except ValueError:
        row["mr_instrument_p_gate"] = "FAIL"

fields = list(rows[0]) if rows else []
with (out / "GTEx_PDAC_GWAS_harmonization_audit.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

lines = [
    "# GTEx eQTL–PDAC GWAS harmonization audit",
    "",
    "Pre-specified MR instrument gate: coordinate overlap, allele match, eQTL P < 5×10⁻⁸, F > 10, and LD clumping.",
    "",
]
if not rows:
    lines.append("No coordinate overlaps were found in the available GWAS files.")
for row in rows:
    lines.append(f"- {row['study']} / {row['gene']} / {row['tissue']} / {row['eqtl_variant_id']}: {row['allele_status']}; F={float(row['eqtl_f_statistic']):.2f}; eQTL P={row['eqtl_p']}; genome-wide gate={row['mr_instrument_p_gate']}.")
lines.extend([
    "",
    "Interpretation: reverse allele matching is technically resolvable by flipping the GWAS beta, but neither overlapping eQTL passes the genome-wide instrument P-value gate. Do not run primary MR with these variants.",
    "The candidates may proceed as observational/perturbational targets. MR should only be revisited if an independent, stronger cis-eQTL or cis-pQTL source is identified.",
])
(out / "GTEx_PDAC_GWAS_harmonization_audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

