from __future__ import annotations

import csv
import json
import subprocess
import sys
from urllib.parse import quote
from pathlib import Path


BASE = Path(r"D:\第二篇大论文")
DRUG = BASE / "analysis" / "13_virtual_perturbation" / "drugreflector"
MAP = DRUG / "drugreflector_compound_mapping.tsv"
REST = Path(r"C:\Users\HUAWEI\.codex\plugins\cache\openai-api-curated\life-science-research\1e285826\skills\pubchem-pug-skill\scripts\rest_request.py")


def call_pubchem(identifier: str, mode: str):
    encoded = quote(identifier, safe="")
    request = {
        "base_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        "path": f"compound/{mode}/{encoded}/property/IUPACName,MolecularFormula,CanonicalSMILES,IsomericSMILES/JSON",
        "record_path": "PropertyTable.Properties",
        "max_items": 3,
        "timeout_sec": 20,
    }
    result = subprocess.run(
        [sys.executable, str(REST)],
        input=json.dumps(request),
        text=True,
        capture_output=True,
        check=False,
    )
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        return {}
    records = payload.get("records", [])
    return records[0] if payload.get("ok") and records else {}


def main():
    with MAP.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = list(rows[0]) + ["pubchem_cid", "pubchem_iupac_name", "pubchem_formula", "pubchem_status"]
    cache = {}
    for index, row in enumerate(rows, start=1):
        name = row.get("cmap_name", "")
        key = row.get("inchi_key", "")
        query_mode = "name" if name and not name.startswith("BRD-") else "inchikey"
        query = name if query_mode == "name" else key
        cache_key = f"{query_mode}:{query}"
        if query and cache_key not in cache:
            cache[cache_key] = call_pubchem(query, query_mode)
        record = cache.get(cache_key, {})
        row["pubchem_cid"] = str(record.get("CID", ""))
        row["pubchem_iupac_name"] = record.get("IUPACName", "")
        row["pubchem_formula"] = record.get("MolecularFormula", "")
        row["pubchem_status"] = "matched" if record else ("no_identifier" if not query else "not_found")
        print(f"{index}/{len(rows)} {row['compound']} {row['pubchem_status']}", flush=True)
    out = DRUG / "drugreflector_compound_mapping_pubchem.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
