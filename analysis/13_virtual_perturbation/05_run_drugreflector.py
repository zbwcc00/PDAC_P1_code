from pathlib import Path
import hashlib
import json
import sys
import pandas as pd

root = Path(r"D:\PDAC_P1")
work = root / "analysis" / "13_virtual_perturbation" / "drugreflector"
repo = root / "software" / "drugreflector" / "drugreflector-main"
sys.path.insert(0, str(repo))
from drugreflector import DrugReflector
checkpoint_dir = repo / "checkpoints"
paths = [checkpoint_dir / f"model_fold_{i}.pt" for i in range(3)]
for path in paths:
    if not path.exists() or path.stat().st_size < 90_000_000:
        raise RuntimeError(f"Checkpoint incomplete: {path} ({path.stat().st_size if path.exists() else 0} bytes)")

long = pd.read_csv(work / "drugreflector_vscores.tsv", sep="\t")
wide = long.pivot_table(index="signature", columns="gene", values="vscore", aggfunc="median")
wide = wide.replace([float("inf"), float("-inf")], 0).fillna(0)
wide.columns = wide.columns.astype(str).str.upper()
wide = wide.loc[:, ~wide.columns.duplicated()]
wide = wide.sort_index()
reverse = -wide
reverse.index = [f"{x}__REVERSE" for x in reverse.index]
inputs = pd.concat([wide, reverse], axis=0)
inputs.to_csv(work / "drugreflector_input_vscores_matrix.tsv", sep="\t")

model = DrugReflector(checkpoint_paths=[str(p) for p in paths])
coverage = model.check_gene_coverage(inputs.columns)
(work / "drugreflector_landmark_coverage.json").write_text(json.dumps(coverage, indent=2), encoding="utf-8")
pred = model.predict(inputs, n_top=None)
pred.to_csv(work / "drugreflector_predictions_full.tsv", sep="\t")

rows = []
for signature in inputs.index:
    rank = pred[("rank", signature)].astype(float)
    logit = pred[("logit", signature)].astype(float)
    prob = pred[("prob", signature)].astype(float)
    tab = pd.DataFrame({"compound": rank.index, "rank": rank.values, "logit": logit.values, "prob": prob.values})
    tab["signature"] = signature
    tab["target"] = signature.split("__", 1)[0]
    tab["direction"] = "reverse" if signature.endswith("__REVERSE") else "state"
    rows.append(tab.sort_values(["rank", "prob"]).head(50))
top = pd.concat(rows, ignore_index=True)
top.to_csv(work / "drugreflector_top50.tsv", sep="\t", index=False)

checks = []
for path in paths:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    checks.append({"file": str(path), "bytes": path.stat().st_size, "md5": h.hexdigest()})
(work / "drugreflector_checkpoint_checksums.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
