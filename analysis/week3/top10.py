from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[2]
path = root / "data" / "airway_scaledcounts.subset.tsv"
counts = pd.read_csv(path, sep="\t", decimal=",", index_col="ensgene")
total = counts.sum(axis=1)
top10 = total.sort_values(ascending=False).head(10)
print(top10)
