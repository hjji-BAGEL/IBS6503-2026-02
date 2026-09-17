"""airway_scaledcounts.subset.tsv에서 dexamethasone 처리(treated)로 발현량이
가장 많이 증가한 유전자 10개를 찾는다. 증가량은 (treated 평균 - control 평균)으로 계산."""

import pandas as pd

DATA_PATH = "data/airway_scaledcounts.subset.tsv"
OUTPUT_PATH = "results/week3/up10_genes.tsv"

df = pd.read_csv(DATA_PATH, sep="\t", decimal=",", index_col="ensgene")

control_cols = [c for c in df.columns if c.endswith("_control")]
treated_cols = [c for c in df.columns if c.endswith("_treated")]

increase = df[treated_cols].mean(axis=1) - df[control_cols].mean(axis=1)
top10 = increase.sort_values(ascending=False).head(10)

print(top10)

top10.rename("mean_increase").to_csv(OUTPUT_PATH, sep="\t")
