"""airway_scaledcounts.subset.tsv에서 전체 샘플 발현량 합이 가장 높은 유전자 10개를 찾는다."""

import pandas as pd

DATA_PATH = "data/airway_scaledcounts.subset.tsv"
OUTPUT_PATH = "results/week3/top10_genes.tsv"

df = pd.read_csv(DATA_PATH, sep="\t", decimal=",", index_col="ensgene")

total_expression = df.sum(axis=1).sort_values(ascending=False)
top10 = total_expression.head(10)

print(top10)

top10.rename("total_counts").to_csv(OUTPUT_PATH, sep="\t")
