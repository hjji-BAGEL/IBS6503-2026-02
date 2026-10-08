# /// script
# dependencies = ["pydeseq2==0.5.2", "pandas<2.2"]
# ///

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

REPO = Path("/mnt/c/Users/현종/OneDrive - 인하대학교/문서/Qoka/IBS6503-2026-02")
out_dir = Path(".")  # 이 run의 결과 폴더 (쓰기 가능). 나중에 results/week5/로 옮김

GENES = {"FKBP5": "ENSG00000096060", "TSC22D3": "ENSG00000157514", "PER1": "ENSG00000179094",
         "ZBTB16": "ENSG00000109906", "TARDBP": "ENSG00000120948", "ANGPTL7": "ENSG00000171819"}

counts = pd.read_csv(REPO / "data" / "airway_scaledcounts.subset.tsv", sep="\t", decimal=",", index_col="ensgene")
counts_df = counts.T.astype(int)

metadata = pd.read_csv(REPO / "data" / "week4" / "samples.tsv", sep="\t", index_col="sample")
metadata = metadata.loc[counts_df.index, ["cell", "dex"]]
metadata["dex"] = pd.Categorical(metadata["dex"], categories=["untrt", "trt"])

counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~dex", quiet=True)
    dds.deseq2()

ds = DeseqStats(dds, contrast=["dex", "trt", "untrt"], quiet=True)
ds.summary()
res = ds.results_df.copy()

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ds.lfc_shrink(coeff="dex[T.trt]")
res_shrunk = ds.results_df.copy()

res_shrunk.sort_values("padj").to_csv(out_dir / "deseq2_dex_results.csv")

sig = res["padj"] < 0.05
sig_shrunk = res_shrunk["padj"] < 0.05

print("=== 유전자 수 ===")
print("필터 후 유전자:", counts_df.shape[1])
print("p < 0.05:   ", (res["pvalue"] < 0.05).sum())
print("padj < 0.05:", sig.sum(), "(증가", (sig & (res.log2FoldChange > 0)).sum(),
      "/ 감소", (sig & (res.log2FoldChange < 0)).sum(), ")")
print("padj NaN (검정 제외):", res["padj"].isna().sum())

print("\n=== Positive control 유전자 (glucocorticoid 반응) ===")
names = ["FKBP5", "TSC22D3", "PER1", "ZBTB16", "TARDBP"]
table = res_shrunk.loc[[GENES[n] for n in names], ["baseMean", "log2FoldChange", "padj"]].set_axis(names)
print(table.round(3))

print("\n=== padj 상위 10개 유전자 ===")
top10 = res_shrunk.sort_values("padj")[["baseMean", "log2FoldChange", "padj"]].head(10)
print(top10.round(3))

fig, ax = plt.subplots(figsize=(5, 4))
up = sig_shrunk & (res_shrunk.log2FoldChange > 0)
down = sig_shrunk & (res_shrunk.log2FoldChange < 0)
ax.scatter(np.log10(res_shrunk.baseMean.clip(lower=1)), res_shrunk.log2FoldChange, s=4, color="lightgray")
ax.scatter(np.log10(res_shrunk.baseMean[up].clip(lower=1)), res_shrunk.log2FoldChange[up], s=5, color="tab:red")
ax.scatter(np.log10(res_shrunk.baseMean[down].clip(lower=1)), res_shrunk.log2FoldChange[down], s=5, color="tab:blue")
ax.axhline(0, color="black", linewidth=0.8)
ax.set_xlabel("log10(baseMean)")
ax.set_ylabel("log2FoldChange (trt / untrt, shrunken)")
ax.set_title("MA plot: dex (trt vs untrt)")
fig.tight_layout()
fig.savefig(out_dir / "ma_plot.png", dpi=150)

print("\n저장: deseq2_dex_results.csv, ma_plot.png (run 폴더 안)")
