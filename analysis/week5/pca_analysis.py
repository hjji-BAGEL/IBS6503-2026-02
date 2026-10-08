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

REPO = Path("/mnt/c/Users/현종/OneDrive - 인하대학교/문서/Qoka/IBS6503-2026-02")
out_dir = Path(".")

counts = pd.read_csv(REPO / "data" / "airway_scaledcounts.subset.tsv", sep="\t", decimal=",", index_col="ensgene")
counts_df = counts.T.astype(int)

metadata = pd.read_csv(REPO / "data" / "week4" / "samples.tsv", sep="\t", index_col="sample")
metadata = metadata.loc[counts_df.index, ["cell", "dex"]]
metadata["dex"] = pd.Categorical(metadata["dex"], categories=["untrt", "trt"])

counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~dex", quiet=True)
    dds.deseq2()

# VST (rnaseqGene의 vst(dds, blind = FALSE)와 동일: 설계를 고려한 분산 추정 사용)
dds.vst(use_design=True)
vst = pd.DataFrame(dds.layers["vst_counts"], index=counts_df.index, columns=counts_df.columns)

# 샘플 간 분산이 큰 유전자 500개 선택 -> 유전자별 평균을 빼고 PCA(SVD)
top = vst.var().sort_values(ascending=False).index[:500]
X = vst[top] - vst[top].mean()
U, S, _ = np.linalg.svd(X.values, full_matrices=False)
pcs = pd.DataFrame(U[:, :2] * S[:2], index=vst.index, columns=["PC1", "PC2"])
ratio = S**2 / (S**2).sum()

fig, ax = plt.subplots(figsize=(5, 4))
for s, row in pcs.iterrows():
    color = "tab:blue" if metadata.loc[s, "dex"] == "trt" else "gray"
    marker = "o" if metadata.loc[s, "cell"] == "N61311" else "s"
    ax.scatter(row.PC1, row.PC2, s=90, color=color, marker=marker)
    ax.annotate(s, (row.PC1, row.PC2), xytext=(6, 6), textcoords="offset points")
ax.set_xlabel(f"PC1 ({ratio[0]:.0%})")
ax.set_ylabel(f"PC2 ({ratio[1]:.0%})")
ax.set_title("PCA (top 500 variable genes, VST)")
fig.tight_layout()
fig.savefig(out_dir / "pca_plot.png", dpi=150)

print("PC1 설명 비율:", f"{ratio[0]:.1%}")
print("PC2 설명 비율:", f"{ratio[1]:.1%}")
print(pcs.round(2))
print("\n저장: pca_plot.png (run 폴더 안)")
