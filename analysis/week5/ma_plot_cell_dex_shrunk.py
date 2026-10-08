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
out_dir = Path(".")

counts = pd.read_csv(REPO / "data" / "airway_scaledcounts.subset.tsv", sep="\t", decimal=",", index_col="ensgene")
counts_df = counts.T.astype(int)

metadata = pd.read_csv(REPO / "data" / "week4" / "samples.tsv", sep="\t", index_col="sample")
metadata = metadata.loc[counts_df.index, ["cell", "dex"]]
metadata["dex"] = pd.Categorical(metadata["dex"], categories=["untrt", "trt"])

counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    dds_cell = DeseqDataSet(counts=counts_df, metadata=metadata, design="~cell + dex", quiet=True)
    dds_cell.deseq2()
    ds_cell = DeseqStats(dds_cell, contrast=["dex", "trt", "untrt"], quiet=True)
    ds_cell.summary()
    ds_cell.lfc_shrink(coeff="dex[T.trt]")

res_shrunk = ds_cell.results_df.copy()

# pydeseq2 내장 MA plot (ds_cell.plot_MA)
fig = ds_cell.plot_MA(s=4, save_path=str(out_dir / "ma_plot_cell_dex_shrunk.png"))

sig = res_shrunk["padj"] < 0.05
print("padj < 0.05 (shrink 후):", sig.sum(),
      "(증가", (sig & (res_shrunk.log2FoldChange > 0)).sum(),
      "/ 감소", (sig & (res_shrunk.log2FoldChange < 0)).sum(), ")")
print("저장: ma_plot_cell_dex_shrunk.png (run 폴더 안)")
