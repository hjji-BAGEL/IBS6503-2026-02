# /// script
# dependencies = ["pydeseq2==0.5.2", "pandas<2.2"]
# ///

import warnings
from pathlib import Path

import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

REPO = Path("/mnt/c/Users/현종/OneDrive - 인하대학교/문서/Qoka/IBS6503-2026-02")
out_dir = Path(".")

GENES = {"FKBP5": "ENSG00000096060", "TSC22D3": "ENSG00000157514", "ANGPTL7": "ENSG00000171819"}

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
res_cell = ds_cell.results_df.copy()

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ds_cell.lfc_shrink(coeff="dex[T.trt]")
res_shrunk = ds_cell.results_df.copy()

res_shrunk.sort_values("padj").to_csv(out_dir / "deseq2_cell_dex_shrunk_results.csv")

names = ["FKBP5", "TSC22D3", "ANGPTL7"]
ids = [GENES[n] for n in names]
compare = pd.DataFrame({"shrink 전": res_cell.loc[ids, "log2FoldChange"].values,
                         "shrink 후": res_shrunk.loc[ids, "log2FoldChange"].values}, index=names)
print(compare.round(3))
print()
print("padj < 0.05 (shrink 후, padj는 shrink로 안 바뀜):", (res_shrunk.padj < 0.05).sum())
print("\n저장: deseq2_cell_dex_shrunk_results.csv (run 폴더 안)")
