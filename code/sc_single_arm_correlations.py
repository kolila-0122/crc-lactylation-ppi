# -*- coding: utf-8 -*-
"""单细胞伪 bulk 的 PPI~writer / PPI~eraser / PPI~WE 单臂相关（产出 Table 1 中 Pelka 与 GSE132465 行）

修复 B1：此前单臂系数仅以字面量形式出现在出图脚本，本脚本为权威计算来源，
产出与 Table 1 逐字一致的数字，并落盘 JSON 供正文引用。

数据：Pelka 上皮 h5ad（CELLxGENE，Pelka et al. Cell 2021）；GSE132465（GEO，Lee et al. Nat Genet 2020）。
两个数据文件均由用户按 data_availability.md 的 accession 下载到本机 DATA_DIR。
"""
import os, json, gzip, warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

warnings.filterwarnings("ignore")

# 数据目录：默认取脚本上级目录的 scRNA/，可用环境变量 SCRNA_DIR 覆盖
DATA_DIR = os.environ.get("SCRNA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scRNA"))

def _results_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.environ.get("RESULTS_DIR") or (
        os.path.join(here, "results") if os.path.isdir(os.path.join(here, "results"))
        else os.path.join(os.path.dirname(here), "results"))

OUT = _results_dir()
os.makedirs(OUT, exist_ok=True)

WRITER = ["EP300", "CREBBP", "KAT8", "KAT2A"]
ERASER = ["SIRT1", "SIRT3", "HDAC1", "HDAC2", "HDAC3"]
OXG = ["MPC1", "MPC2", "PDHA1", "PDHB", "DLAT", "DLD", "PDHX"]
EFFG = ["LDHA", "LDHB", "SLC16A1", "SLC16A3", "SLC16A7"]
ALLG = WRITER + ERASER + OXG + EFFG


def zscore(x):
    x = np.asarray(x, float)
    m, sd = np.nanmean(x), np.nanstd(x)
    return (x - m) / sd if sd > 0 else np.zeros_like(x)


def score_from_logexpr(E):
    Z = E.apply(zscore, axis=0)
    w = [g for g in WRITER if g in Z.columns]
    e = [g for g in ERASER if g in Z.columns]
    ox = [g for g in OXG if g in Z.columns]
    ef = [g for g in EFFG if g in Z.columns]
    out = pd.DataFrame(index=E.index)
    out["writer"] = Z[w].mean(axis=1)
    out["eraser"] = Z[e].mean(axis=1)
    out["WE"] = out["writer"] - out["eraser"]
    out["PPI"] = Z[ox].mean(axis=1) - Z[ef].mean(axis=1)
    return out


def arm_corr(S):
    r_w, p_w = stats.spearmanr(S["PPI"], S["writer"])
    r_e, p_e = stats.spearmanr(S["PPI"], S["eraser"])
    r_we, p_we = stats.spearmanr(S["PPI"], S["WE"])
    return {"writer_rho": float(r_w), "writer_p": float(p_w),
            "eraser_rho": float(r_e), "eraser_p": float(p_e),
            "WE_rho": float(r_we), "WE_p": float(p_we), "n": int(len(S))}


def pelka_pseudobulk():
    p = os.path.join(DATA_DIR, "pelka_epithelial.h5ad")
    a = ad.read_h5ad(p, backed="r")
    sym = a.var["feature_name"].astype(str)
    sel, seen = [], set()
    for i, s in enumerate(sym):
        if s in ALLG and s not in seen:
            sel.append(i); seen.add(s)
    sub = a[:, sel].to_memory()
    X = sub.X.toarray() if not isinstance(sub.X, np.ndarray) else np.asarray(sub.X)
    genes = [sym.iloc[i] for i in sel]
    donor = a.obs["donor_id"].astype(str).values
    disease = a.obs["disease"].astype(str).values
    a.file.close()
    Xd = pd.DataFrame(X, columns=genes)
    Xd["donor"] = donor; Xd["disease"] = disease
    return Xd


def gse132465_pseudobulk():
    # 注释文件列序（GEO 官方 6 列）：Index, Patient, Class, Sample, Cell_type, Cell_subtype
    ann = {}
    with gzip.open(os.path.join(DATA_DIR, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"), "rt") as f:
        next(f)
        for line in f:
            p = line.rstrip().split("\t")
            ann[p[0]] = {"patient": p[1], "class": p[2], "cell_type": p[4]}
    f2 = gzip.open(os.path.join(DATA_DIR, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"), "rt")
    header = next(f2).rstrip().split("\t")
    cell_ids = header[1:]
    valid = {c: i for i, c in enumerate(cell_ids) if c in ann}
    gene_data = {}
    for line in f2:
        parts = line.rstrip().split("\t")
        if parts[0] in set(ALLG):
            gene_data[parts[0]] = parts
    f2.close()
    valid_cells = sorted(valid.keys())
    idx = [valid[c] + 1 for c in valid_cells]
    expr = pd.DataFrame({g: [float(parts[i]) for i in idx] for g, parts in gene_data.items()},
                        index=valid_cells)
    meta = pd.DataFrame([ann[c] for c in valid_cells], index=valid_cells)
    E = np.log1p(expr).join(meta)
    return E


results = {}

# Pelka 肿瘤上皮 donor 伪 bulk
pelka = pelka_pseudobulk()
tum = pelka[pelka["disease"] == "colon adenocarcinoma"]
pb = tum.groupby("donor")[ALLG].mean()
pb = pb[pb.sum(axis=1) > 0]
results["pelka_tumour_epithelium"] = arm_corr(score_from_logexpr(pb))
print("Pelka tumour:", results["pelka_tumour_epithelium"])

# GSE132465 肿瘤上皮 patient 伪 bulk
gse = gse132465_pseudobulk()
tum_epi = gse[(gse["class"] == "Tumor") & (gse["cell_type"] == "Epithelial cells")]
pb2 = tum_epi.groupby("patient")[ALLG].mean()
results["gse132465_tumour_epithelium"] = arm_corr(score_from_logexpr(pb2))
print("GSE132465 tumour:", results["gse132465_tumour_epithelium"])

with open(os.path.join(OUT, "sc_single_arm_correlations.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, indent=1, ensure_ascii=False)
print("DONE -> results/sc_single_arm_correlations.json")
