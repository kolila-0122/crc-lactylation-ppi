# -*- coding: utf-8 -*-
"""统一出图：从 results/ 与数据文件读，产出 Fig1–Fig5（输出名与包内图名一致）"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def _root_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    # 同级有 results/figures → 根=同级（开发布局）；否则根=上级（包内 code/ 布局）
    if os.path.isdir(os.path.join(here, "results")) or os.path.isdir(os.path.join(here, "figures")):
        return here
    return os.path.dirname(here)

ROOT = _root_dir()
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "figure.dpi": 300, "savefig.bbox": "tight",
                     "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})
POS = "#2b6cb0"; NEG = "#c53030"
def col(v): return POS if v > 0 else NEG

def save_fig(name):
    """输出 PNG（预览）+ TIFF（投稿用，300 dpi，LZW 压缩）"""
    plt.savefig(os.path.join(FIG, name + ".png"))
    plt.savefig(os.path.join(FIG, name + ".tif"), dpi=300, format="tiff",
                pil_kwargs={"compression": "tiff_lzw"})

def load(p):
    """先找 results/（包内布局），再找根目录（开发布局）"""
    for d in (RES, ROOT):
        f = os.path.join(d, p)
        if os.path.exists(f):
            return json.load(open(f, encoding="utf-8"))
    raise FileNotFoundError(p)

def load_csv(p):
    for d in (RES, ROOT):
        f = os.path.join(d, p)
        if os.path.exists(f):
            return pd.read_csv(f, index_col=0)
    raise FileNotFoundError(p)

tcga = load("tcga_bulk_scores.json")
sc = load("sc_single_arm_correlations.json")
gse = load("gse39582_WE_validation.json")

# ===== Fig1: TCGA bulk discovery =====
labels1 = ["Writer", "Eraser", "WE difference", "Target genes"]
rhos1 = [tcga["writer_rho"], tcga["eraser_rho"], tcga["WE_rho"], tcga["target_rho"]]
ps1 = [tcga["writer_p"], tcga["eraser_p"], tcga["WE_p"], tcga["target_p"]]
fig, ax = plt.subplots(figsize=(7.6, 3.6))
bars = ax.barh(labels1, rhos1, color=[col(v) for v in rhos1], alpha=0.88)
ax.axvline(0, color="black", lw=0.8)
for b, r, p in zip(bars, rhos1, ps1):
    ax.text(r + (0.012 if r >= 0 else -0.012), b.get_y() + b.get_height()/2,
            "rho=%+.3f\np=%.1e" % (r, p), va="center", ha="left" if r >= 0 else "right", fontsize=8)
ax.set_xlim(-0.46, 0.50)
ax.set_xlabel("Spearman rho (PPI vs score)")
ax.set_title("PPI vs lactylation-related scores (TCGA COADREAD, n=592)", fontsize=10)
save_fig("Fig1_bulk_discovery")
plt.close()

# ===== Fig2: cross-cohort =====
cohorts = ["TCGA\n(n=592)", "GSE39582\n(n=585)", "Pelka tumour\n(n=62)", "GSE132465\n(n=23)"]
writer = [tcga["writer_rho"], gse["writer"]["rho"], sc["pelka_tumour_epithelium"]["writer_rho"], sc["gse132465_tumour_epithelium"]["writer_rho"]]
eraser = [tcga["eraser_rho"], gse["eraser"]["rho"], sc["pelka_tumour_epithelium"]["eraser_rho"], sc["gse132465_tumour_epithelium"]["eraser_rho"]]
we = [tcga["WE_rho"], gse["WE"]["rho"], sc["pelka_tumour_epithelium"]["WE_rho"], sc["gse132465_tumour_epithelium"]["WE_rho"]]
x = np.arange(len(cohorts)); w = 0.26
fig, ax = plt.subplots(figsize=(8.4, 4.6))
ax.bar(x - w, writer, w, color=[col(v) for v in writer], alpha=0.85, label="Writer")
ax.bar(x, eraser, w, color=[col(v) for v in eraser], alpha=0.85, label="Eraser")
ax.bar(x + w, we, w, color=[col(v) for v in we], alpha=0.6, label="WE difference")
ax.axhline(0, color="black", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(cohorts, fontsize=9)
ax.set_ylabel("Spearman rho (PPI vs score)")
ax.set_title("PPI correlation with lactylation writer and eraser scores across cohorts", fontsize=11)
ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.set_ylim(-0.5, 0.55)
save_fig("Fig2_cross_cohort_replication")
plt.close()

# ===== Fig3: cell-type =====
pb = load("sc_lactylation_pseudobulk_results.json")
ct = pb["pelka_celltype_pseudobulk"]
cts = list(ct["PPI"].keys())
ppi_v = [ct["PPI"][c] for c in cts]
we_v = [ct["WE_ratio"][c] for c in cts]
fig, ax = plt.subplots(figsize=(6.4, 5.2))
ax.axhline(0, color="grey", lw=0.6, ls="--"); ax.axvline(0, color="grey", lw=0.6, ls="--")
for c, xx, yy in zip(cts, ppi_v, we_v):
    ax.scatter(xx, yy, s=90, color=POS if xx > 0 else NEG, zorder=3)
    ax.annotate(c, (xx, yy), textcoords="offset points", xytext=(6, 4), fontsize=9)
ax.set_xlabel("PPI (cell-type pseudobulk, z)")
ax.set_ylabel("WE difference (z)")
ax.set_title("Cell-type pyruvate partitioning vs lactylation writer-minus-eraser (Pelka, 7 cell types)", fontsize=10)
save_fig("Fig3_celltype_WE")
plt.close()

# ===== Fig4: spatial =====
sp = load_csv("spatial_lactylation_results.csv")
sp = sp.sort_values("PPI_WE_rho")
fig, ax = plt.subplots(figsize=(8.0, 3.6))
ax.bar(range(len(sp)), sp["PPI_WE_rho"], color=[col(v) for v in sp["PPI_WE_rho"]], alpha=0.85)
ax.axhline(0, color="black", lw=0.8)
ax.axhline(sp["PPI_WE_rho"].mean(), color="grey", lw=1, ls="--",
           label="mean rho = %.3f" % sp["PPI_WE_rho"].mean())
for i, (idx, row) in enumerate(sp.iterrows()):
    ax.text(i, row["PPI_WE_rho"] + (0.004 if row["PPI_WE_rho"] >= 0 else -0.004),
            "p=%.2g" % row["PPI_WE_p"], ha="center",
            va="bottom" if row["PPI_WE_rho"] >= 0 else "top", fontsize=6, rotation=90)
ax.set_xticks(range(len(sp))); ax.set_xticklabels(sp.index, rotation=60, fontsize=7)
ax.set_ylabel("Spot-level Spearman rho (PPI vs WE difference)")
ax.set_title("Spatial PPI-WE-difference coupling across 13 Visium slices", fontsize=10)
ax.legend(fontsize=8, frameon=False)
save_fig("Fig4_spatial_WE")
plt.close()

# ===== Fig5: single-gene =====
sg = tcga["single_gene"]
gse_sg = gse["single_gene"]
genes = ["EP300", "CREBBP", "KAT8", "KAT2A", "SIRT1", "SIRT3", "HDAC1", "HDAC2", "HDAC3"]
tcga_r = [sg[g]["rho"] for g in genes]
gse_r = [gse_sg[g]["rho"] for g in genes]
x = np.arange(len(genes)); w = 0.36
fig, ax = plt.subplots(figsize=(9.2, 4.6))
ax.bar(x - w/2, tcga_r, w, color=[col(v) for v in tcga_r], alpha=0.85, label="TCGA (n=592)")
ax.bar(x + w/2, gse_r, w, color=[col(v) for v in gse_r], alpha=0.5, label="GSE39582 (n=585)")
ax.axhline(0, color="black", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Spearman rho (PPI vs gene)")
ax.set_title("Per-gene PPI correlation with lactylation writer and eraser genes", fontsize=10.5)
ax.legend(fontsize=9, frameon=False, loc="upper right")
ax.axvspan(3.5, 8.5, color="#f0f4f8", alpha=0.5, zorder=0)
ax.text(6, ax.get_ylim()[1]*0.92, "erasers", ha="center", fontsize=8, color="#4a5568")
ax.text(1.5, ax.get_ylim()[1]*0.92, "writers", ha="center", fontsize=8, color="#4a5568")
save_fig("Fig5_single_gene")
plt.close()

print("ALL FIGS DONE:", sorted(os.listdir(FIG)))
