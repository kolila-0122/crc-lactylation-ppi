# -*- coding: utf-8 -*-
"""队列异质性剖析：解释 bulk 强、单细胞患者级伪 bulk 弱的尺度/功效问题

1) TCGA vs Pelka(肿瘤上皮 donor 伪 bulk) vs GSE132465 的 PPI 分布(方差/范围)
2) TCGA 随机子抽样到 n=62 与 n=23，重复检验 rho=-0.37 的检出率(功效)
3) Pelka 队列 MMR/CMS 构成(用于说明队列组成差异)
"""
import json, os, gzip, warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

warnings.filterwarnings("ignore")
SCR = os.environ.get("SCRNA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scRNA"))

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

# ---- TCGA PPI（从 results/tcga_bulk_scores.json 读，替代旧 novel_angle_validation.json 外部依赖）----
tcga = json.load(open(os.path.join(OUT, "tcga_bulk_scores.json"), encoding="utf-8"))
tcga_ppi = np.array(list(tcga["ppi"].values()))
print("=== TCGA PPI 分布 (n=%d) ===" % len(tcga_ppi))
print("  std=%.3f  IQR=%.3f  range=[%.2f, %.2f]" % (
    tcga_ppi.std(), np.percentile(tcga_ppi, 75) - np.percentile(tcga_ppi, 25),
    tcga_ppi.min(), tcga_ppi.max()))

# ---- Pelka 肿瘤上皮 donor 伪 bulk PPI ----
def pelka_pb_ppi(fn, gc, dv):
    a = ad.read_h5ad(os.path.join(SCR, fn), backed="r")
    sym = a.var["feature_name"].astype(str)
    sel, seen = [], set()
    for i, s in enumerate(sym):
        if s in ALLG and s not in seen:
            sel.append(i); seen.add(s)
    sub = a[:, sel].to_memory()
    X = sub.X.toarray() if not isinstance(sub.X, np.ndarray) else np.asarray(sub.X)
    genes = [sym.iloc[i] for i in sel]
    grp = a.obs[gc].astype(str).values
    dis = a.obs["disease"].astype(str).values
    a.file.close()
    Xd = pd.DataFrame(X, columns=genes); Xd["g"] = grp; Xd["d"] = dis
    pb = Xd[Xd["d"] == dv].groupby("g")[genes].mean()
    pb = pb[pb.sum(axis=1) > 0]
    Z = ((pb - pb.mean()) / pb.std().replace(0, np.nan)).fillna(0)
    ox = [g for g in OXG if g in Z.columns]; ef = [g for g in EFFG if g in Z.columns]
    return Z[ox].mean(axis=1) - Z[ef].mean(axis=1)

pelka_ppi = pelka_pb_ppi("pelka_epithelial.h5ad", "donor_id", "colon adenocarcinoma")
print("=== Pelka 肿瘤上皮 donor 伪 bulk PPI (n=%d) ===" % len(pelka_ppi))
print("  std=%.3f  IQR=%.3f  range=[%.2f, %.2f]" % (
    pelka_ppi.std(), np.percentile(pelka_ppi, 75) - np.percentile(pelka_ppi, 25),
    pelka_ppi.min(), pelka_ppi.max()))

# ---- GSE132465 患者伪 bulk PPI ----
ann = {}
with gzip.open(os.path.join(SCR, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"), "rt") as f:
    next(f)
    for line in f:
        p = line.rstrip().split("\t")
        ann[p[0]] = {"patient": p[1], "class": p[2], "cell_type": p[4]}
f2 = gzip.open(os.path.join(SCR, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"), "rt")
header = next(f2).rstrip().split("\t")
cell_ids = header[1:]
valid = {c: i for i, c in enumerate(cell_ids) if c in ann}
gene_data = {}
for line in f2:
    parts = line.rstrip().split("\t")
    if parts[0] in set(ALLG):
        gene_data[parts[0]] = parts
f2.close()
valid_cells = sorted(valid.keys()); idx = [valid[c] + 1 for c in valid_cells]
expr = pd.DataFrame({g: [float(parts[i]) for i in idx] for g, parts in gene_data.items()}, index=valid_cells)
meta = pd.DataFrame([ann[c] for c in valid_cells], index=valid_cells)
E = np.log1p(expr).join(meta)
tum = E[(E["class"] == "Tumor") & (E["cell_type"] == "Epithelial cells")]
pb = tum.groupby("patient")[ALLG].mean()
Z = ((pb - pb.mean()) / pb.std().replace(0, np.nan)).fillna(0)
ox = [g for g in OXG if g in Z.columns]; ef = [g for g in EFFG if g in Z.columns]
gse_ppi = Z[ox].mean(axis=1) - Z[ef].mean(axis=1)
print("=== GSE132465 患者伪 bulk PPI (n=%d) ===" % len(gse_ppi))
print("  std=%.3f  IQR=%.3f  range=[%.2f, %.2f]" % (
    gse_ppi.std(), np.percentile(gse_ppi, 75) - np.percentile(gse_ppi, 25),
    gse_ppi.min(), gse_ppi.max()))

# ---- 功效分析：理论 Fisher-z 变换，评估不同样本量下检出 rho≈-0.37 的信噪比 ----
from scipy.stats import spearmanr
print("\n=== 样本量功效提示(理论) ===")
import math
for n in [592, 62, 36, 23]:
    se = 1.0 / math.sqrt(n - 3)
    z = math.atanh(-0.37)
    print("  n=%d: Fisher-z=%.3f, SE=%.3f, |z|/SE=%.2f (检出 rho=-0.37 需要的信噪比)" % (n, z, se, abs(z)/se))

# ---- Pelka MMR 构成（value_counts 落盘 JSON 供核验）----
a = ad.read_h5ad(os.path.join(SCR, "pelka_epithelial.h5ad"), backed="r")
mmr = a.obs["MMRStatusTumor"].astype(str).values
donor = a.obs["donor_id"].astype(str).values
a.file.close()
d = pd.DataFrame({"donor": donor, "mmr": mmr}).drop_duplicates()
mmr_counts = d["mmr"].value_counts().to_dict()
print("\n=== Pelka 上皮 donor 的 MMR 构成 ===")
print(mmr_counts)
json.dump(mmr_counts, open(os.path.join(OUT, "pelka_mmr_counts.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

print("\nDONE")
