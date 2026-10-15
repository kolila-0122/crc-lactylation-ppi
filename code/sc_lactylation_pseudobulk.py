# -*- coding: utf-8 -*-
"""子课题② 单细胞层乳酸化 writer/eraser —— 伪 bulk 版（修正全局 z 塌缩）

方法：细胞聚合到 donor(患者)/细胞类型级伪 bulk → 跨样本 z-score → PPI 与 WE difference。
本脚本产出细胞类型伪 bulk 的 PPI/WE 结果（供 Fig3）；患者级单臂相关见 sc_single_arm_correlations.py。
"""
import os, json, gzip, warnings
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

def score_from_logexpr(E):
    """E: DataFrame (samples x genes), 已经是 log 表达。返回 PPI/writer/eraser/WE_ratio (跨样本 z)"""
    Z = ((E - E.mean()) / E.std().replace(0, np.nan)).fillna(0)
    w = [g for g in WRITER if g in Z.columns]
    e = [g for g in ERASER if g in Z.columns]
    ox = [g for g in OXG if g in Z.columns]
    ef = [g for g in EFFG if g in Z.columns]
    S = pd.DataFrame(index=E.index)
    S["writer"] = Z[w].mean(axis=1)
    S["eraser"] = Z[e].mean(axis=1)
    S["WE_ratio"] = S["writer"] - S["eraser"]
    S["PPI"] = Z[ox].mean(axis=1) - Z[ef].mean(axis=1)
    return S

results = {}

# ===================== Pelka 上皮：患者级伪 bulk =====================
a = ad.read_h5ad(os.path.join(SCR, "pelka_epithelial.h5ad"), backed="r")
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
mmr = a.obs["MMRStatusTumor"].astype(str).values if "MMRStatusTumor" in a.obs.columns else None
a.file.close()

Xd = pd.DataFrame(X, columns=genes)
Xd["donor"] = donor
Xd["disease"] = disease
if mmr is not None:
    Xd["MMR"] = mmr

# 肿瘤上皮 donor 伪 bulk
tum = Xd[Xd["disease"] == "colon adenocarcinoma"]
donor_pb = tum.groupby("donor")[genes].mean()
donor_pb = donor_pb[donor_pb.sum(axis=1) > 0]  # 去掉空 donor
S_donor = score_from_logexpr(donor_pb)
r_donor, p_donor = stats.spearmanr(S_donor["PPI"], S_donor["WE_ratio"])
results["pelka_tumor_donor_pseudobulk"] = {
    "n_donors": int(len(donor_pb)), "PPI_WE_rho": float(r_donor), "p": float(p_donor)}
print(f"[Pelka 肿瘤上皮 donor 伪 bulk] n={len(donor_pb)} PPI~WE rho={r_donor:+.3f} p={p_donor:.3e}")

# 正常上皮 donor 伪 bulk（对照）
norm = Xd[Xd["disease"] == "normal"]
norm_pb = norm.groupby("donor")[genes].mean()
if len(norm_pb) > 5:
    S_norm = score_from_logexpr(norm_pb)
    r_norm, p_norm = stats.spearmanr(S_norm["PPI"], S_norm["WE_ratio"])
    results["pelka_normal_donor_pseudobulk"] = {"n_donors": int(len(norm_pb)),
                                                 "PPI_WE_rho": float(r_norm), "p": float(p_norm)}
    print(f"[Pelka 正常上皮 donor 伪 bulk] n={len(norm_pb)} PPI~WE rho={r_norm:+.3f} p={p_norm:.3e}")

# MMR 分组伪 bulk（MMRd vs MMRp vs Normal 各一个 bulk 点，跨基因 z 对比）
if mmr is not None:
    mmr_groups = {}
    for m in ["Normal", "MMRp", "MMRd"]:
        sub_m = Xd[Xd["MMR"] == m]
        if len(sub_m) > 100:
            mmr_groups[m] = sub_m[genes].mean()
    if len(mmr_groups) == 3:
        mmr_E = pd.DataFrame(mmr_groups).T
        S_mmr = score_from_logexpr(mmr_E)
        mmr_we = {m: float(S_mmr.loc[m, "WE_ratio"]) for m in mmr_E.index}
        mmr_ppi = {m: float(S_mmr.loc[m, "PPI"]) for m in mmr_E.index}
        results["pelka_MMR_pseudobulk"] = {"WE_ratio": mmr_we, "PPI": mmr_ppi}
        print(f"[Pelka MMR 伪 bulk] WE_ratio={ {k: round(v,3) for k,v in mmr_we.items()} } PPI={ {k: round(v,3) for k,v in mmr_ppi.items()} }")

# ===================== Pelka 细胞类型伪 bulk =====================
FILES = [("pelka_epithelial.h5ad", "Epithelial"), ("pelka_tnk.h5ad", "T/NK"),
         ("pelka_myeloid.h5ad", "Myeloid"), ("pelka_b.h5ad", "B"),
         ("pelka_plasma.h5ad", "Plasma"), ("pelka_stromal.h5ad", "Stromal"),
         ("pelka_mast.h5ad", "Mast")]
ct_pb = {}
for fn, ct in FILES:
    p = os.path.join(SCR, fn)
    if not os.path.exists(p) or os.path.getsize(p) < 5e6:
        continue
    try:
        aa = ad.read_h5ad(p, backed="r")
    except Exception:
        continue
    s = aa.var["feature_name"].astype(str)
    s2, seen2 = [], set()
    for i, gg in enumerate(s):
        if gg in ALLG and gg not in seen2:
            s2.append(i); seen2.add(gg)
    sb = aa[:, s2].to_memory()
    Xb = sb.X.toarray() if not isinstance(sb.X, np.ndarray) else np.asarray(sb.X)
    gnames = [s.iloc[i] for i in s2]
    ct_pb[ct] = pd.Series(Xb.mean(axis=0), index=gnames)
    aa.file.close()

ct_E = pd.DataFrame(ct_pb).T  # 7 类型 × 基因
S_ct = score_from_logexpr(ct_E)
ct_we = {ct: float(S_ct.loc[ct, "WE_ratio"]) for ct in ct_E.index}
ct_ppi = {ct: float(S_ct.loc[ct, "PPI"]) for ct in ct_E.index}
results["pelka_celltype_pseudobulk"] = {"WE_ratio": ct_we, "PPI": ct_ppi}
print(f"[Pelka 细胞类型伪 bulk]")
for ct in ct_E.index:
    print(f"   {ct:12s} PPI={ct_ppi[ct]:+.3f} WE_ratio={ct_we[ct]:+.3f}")

# ===================== GSE132465 患者级伪 bulk =====================
ann = {}
with gzip.open(os.path.join(SCR, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"), "rt") as f:
    next(f)
    for line in f:
        parts = line.rstrip().split("\t")
        ann[parts[0]] = {"patient": parts[1], "class": parts[2], "cell_type": parts[4]}

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

valid_cells = sorted(valid.keys())
idx = [valid[c] + 1 for c in valid_cells]
expr = pd.DataFrame({g: [float(parts[i]) for i in idx] for g, parts in gene_data.items()},
                    index=valid_cells)
meta = pd.DataFrame([ann[c] for c in valid_cells], index=valid_cells)
E = np.log1p(expr)
E = E.join(meta)
tum_epi = E[(E["class"] == "Tumor") & (E["cell_type"] == "Epithelial cells")]
patient_pb = tum_epi.groupby("patient")[ALLG].mean()
S_pat = score_from_logexpr(patient_pb)
r_pat, p_pat = stats.spearmanr(S_pat["PPI"], S_pat["WE_ratio"])
results["gse_tumor_patient_pseudobulk"] = {"n_patients": int(len(patient_pb)),
                                            "PPI_WE_rho": float(r_pat), "p": float(p_pat)}
print(f"[GSE132465 肿瘤上皮 patient 伪 bulk] n={len(patient_pb)} PPI~WE rho={r_pat:+.3f} p={p_pat:.3e}")

json.dump(results, open(os.path.join(OUT, "sc_lactylation_pseudobulk_results.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False, default=str)
print("\nPSEUDOBULK DONE -> sc_lactylation_pseudobulk_results.json")
