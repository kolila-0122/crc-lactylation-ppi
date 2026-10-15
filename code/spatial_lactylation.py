# -*- coding: utf-8 -*-
"""子课题② 空间层乳酸化 writer/eraser 评分补跑（Visium 13 张有效切片）

核心假设：低 PPI 生态位（外排型）= 乳酸化 writer 高负荷（WE_ratio 高）
验证：(1) spot 级 PPI vs WE_ratio 负相关；(2) WE_ratio 空间自相关；(3) 低 PPI spot 的 WE_ratio 更高。
"""
import os, json, glob, warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
SPD = os.environ.get("SPATIAL_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scRNA", "spatial"))

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
ALLG = sorted(set(WRITER + ERASER + OXG + EFFG))

def morans_i(x, coords, n_perm=99, seed=42):
    tree = cKDTree(coords)
    d, _ = tree.query(coords, k=2)
    eps = np.median(d[:, 1]) * 1.5
    W = tree.query_pairs(eps, output_type="ndarray")
    i_, j_ = W[:, 0], W[:, 1]
    xc = x - x.mean(); n = len(x)
    S0 = len(i_) * 2.0
    I = (n / S0) * (2 * np.sum(xc[i_] * xc[j_])) / np.sum(xc ** 2)
    rng = np.random.RandomState(seed)
    perm = np.array([(n / S0) * (2 * np.sum((x[rng.permutation(n)] - x.mean())[i_] *
                                            (x[rng.permutation(n)] - x.mean())[j_])) / np.sum((x - x.mean()) ** 2)
                     for _ in range(n_perm)])
    p = (np.sum(np.abs(perm) >= abs(I)) + 1) / (n_perm + 1)
    return I, p

per_slice = {}
for f in sorted(glob.glob(os.path.join(SPD, "*.h5ad"))):
    sid = os.path.basename(f).split(".")[0][:8]
    if os.path.getsize(f) < 1e6:
        continue
    try:
        a = ad.read_h5ad(f)
    except Exception:
        print("skip broken", sid); continue
    a = a[a.obs["in_tissue"].values == 1]
    sym = a.var["feature_name"].astype(str)
    sel, seen = [], set()
    for i, s in enumerate(sym):
        if s in ALLG and s not in seen:
            sel.append(i); seen.add(s)
    genes = [sym.iloc[i] for i in sel]
    Xs = a.X[:, sel]
    Xs = Xs.toarray() if not isinstance(Xs, np.ndarray) else Xs
    lib = np.asarray(a.X.sum(axis=1)).ravel()
    Xn = np.log1p(Xs / np.maximum(lib, 1)[:, None] * 1e4)
    Xd = pd.DataFrame(Xn, columns=genes)
    Z = ((Xd - Xd.mean()) / Xd.std().replace(0, np.nan)).fillna(0)

    w = [g for g in WRITER if g in Z.columns]
    e = [g for g in ERASER if g in Z.columns]
    ox = [g for g in OXG if g in Z.columns]
    ef = [g for g in EFFG if g in Z.columns]
    writer = Z[w].mean(axis=1).values
    eraser = Z[e].mean(axis=1).values
    WE = writer - eraser
    PPI = Z[ox].mean(axis=1).values - Z[ef].mean(axis=1).values

    coords = a.obsm["spatial"]
    r, pv = stats.spearmanr(PPI, WE)
    I, pI = morans_i(WE, coords)
    med = np.median(PPI)
    low = PPI < med
    u, up = stats.mannwhitneyu(WE[low], WE[~low])
    per_slice[sid] = {
        "n_spots": int(a.n_obs), "PPI_WE_rho": float(r), "PPI_WE_p": float(pv),
        "WE_morans_I": float(I), "WE_morans_p": float(pI),
        "WE_lowPPI_median": float(np.median(WE[low])), "WE_highPPI_median": float(np.median(WE[~low])),
        "WE_low_vs_high_mw_p": float(up),
    }
    print(f"{sid}: n={a.n_obs} PPI~WE rho={r:+.3f} p={pv:.2e} | WE MoranI={I:+.3f} p={pI:.2e} | lowPPI_WE={np.median(WE[low]):+.3f} highPPI_WE={np.median(WE[~low]):+.3f} mw_p={up:.2e}", flush=True)
    a.file.close() if hasattr(a, "file") else None

res = pd.DataFrame(per_slice).T
res.to_csv(os.path.join(OUT, "spatial_lactylation_results.csv"), encoding="utf-8-sig")

summary = {
    "n_slices": int(len(res)),
    "PPI_WE_rho_mean": float(res["PPI_WE_rho"].mean()),
    "PPI_WE_rho_neg_frac": float((res["PPI_WE_rho"] < 0).mean()),
    "PPI_WE_rho_significant_frac": float((res["PPI_WE_p"] < 0.05).mean()),
    "WE_morans_sig_frac": float((res["WE_morans_p"] < 0.05).mean()),
    "WE_lowPPI_higher_frac": float((res["WE_lowPPI_median"] > res["WE_highPPI_median"]).mean()),
    "WE_lowPPI_median_mean": float(res["WE_lowPPI_median"].mean()),
    "WE_highPPI_median_mean": float(res["WE_highPPI_median"].mean()),
    "per_slice": {sid: {k: float(v) if isinstance(v, (int, float)) else v for k, v in row.items()}
                  for sid, row in per_slice.items()},
}
json.dump(summary, open(os.path.join(OUT, "spatial_lactylation_summary.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False, default=str)

print("\n=== META SUMMARY ===")
print("PPI~WE rho mean:", round(summary["PPI_WE_rho_mean"], 3),
      "| 负相关切片占比:", summary["PPI_WE_rho_neg_frac"],
      "| 显著切片占比:", summary["PPI_WE_rho_significant_frac"])
print("WE 空间自相关显著占比:", summary["WE_morans_sig_frac"])
print("低 PPI spot WE 更高占比:", summary["WE_lowPPI_higher_frac"])
print("低 PPI WE 均值:", round(summary["WE_lowPPI_median_mean"], 3),
      "vs 高 PPI WE 均值:", round(summary["WE_highPPI_median_mean"], 3))
print("SPATIAL LACTYLATION DONE")
