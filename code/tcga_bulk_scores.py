# -*- coding: utf-8 -*-
"""TCGA COADREAD bulk 复合评分（产出 PPI、writer、eraser、WE 及相关系数）

修复 B1：替代对旧项目目录 novel_angle_validation.json 的外部依赖，本脚本为权威计算来源。
从 cBioPortal REST API 抓取 PPI 双臂基因 + 乳酸化 writer/eraser 基因的表达，计算复合评分。
Entrez Gene ID 为稳定公开标识，硬编码于下方 GENE_ENTREZ。

输出：results/tcga_bulk_scores.json（每样本 PPI/writer/eraser/WE + 各轴相关系数）。
"""
import os, json, time, urllib.parse, urllib.request
import numpy as np
from scipy import stats

CB = "https://www.cbioportal.org/api/"
STUDY = "coadread_tcga_pan_can_atlas_2018"

def _results_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.environ.get("RESULTS_DIR") or (
        os.path.join(here, "results") if os.path.isdir(os.path.join(here, "results"))
        else os.path.join(os.path.dirname(here), "results"))

OUT = _results_dir()
os.makedirs(OUT, exist_ok=True)

OXG = ["MPC1", "MPC2", "PDHA1", "PDHB", "DLAT", "DLD", "PDHX"]
EFFG = ["LDHA", "LDHB", "SLC16A1", "SLC16A3", "SLC16A7"]
WRITER = ["EP300", "CREBBP", "KAT8", "KAT2A"]
ERASER = ["SIRT1", "SIRT3", "HDAC1", "HDAC2", "HDAC3"]
TARGET = ["LDHA", "PKM", "HIF1A", "MDH2", "ALDH2", "FASN", "GCLC", "METTL3", "EEF1A2"]

GENE_ENTREZ = {
    "MPC1": 51660, "MPC2": 25874, "PDHA1": 5160, "PDHB": 5162, "DLAT": 1737, "DLD": 1738, "PDHX": 8050,
    "LDHA": 3939, "LDHB": 3945, "SLC16A1": 6566, "SLC16A3": 9123, "SLC16A7": 9194,
    "EP300": 2033, "CREBBP": 1387, "KAT8": 84148, "KAT2A": 2648,
    "SIRT1": 23411, "SIRT3": 23410, "HDAC1": 3065, "HDAC2": 3066, "HDAC3": 8841,
    "PKM": 5315, "HIF1A": 3091, "MDH2": 4191, "ALDH2": 217,
    "FASN": 2194, "GCLC": 2729, "METTL3": 56339, "EEF1A2": 1917,
}


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (tcga-bulk-scores)"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read().decode("utf-8", "replace")


def fetch_gene(mrna, entrez_id):
    url = CB + "molecular-profiles/%s/molecular-data?" % mrna + urllib.parse.urlencode(
        {"sampleListId": mrna, "entrezGeneId": entrez_id, "pageSize": 100000, "projection": "SUMMARY"})
    rows = json.loads(get(url))
    return {r["sampleId"]: r["value"] for r in rows if "sampleId" in r and "value" in r}


def zscore(x):
    x = np.asarray(x, float)
    m, sd = np.nanmean(x), np.nanstd(x)
    return (x - m) / sd if sd > 0 else np.zeros_like(x)


def main():
    profs = json.loads(get(CB + "studies/%s/molecular-profiles" % STUDY))
    mrna = next(p["molecularProfileId"] for p in profs if p["molecularProfileId"].endswith("_rna_seq_v2_mrna"))
    print("mrna profile:", mrna)

    gene_expr = {}
    for g, eid in GENE_ENTREZ.items():
        gene_expr[g] = fetch_gene(mrna, eid)
        time.sleep(0.15)
    samples = sorted(gene_expr["MPC1"].keys())
    print("samples:", len(samples))

    def arm_score(genes):
        mat = np.array([[np.log1p(float(gene_expr[g].get(s, np.nan))) for s in samples] for g in genes])
        z = np.apply_along_axis(zscore, axis=1, arr=mat)
        return np.nanmean(z, axis=0)

    ppi = arm_score(OXG) - arm_score(EFFG)
    writer = arm_score(WRITER)
    eraser = arm_score(ERASER)
    we = writer - eraser

    def corr(a, b):
        return stats.spearmanr(a, b)

    # 单基因 rho（PPI vs 每个 writer/eraser 基因）
    single_gene = {}
    for g in WRITER + ERASER:
        zg = zscore(np.array([np.log1p(float(gene_expr[g].get(s, np.nan))) for s in samples]))
        r, p = corr(ppi, zg)
        single_gene[g] = {"rho": float(r), "p": float(p)}

    # lactylation target-gene score
    target_score = arm_score(TARGET)
    tr, tp = corr(ppi, target_score)

    res = {
        "n_samples": int(len(samples)),
        "writer_rho": float(corr(ppi, writer)[0]), "writer_p": float(corr(ppi, writer)[1]),
        "eraser_rho": float(corr(ppi, eraser)[0]), "eraser_p": float(corr(ppi, eraser)[1]),
        "WE_rho": float(corr(ppi, we)[0]), "WE_p": float(corr(ppi, we)[1]),
        "target_rho": float(tr), "target_p": float(tp),
        "single_gene": single_gene,
        "ppi": {s: float(ppi[i]) for i, s in enumerate(samples)},
    }
    with open(os.path.join(OUT, "tcga_bulk_scores.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, ensure_ascii=False)
    print("TCGA bulk:", {k: round(res[k], 3) for k in ["writer_rho", "eraser_rho", "WE_rho", "target_rho"]})
    for g in WRITER + ERASER:
        print("  %-8s rho=%+.3f p=%.2e" % (g, single_gene[g]["rho"], single_gene[g]["p"]))
    print("DONE -> results/tcga_bulk_scores.json")


if __name__ == "__main__":
    main()
