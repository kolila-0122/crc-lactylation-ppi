# -*- coding: utf-8 -*-
"""GSE39582 独立队列验证：PPI ~ 乳酸化 writer/eraser 相关（补做，2026-10-13）

复用 ext_survival.py 的 series matrix 解析逻辑，补充 writer/eraser 基因的探针映射与相关分析。
"""
import os, gzip, re, json, time, urllib.request
import numpy as np
from scipy import stats

def _results_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.environ.get("RESULTS_DIR") or (
        os.path.join(here, "results") if os.path.isdir(os.path.join(here, "results"))
        else os.path.join(os.path.dirname(here), "results"))

# GPL570 注释与 series matrix 目录：用 DATA_DIR 环境变量指定（默认脚本上级 data/）
GEO_DIR = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
OUT = _results_dir()

WRITER = ["EP300", "CREBBP", "KAT8", "KAT2A"]
ERASER = ["SIRT1", "SIRT3", "HDAC1", "HDAC2", "HDAC3"]
OXG = ["MPC1", "MPC2", "PDHA1", "PDHB", "DLAT", "DLD", "PDHX"]
EFFG = ["LDHA", "LDHB", "SLC16A1", "SLC16A3", "SLC16A7"]
ALLG = WRITER + ERASER + OXG + EFFG

def download(url, outpath, retries=4):
    for attempt in range(retries):
        try:
            # 完整覆盖下载，不依赖断点续传（避免截断文件续传错位）
            if os.path.exists(outpath):
                os.remove(outpath)
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=180) as resp:
                with open(outpath, "wb") as f:
                    while True:
                        chunk = resp.read(65536)
                        if not chunk:
                            break
                        f.write(chunk)
            # 校验 gzip 完整性
            with gzip.open(outpath, "rb") as g:
                while g.read(1 << 20):
                    pass
            print("  downloaded + verified:", os.path.getsize(outpath), "bytes")
            return True
        except Exception as e:
            print("  attempt %d failed: %s" % (attempt + 1, e))
            if os.path.exists(outpath):
                os.remove(outpath)
            time.sleep(4)
    return False

# ---- GPL570 注释 ----
gpl_file = os.path.join(GEO_DIR, "GPL570.annot.gz")
def gzip_ok(p):
    try:
        with gzip.open(p, "rb") as g:
            while g.read(1 << 20):
                pass
        return True
    except Exception:
        return False

if not os.path.exists(gpl_file) or not gzip_ok(gpl_file):
    print("downloading GPL570 annotation...")
    ok = download("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz", gpl_file)
    print("download ok:", ok)
else:
    print("GPL570 annotation exists:", os.path.getsize(gpl_file))

print("parsing GPL570 annotation...")
probe_to_gene = {}
with gzip.open(gpl_file, "rt", encoding="utf-8", errors="replace") as f:
    header_seen = False
    id_col = gene_col = -1
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("!") or line.startswith("#") or line.startswith("^"):
            continue
        parts = line.split("\t")
        if not header_seen:
            for j, h in enumerate(parts):
                hl = h.strip('"').strip().lower()
                if hl == "id":
                    id_col = j
                elif hl == "gene symbol":
                    gene_col = j
            header_seen = True
            continue
        if gene_col >= 0 and id_col >= 0 and len(parts) > max(id_col, gene_col):
            pid = parts[id_col].strip('"')
            gene = parts[gene_col].strip('"').split("///")[0].strip()
            if gene and gene not in ("---", ""):
                probe_to_gene[pid] = gene
print("mapped probes:", len(probe_to_gene))

# ---- series matrix ----
print("parsing series matrix...")
matrix_file = os.path.join(GEO_DIR, "GSE39582_series_matrix.txt.gz")
with gzip.open(matrix_file, "rt", encoding="utf-8", errors="replace") as f:
    content = f.read()
lines = content.split("\n")
header_line = None
data_lines = []
in_data = False
for line in lines:
    if line.startswith("!series_matrix_table_begin"):
        in_data = True
    elif in_data and line.startswith('"ID_REF"'):
        header_line = line.strip().split("\t")
    elif in_data and line.startswith('"') and header_line:
        data_lines.append(line.strip().split("\t"))
    elif in_data and line.startswith("!series_matrix_table_end"):
        break
col_names = [h.strip('"') for h in header_line[1:]] if header_line else []
print("samples:", len(col_names), "| data lines:", len(data_lines))

expr = {}
for parts in data_lines:
    pid = parts[0].strip('"')
    expr[pid] = [float(v) if v not in ("", "NA", "null") else np.nan for v in parts[1:]]

# 基因 -> 探针（选 IQR 最大）
gene_to_probes = {}
for p, g in probe_to_gene.items():
    gene_to_probes.setdefault(g, []).append(p)

gene_probe = {}
for g in ALLG:
    probes = [p for p in gene_to_probes.get(g, []) if p in expr]
    best, best_iqr = None, -1
    for p in probes:
        v = np.array(expr[p]); v = v[~np.isnan(v)]
        if len(v) > 10:
            iqr = np.percentile(v, 75) - np.percentile(v, 25)
            if iqr > best_iqr:
                best_iqr, best = iqr, p
    gene_probe[g] = best
    print("  %-8s -> %s" % (g, best or "MISSING"))

# z-score 每基因跨样本
gene_z = {}
for g, p in gene_probe.items():
    if p:
        v = np.array(expr[p])
        mu, sd = np.nanmean(v), np.nanstd(v)
        gene_z[g] = (v - mu) / sd if sd > 0 else np.zeros_like(v)

def score(genes):
    gs = [g for g in genes if g in gene_z]
    if not gs:
        return None
    return np.nanmean(np.array([gene_z[g] for g in gs]), axis=0)

ppi = score(OXG) - score(EFFG)
writer = score(WRITER)
eraser = score(ERASER)
we = writer - eraser

def corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return stats.spearmanr(a[ok], b[ok])

res = {}
for name, s in [("writer", writer), ("eraser", eraser), ("WE", we)]:
    r, p = corr(ppi, s)
    res[name] = {"rho": float(r), "p": float(p)}
    print("GSE39582 PPI ~ %s: rho=%+.3f p=%.2e" % (name, r, p))

# 单基因 rho（GSE39582）
single_gene = {}
for g in WRITER + ERASER:
    if g in gene_z:
        r, p = corr(ppi, gene_z[g])
        single_gene[g] = {"rho": float(r), "p": float(p)}
res["single_gene"] = single_gene

# 也输出 writer/eraser 基因的覆盖
res["coverage"] = {g: (gene_probe[g] is not None) for g in ALLG}
res["n_samples"] = int((np.isfinite(ppi)).sum())
res["writer_genes_found"] = [g for g in WRITER if g in gene_z]
res["eraser_genes_found"] = [g for g in ERASER if g in gene_z]

json.dump(res, open(os.path.join(OUT, "gse39582_WE_validation.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
print("\nDONE -> gse39582_WE_validation.json")
print(json.dumps(res, indent=1))
