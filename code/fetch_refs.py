# -*- coding: utf-8 -*-
"""用 NCBI E-utilities 检索核实子课题②关键文献的真实 PMID/标题/期刊/年份"""
import json, time, urllib.parse, urllib.request, os

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = os.path.dirname(os.path.abspath(__file__))
UA = {"User-Agent": "Mozilla/5.0 (subproject2-ref-check)"}

def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")

def esearch(term, retmax=4):
    p = {"db": "pubmed", "term": term, "retmax": str(retmax), "retmode": "json", "sort": "relevance"}
    return json.loads(get(BASE + "esearch.fcgi?" + urllib.parse.urlencode(p)))["esearchresult"]

def esummary(ids):
    if not ids:
        return {}
    p = {"db": "pubmed", "id": ",".join(ids), "retmode": "json"}
    return json.loads(get(BASE + "esummary.fcgi?" + urllib.parse.urlencode(p)))["result"]

QUERIES = [
    ("乳酸化发现(2019 Nature)", "histone lactylation[Title] AND metabolic regulation[Title/Abstract] AND 2019[dp]"),
    ("铜死亡发现(2022 Science)", "copper induces cell death[Title] AND lipoylated[Title/Abstract]"),
    ("p300/CBP 乳酸化 writer", "p300[Title/Abstract] AND lactylation[Title] AND (writer[Title/Abstract] OR acetyltransferase[Title/Abstract])"),
    ("KAT8 乳酸化(PNAS)", "KAT8[Title/Abstract] AND lactylation[Title/Abstract]"),
    ("METTL3 乳酸化(Mol Cell)", "METTL3[Title/Abstract] AND lactylation[Title/Abstract]"),
    ("PKM2 乳酸化", "PKM2[Title/Abstract] AND lactylation[Title/Abstract]"),
    ("HDAC/SIRT 乳酸化 eraser", "lactylation[Title] AND (deacetylation[Title/Abstract] OR eraser[Title/Abstract] OR SIRT[Title/Abstract] OR HDAC[Title/Abstract])"),
    ("乳酸化 免疫 肿瘤 综述", "lactylation[Title] AND (tumor[Title/Abstract] OR cancer[Title/Abstract] OR immunity[Title/Abstract]) AND review[pt]"),
    ("HIF-1α 铜死亡抵抗(Cancer Cell)", "HIF-1alpha[Title/Abstract] AND cuproptosis[Title/Abstract]"),
    ("乳酸化 结直肠癌", "lactylation[Title/Abstract] AND colorectal cancer[Title/Abstract]"),
    ("Warburg 效应经典", "Warburg[Title] AND metabolism of tumors[Title/Abstract]"),
    ("线粒体丙酮酸载体 肿瘤", "mitochondrial pyruvate carrier[Title/Abstract] AND cancer[Title/Abstract] AND review[pt]"),
]

results = []
for label, q in QUERIES:
    try:
        r = esearch(q)
        ids = r.get("idlist", [])
        summ = esummary(ids)
        rows = []
        for pid in ids:
            s = summ.get(pid, {})
            rows.append({
                "pmid": pid,
                "title": s.get("title", ""),
                "journal": s.get("fulljournalname", s.get("source", "")),
                "year": s.get("pubdate", "")[:4],
            })
        results.append({"label": label, "query": q, "count": r.get("count", 0), "hits": rows})
        print(f"[{label}] count={r.get('count')}")
        for x in rows:
            print(f"   PMID {x['pmid']} | {x['year']} | {x['journal']} | {x['title'][:100]}")
    except Exception as e:
        print(f"[{label}] ERROR: {e}")
    time.sleep(0.4)

with open(os.path.join(OUT, "refs_verified.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=1)
print("\nDONE -> refs_verified.json")
