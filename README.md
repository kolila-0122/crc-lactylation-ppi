# Pyruvate partitioning and lactylation writer/eraser machinery in colorectal cancer

Analysis code and intermediate results for the manuscript:

> **Pyruvate partitioning correlates with lactylation erasers, with sign-unstable writer correlations, in colorectal cancer cohorts**

This repository contains every script needed to reproduce all reported correlations, together with the intermediate result files (`results/`) that are the authoritative source for the numbers in the manuscript.

---

## Data sources

All datasets are public. No primary data are hosted here (they are large and should be obtained from the original repositories).

| Dataset | Accession / source | Used for |
|---|---|---|
| TCGA COADREAD (n = 592) | cBioPortal (PanCancer Atlas 2018), study `coadread_tcga_pan_can_atlas_2018` | bulk discovery cohort |
| GSE39582 (n = 585) | GEO `GSE39582` (Affymetrix GPL570) | independent bulk replication |
| Pelka CRC atlas | CELLxGENE, Pelka et al. *Cell* 2021 (DOI 10.1016/j.cell.2021.08.003) | single-cell pseudobulk |
| GSE132465 (SMC) | GEO `GSE132465` (Lee et al. *Nat Genet* 2020) | single-cell pseudobulk |
| CMS spatial atlas | CELLxGENE, *npj Precis Oncol* 2024 (DOI 10.1038/s41698-023-00488-4) | 13 Visium slices |

---

## Requirements

- Python 3.13
- Packages: `numpy`, `pandas`, `scipy`, `anndata`, `matplotlib`, `python-docx`

---

## Configuration

Scripts resolve paths in this order (no absolute paths are hard-coded):

- `results/` — first the `results/` folder next to the script, then the parent folder (so the code runs both in this repository layout and inside the manuscript submission package).
- Original data are located through environment variables (set them to wherever you downloaded the data):
  - `SCRNA_DIR` — Pelka `.h5ad` and GSE132465 matrix files
  - `SPATIAL_DIR` — Visium slice directory
  - `DATA_DIR` — GSE39582 annotation (`GPL570.annot.gz`) and series matrix
  - `RESULTS_DIR` — optional override for the output directory

---

## Scripts

| Script | Purpose | Output |
|---|---|---|
| `tcga_bulk_scores.py` | Fetches TCGA COADREAD expression from the cBioPortal REST API and computes PPI, writer, eraser and WE-difference scores (plus per-gene correlations) | `results/tcga_bulk_scores.json` |
| `sc_single_arm_correlations.py` | Single-cell pseudobulk PPI–writer / PPI–eraser / PPI–WE correlations (Pelka and GSE132465) | `results/sc_single_arm_correlations.json` |
| `gse39582_we.py` | GSE39582 bulk replication (probe mapping via GPL570, arm scores, per-gene correlations) | `results/gse39582_WE_validation.json` |
| `sc_lactylation_pseudobulk.py` | Cell-type and MMR-stratified pseudobulk scores | `results/sc_lactylation_pseudobulk_results.json` |
| `spatial_lactylation.py` | Spot-level PPI ~ WE-difference correlation across 13 Visium slices | `results/spatial_lactylation_results.csv` |
| `cohort_heterogeneity.py` | PPI distribution, power analysis, and Pelka MMR composition | `results/pelka_mmr_counts.json` |
| `make_all_figures.py` | Regenerates Fig1–Fig5 as 300 dpi PNG **and** LZW-TIFF | `figures/Fig1–5.{png,tif}` |
| `fetch_refs.py` | PubMed E-utilities verification of all 47 references | `refs_verified.json` |

**Run order for a full reproduction:** `tcga_bulk_scores.py` → `sc_single_arm_correlations.py` → `gse39582_we.py` → `sc_lactylation_pseudobulk.py` → `spatial_lactylation.py` → `cohort_heterogeneity.py` → `make_all_figures.py`.

> Note: `tcga_bulk_scores.py` and `sc_single_arm_correlations.py` need internet access (cBioPortal API) and the local data files respectively. `make_all_figures.py` only needs the `results/` files, and regenerates all five figures byte-for-byte.

---

## Notes on interpretation

- The writer and eraser scores are reported **separately**: the eraser correlation replicates (carried by HDAC1/HDAC2), whereas the hypothesised writer correlation is **not** reproduced across cohorts (it reverses sign in single-cell tumour epithelium). Both outcomes are reported as equivalent results in the manuscript.
- Scores are transcriptional proxies; they do not measure pyruvate flux or histone-lactylation abundance. No wet-lab validation is included.

## License

Code released for academic reuse. See the manuscript for citation details.
