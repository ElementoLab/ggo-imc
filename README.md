# ggo-imc

**Spatial profiling of early-stage lung adenocarcinoma reveals patterns of immunomodulation and epithelial plasticity**

Kim, Ravichandran, Yoffe et al. — *Cancer Cell* (2025)

This repository contains code to reproduce all main figures from the manuscript. The study uses Imaging Mass Cytometry (IMC) to profile 2.24 million cells across 122 early-stage lung adenocarcinoma specimens, identifying two progression trajectories (immune-inflamed and fibrotic) and a diagnostic blind spot in which 20.5% of fibrotic tumors are radiologically misclassified.

> **Citation:**
> Kim J, Ravichandran S, Yoffe L, et al. Spatial profiling of early-stage lung adenocarcinoma reveals patterns of immunomodulation and epithelial plasticity. *Cancer Cell*. 2025. *(citation to be updated upon publication)*

---

## Data

Processed and raw data are available on Zenodo (access may be restricted prior to publication):

- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14822106.svg)](https://doi.org/10.5281/zenodo.14822106) — Processed Data
- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14822106.svg)](https://doi.org/10.5281/zenodo.14822106) — Raw Data

---

## Environment Setup

```bash
conda create -n ggo-imc python=3.9 -y
conda activate ggo-imc
pip install -r requirements.txt
```

---

## Reproducing Figures

Run from the repository root using Snakemake. The pipeline auto-detects the container runtime (Docker on macOS, Apptainer on Linux/HPC).

**All figures:**
```bash
snakemake -d . -s Snakefile all
```

**Individual figures:**
```bash
snakemake -d . -s Snakefile figure1
snakemake -d . -s Snakefile figure2
snakemake -d . -s Snakefile figure3
snakemake -d . -s Snakefile figure4
snakemake -d . -s Snakefile figure5
```

**Reset sentinel files:**
```bash
snakemake -d . -s Snakefile clean
```

**Local conda (no container):**
```bash
conda activate ggo_imc
SC_TOOLS_RUNTIME=none snakemake -d . -s Snakefile figure1
```

> **Note:** `snakemake figure5` runs `patient_group.py` and `roi_pca_plot_group.py` (Python, no R required).
> The `patient` rule (`scripts/asd.R`) is a separate target not included in `figure5` or `all` — it requires R and is not yet committed to this repo. The R-dependent panels of Figure 5 cannot be reproduced from this repo until `asd.R` is added.

---

## Figure-to-Script Mapping

| Figure | Content | Script(s) |
|--------|---------|-----------|
| **Figure 1** | Cohort overview, cell type heatmap (1a-d), ROI PCA archetypes (1e) | `celltype_heatmap_info.py`, `roi_pca_plot.py` |
| **Figure 2** | Immune dynamics: lymphocyte abundance across stages (2a-h), myeloid/macrophage polarization (2i-m) | `celltype_differential_abundance.py`, `t_cell_analysis.py`, `myeloid_analysis.py` |
| **Figure 3** | Stromal expansion, epithelial remodeling, EMP (3a-g) | `celltype_differential_abundance.py`, `epithelial_characterization.py` |
| **Figure 4** | UTAG microenvironments, TLS, tumor-stroma interface, cell-cell interactions (4a-f) | `ue_analysis.py` |
| **Figure 5** | Patient risk groups via hierarchical clustering, fibrotic trajectory, diagnostic gap (5a-e) | `roi_pca_plot_group.py`, `patient_group.py`, `asd.R` |

---

## Verifying Reproducibility

A checksum manifest of all pipeline outputs is provided in `figures_checksums.md5`. After running the pipeline:

```bash
python3 -c "
import hashlib
from pathlib import Path
exts = {'.pdf', '.png', '.svg'}
for line in open('figures_checksums.md5'):
    md5, path = line.strip().split(' ', 1)
    p = Path(path)
    if p.exists():
        actual = hashlib.md5(p.read_bytes()).hexdigest()
        status = 'OK' if actual == md5 else 'CHANGED'
        if status != 'OK':
            print(f'{status}: {path}')
"
```

No output means all figures are identical to the verified run.

---

## Repository Structure

```
scripts/
  # Pipeline — Figures 1-5
  celltype_heatmap_info.py           # Figure 1: cell type heatmap
  roi_pca_plot.py                    # Figure 1: ROI PCA archetypes
  celltype_differential_abundance.py # Figures 2, 3: immune/stromal density
  t_cell_analysis.py                 # Figure 2: lymphocyte functional states
  myeloid_analysis.py                # Figure 2: myeloid/macrophage polarization
  epithelial_characterization.py     # Figure 3: epithelial phenotypes, EMP
  ue_analysis.py                     # Figure 4: UTAG microenvironments, TLS
  roi_pca_plot_group.py              # Figure 5: patient group PCA overlay
  patient_group.py                   # Figure 5: patient risk stratification
  # Utilities
  download_yaml.py                   # Download data from Zenodo/Box
  concat_anndata.py                  # Concatenate per-sample AnnData
  label_metadata_anndata.py          # Attach clinical metadata
  generate_manifests.py              # Generate batch manifests
  # Benchmarks
  benchmark_integration.py           # Integration method comparison
  benchmark_segmentation.py          # Segmentation benchmark
  # Archived
  exploratory/                       # Exploratory scripts (not in pipeline)
```
