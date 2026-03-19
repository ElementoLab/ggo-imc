# Technology Stack

**Analysis Date:** 2026-03-19

## Languages

**Primary:**
- Python 3.9+ - All pipeline scripts, data processing, visualization
- Snakemake (Python-based workflow) - Orchestration of figure generation pipeline

**Secondary:**
- R - Required for `scripts/asd.R` (hierarchical clustering, patient stratification; not yet committed to repo)
- Shell/Bash - Container orchestration via `scripts/run_container.sh`

## Runtime

**Environment:**
- Python 3.9 (specified in README; pyproject.toml requires >=3.10)
- Conda for environment management

**Package Manager:**
- pip (via conda environment)
- Lockfile: requirements.txt (present but minimal; most dependencies come from parent sc-tools project)

## Frameworks

**Core Data Science:**
- scanpy 1.x - Single-cell/spatial data analysis
- AnnData - In-memory representation of annotated data matrices
- pandas - Tabular data manipulation
- numpy - Numerical computations
- scipy - Scientific computing (ndi, spatial operations)

**Domain-Specific:**
- imc-analysis - Custom IMC (Imaging Mass Cytometry) analysis utilities (external package)
- squidpy 1.6.1 - Spatial data analysis for cells, neighborhoods, and graphs
- scimap - Single-cell image analysis platform (cell phenotyping, segmentation)
- harmonypy - Batch integration via Harmony algorithm
- leidenalg 0.10.2 - Community detection (Leiden clustering algorithm)

**Visualization:**
- matplotlib - Figure generation and saving
- seaborn - Statistical graphics and plotting helpers
- adjustText - Label placement optimization on scatter plots

**Testing:**
- pytest - Test discovery and execution (config: `pyproject.toml` [tool.pytest.ini_options])

**Utilities:**
- tqdm - Progress bars for long-running operations
- requests - HTTP requests (data download from Box)
- tifffile - TIFF image I/O for microscopy data
- ipython 8.18.1 - Interactive shell support

## Key Dependencies

**Critical:**
- imc-analysis - Provides IMC-specific tools (`imc.utils.parse_yaml`, `imc.pl.celltype_heatmap`, `imc.tl.celltype_density`)
- scanpy - Core single-cell analysis; used in every figure script
- AnnData - Data structure for storing preprocessed panels (PANEL_G, PANEL_H)
- squidpy 1.6.1 - UTAG microenvironment analysis (Figure 4: `ue_analysis.py`)
- harmonypy - Batch integration for combining panels
- leidenalg 0.10.2 - Graph-based clustering for cell types and communities

**Infrastructure:**
- pandas - Metadata merging, obs/var manipulation
- numpy - Array operations, statistical functions
- scipy - Neighbor calculations, morphology operations
- matplotlib/seaborn - Multi-panel figure generation with publication-quality output

## Configuration

**Environment:**
- `.env` files: Not detected; secrets/credentials managed at container level
- Snakemake config: `config.yaml`
  - `repo_root`: Path to sc-tools root (../../..)
  - `project_rel`: Relative path from repo root to project (projects/imc/ggo-imc)
  - `container_sif`: Path to sc_tools Apptainer/Docker image (containers/sc_tools.sif)

**Data Configuration:**
- Metadata YAML: `metadata/ggo_config.yml` (downloaded via `scripts/download_yaml.py`)
  - Downloaded from Box: https://wcm.box.com/shared/static/mdntp2xf9tjobxeidkw93mg8jysb7nh9.yml
  - Contains panel names (PANEL_G, PANEL_H), AnnData paths, backup URLs, cell type marker groups
  - Parsed by `utils.load_config()` at script initialization

**Build:**
- `pyproject.toml` - Project metadata, pytest config, no additional build steps
- `requirements.txt` - Minimal project-specific requirements (most deps from parent sc-tools)
- Snakefile - Workflow rules for figures 1-5, data download, cleanup

## Platform Requirements

**Development:**
- macOS or Linux with Python 3.9+
- Conda or similar environment manager
- Either Docker (macOS) or Apptainer (Linux/HPC)
- For Figure 5 R component: R with renv/Bioconductor packages in same conda env

**Production:**
- Apptainer container (sc_tools.sif) for consistent reproducibility
- Zenodo/Box access for downloading processed and raw data
- Disk space for 122 specimens × 2.24M cells (~multiple GB per panel)

## Container Strategy

**Runtime Auto-Detection:**
- `scripts/run_container.sh` selects container based on platform:
  - Docker on macOS
  - Apptainer on Linux/HPC
- Enables reproducible execution without local package conflicts
- sc_tools image includes all Python dependencies except R

**Local Alternative:**
- `SC_TOOLS_RUNTIME=none snakemake` bypasses containerization and uses local conda env
- Required for R-dependent rules (patient rule runs Rscript directly)

---

*Stack analysis: 2026-03-19*
