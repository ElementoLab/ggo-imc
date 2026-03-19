# Technology Stack

**Analysis Date:** 2026-03-19

## Languages

**Primary:**
- Python 3.10 - Main pipeline language for analysis scripts
- R - For statistical analysis (asd.R, not currently in repository)

**Secondary:**
- YAML - Configuration format for project metadata

## Runtime

**Environment:**
- Python 3.10.15 (or higher, per pyproject.toml requires-python = ">=3.10")

**Package Manager:**
- pip (via setuptools, uv for dependency management)
- Lockfile: requirements.txt present

**Execution:**
- Snakemake workflow orchestration
- Docker/Apptainer container runtime for reproducibility
- Native conda environment (ggo_imc) for R dependencies

## Frameworks

**Core Bioinformatics:**
- `imc-analysis` (custom IMC analysis package) - Core analysis utilities and plotting
- `scanpy` (1.11.5) - Single-cell analysis and AnnData manipulation
- `anndata` (0.11.4) - Annotated data matrix format for storing analysis results

**Batch Correction & Integration:**
- `harmonypy` - Batch effect correction (Harmony integration)
- `leidenalg` (0.10.2) - Community detection clustering algorithm
- `squidpy` (1.6.1) - Spatial transcriptomics analysis

**Spatial & Domain Analysis:**
- `scimap` - Spatial cell interaction mapping

**Build/Workflow:**
- Snakemake - Workflow orchestration and reproducibility
- `setuptools` (>=61) - Package building

**Testing:**
- pytest (via tool.pytest.ini_options in pyproject.toml)
- Test directory at `tests/`

## Key Dependencies

**Critical:**
- `scanpy` (1.11.5) - Core single-cell data structure and preprocessing
- `anndata` (0.11.4) - Data format used throughout pipeline for HDF5-based analysis
- `imc-analysis` - Custom wrapper providing IMC-specific utilities and analysis methods
- `pandas` - Data manipulation and clinical metadata handling
- `numpy` (2.0.2) - Numerical computation foundation

**Scientific Computing:**
- `scipy` (1.15.3) - Statistical tests and scientific functions
- `matplotlib` (3.10.7) - Figure generation (PDF/raster output)
- `seaborn` (0.13.2) - Statistical visualization and heatmaps
- `tqdm` - Progress bars for long-running operations
- `requests` - HTTP downloads for configuration files

**Data Processing:**
- `harmonypy` - Harmony batch correction
- `leidenalg` (0.10.2) - Graph-based clustering
- `squidpy` (1.6.1) - Spatial analysis operations
- `scimap` - Spatial cell interaction scoring

**Image Processing:**
- `tifffile` - TIFF image reading/writing (for benchmark scripts)
- `scikit-image` - Image processing utilities

**Specialized:**
- `IPython` (8.18.1) - Interactive environment (pinned version)
- `geopandas` (1.1.2) - Geospatial operations (spatial cell coordinates)

## Configuration

**Environment:**
- `config.yaml` - Project-level Snakemake config specifying paths and container image
- `metadata/ggo_config.yml` - Downloadable configuration with panel-specific file paths, color palettes, and clinical metadata mappings
- Environment variables: None detected; configuration is YAML-based

**Build:**
- `pyproject.toml` - Modern Python project configuration
  - Build backend: setuptools
  - Python version requirement: >=3.10
  - pytest configuration with pythonpath
  - No additional project dependencies (inherits from sc-tools base + deconvolution extra)

**Container:**
- `containers/sc_tools.sif` - Apptainer/Singularity image containing sc-tools base + dependencies
- Used via `scripts/run_container.sh` wrapper for reproducible execution

## Platform Requirements

**Development:**
- Python 3.10+
- Snakemake installed
- Docker (macOS) or Apptainer/Singularity (Linux) for containerized execution
- Conda with ggo_imc environment for R support
- R and renv/Bioconductor packages (for asd.R statistical analysis)

**Production/Execution:**
- Container runtime (Docker or Apptainer)
- Internet access for downloading:
  - Config file: `https://wcm.box.com/shared/static/mdntp2xf9tjobxeidkw93mg8jysb7nh9.yml`
  - Backup h5ad files from Box (various URLs for Panel G/H data)
- ~200GB+ disk space for processed image data and analysis outputs

## Containerization

**Image:** `containers/sc_tools.sif` (Apptainer/Singularity format)

**Execution:**
- Command wrapper: `scripts/run_container.sh` (auto-detects Docker vs. Apptainer)
- Usage: `cd {repo_root} && ./scripts/run_container.sh {project_path} python {script}`
- Note: R scripts run natively in conda env (not containerized)

## Data Format

**Primary:**
- HDF5 via AnnData (.h5ad) - All single-cell data, expression matrices, metadata, UMAP coordinates
- CSV - ROI metadata (areas, mean intensity, variance, spillover, masks)
- TIFF - Source microscopy images and segmentation masks

**Secondary:**
- XLSX - Clinical annotations (De-identified clinical data)

---

*Stack analysis: 2026-03-19*
