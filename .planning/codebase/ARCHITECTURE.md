# Architecture

**Analysis Date:** 2026-03-19

## Pattern Overview

**Overall:** Pipeline-driven analysis pattern with modular analysis scripts orchestrated by Snakemake.

**Key Characteristics:**
- Declarative workflow definition (Snakemake) coordinating independent analysis scripts
- Shared utilities layer (`utils.py`) providing common operations across all scripts
- Figure-per-rule organization: each Snakemake rule produces one or more complete figures
- Config-driven data loading and artifact paths (YAML-based metadata)
- AnnData-centric data structures (single-cell/cell-type abundance data)
- External dependency on `imc_analysis` package for domain-specific visualization and statistics

## Layers

**Orchestration Layer:**
- Purpose: Define and coordinate execution of analysis tasks
- Location: `Snakefile`
- Contains: 12 rules organized by figure (figure1–figure5 + utility rules)
- Depends on: Script files, Snakemake configuration
- Used by: Command-line invocation via `snakemake`

**Analysis Scripts Layer:**
- Purpose: Implement individual analyses that produce one or more figures
- Location: `scripts/*.py` (9 pipeline scripts for figures 1–5)
- Contains: Standalone analysis scripts that load data, compute statistics, and generate visualizations
- Depends on: Config, AnnData files, `utils.py`, `imc_analysis` and `scanpy` libraries
- Used by: Snakemake rules

**Shared Utilities Layer:**
- Purpose: Provide common operations to avoid duplication across analysis scripts
- Location: `scripts/utils.py`
- Contains: Config loading, AnnData loading, figure saving, directory creation
- Depends on: `scanpy`, `imc_analysis`, `matplotlib`
- Used by: All 9 pipeline scripts

**Configuration Layer:**
- Purpose: Centralize data paths, metadata, color schemes, and analysis parameters
- Location: `metadata/ggo_config.yml` (downloaded at runtime) + `config.yaml` (Snakemake config)
- Contains: Panel definitions, AnnData paths, condition labels, color palettes, marker lists
- Depends on: Box/Zenodo for backup downloads
- Used by: All analysis scripts via `load_config()`

**Test Layer:**
- Purpose: Validate utility functions and pure business logic
- Location: `tests/*.py`
- Contains: Unit tests for `utils.py` functions and `patient_group.py` functions
- Depends on: pytest, pandas
- Used by: pytest runner

## Data Flow

**Figure Generation Pipeline:**

1. **Download** (`download_yaml.py` rule)
   - Downloads `ggo_config.yml` from Box/Zenodo
   - Provides metadata for all downstream analysis

2. **Load Data** (all analysis scripts)
   - Each script calls `load_config()` → loads YAML metadata
   - Calls `load_panels()` or `load_single_panel()` → reads AnnData from disk or backup URL
   - Optionally reads pre-computed subsets (e.g., myeloid cells, lymphocytes)

3. **Compute** (analysis-specific)
   - Cell type density calculations (`imc.tl.celltype_density`)
   - PCA and dimensionality reduction (`scanpy.tl.pca`, `scanpy.tl.umap`)
   - Statistical tests (Mann-Whitney U via `imc.tl.grouped_mwu_test`)
   - Custom aggregations (cell proportions, marker co-expression)

4. **Visualize** (analysis-specific)
   - Generate matplotlib figures using `scanpy.pl` and `imc_analysis.pl` functions
   - Heatmaps, UMAPs, PCA plots, violin plots, density plots
   - Custom scatter plots and regression plots

5. **Save** (all analysis scripts)
   - Figures saved via `save_figure()` → creates directories and saves PDFs/PNGs
   - Output collected in `figures/figure{1-5}/` hierarchy
   - Intermediate AnnData objects cached in `processed/` or metadata/

**State Management:**
- No stateful services; all state is file-based (AnnData HDF5 objects, PNG/PDF files)
- Sentinel `.done` files used by Snakemake to track rule completion
- Reproducibility via checksum verification (`figures_checksums.md5`)

## Key Abstractions

**Config Metadata Dict:**
- Purpose: Centralized source of truth for paths, colors, markers
- Examples: `scripts/celltype_heatmap_info.py`, `scripts/roi_pca_plot.py`
- Pattern: `metadata = load_config()` → dictionary access like `metadata[panel]['AnnData']['phenotyped_umap_name']`
- Used for: AnnData paths, panel definitions, condition labels, color palettes

**AnnData (Anndata Objects):**
- Purpose: Represent cell-type abundance or single-cell data with metadata
- Examples: `PANEL_G` and `PANEL_H` phenotyped UMAPs, lymphocyte/myeloid subsets
- Pattern: Loaded via `load_panels()` or `load_single_panel()`, manipulated using `scanpy` API
- Contains: `.X` (expression/abundance matrix), `.obs` (cell/ROI metadata), `.var` (markers/cell types), `.uns` (color palettes, parameters)

**Condition Keys:**
- Purpose: Grouping/stratification variables for density calculations and comparisons
- Examples: 'pathology', 'radio' (radiology), 'Group' (patient risk group)
- Pattern: Stored in `metadata['CONDITIONS']` or hardcoded in scripts
- Used for: Subsetting, statistical tests, conditional visualization

**Density Matrices:**
- Purpose: Summarize cell-type abundance per ROI and condition
- Pattern: Output of `imc.tl.celltype_density()` → AnnData where rows=ROIs, cols=cell types, values=counts/density
- Used for: Plotting density distributions, statistical testing (Mann-Whitney U), correlation analysis

## Entry Points

**Snakemake Orchestration:**
- Location: `Snakefile` (root)
- Triggers: User runs `snakemake -d . -s Snakefile [target]`
- Responsibilities: Define rules, dependency management, rule execution via `run_container()`

**Analysis Scripts:**
- Location: `scripts/*.py` (9 pipeline scripts)
- Triggers: Invoked by Snakemake rules or directly for development/debugging
- Responsibilities: Load data, compute figures, save outputs

**Download Script:**
- Location: `scripts/download_yaml.py`
- Triggers: `snakemake rule download` or manual execution
- Responsibilities: Download config from Box/Zenodo and save to `metadata/ggo_config.yml`

## Error Handling

**Strategy:** Minimal error handling; relies on exceptions to fail fast and fail loudly.

**Patterns:**
- Config validation in `utils.load_panels()` raises `ValueError` if required keys are missing
- File existence checks in `scripts/myeloid_analysis.py` (`os.path.exists()`) to decide whether to recompute or load cached UMAP
- Backup URL fallback in `scanpy.read()` (via `backup_url` parameter) for remote data availability
- Test validation via `pytest` for utility functions and business logic

## Cross-Cutting Concerns

**Logging:**
- Framework: `scanpy.logging` and print statements
- Approach: Informational messages printed to stdout (e.g., "Reading {path}...", "Plot celltype heatmap")
- No structured logging; traces execution but does not persist logs

**Validation:**
- Config keys validated in `utils.load_panels()` and `utils.load_single_panel()`
- Panel names and AnnData keys checked against config structure
- Test cases validate critical functions in `test_units.py`

**Authentication:**
- No authentication required for local execution
- Zenodo/Box URLs accessed anonymously for backup data
- Snakemake `run_container()` handles container runtime selection (Docker/Apptainer)

---

*Architecture analysis: 2026-03-19*
