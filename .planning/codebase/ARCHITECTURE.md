# Architecture

**Analysis Date:** 2026-03-19

## Pattern Overview

**Overall:** Snakemake-orchestrated data analysis pipeline with compartmentalized figure-generation scripts

**Key Characteristics:**
- Figure-centric design: each manuscript figure maps 1:1 to analysis phases
- Modular scripts: each script loads data independently and produces figure outputs
- Shared utility layer: common loading, saving, and config patterns extracted to `scripts/utils.py`
- Configuration-driven: single YAML config file (`metadata/ggo_config.yml`) provides data paths and metadata
- Container-agnostic: runtime abstraction via `run_container.sh` supports Docker (macOS), Apptainer (Linux/HPC), or local conda

## Layers

**Pipeline Orchestration:**
- Purpose: Coordinate execution order and sentinel-based dependency tracking
- Location: `Snakefile`
- Contains: Figure-level rules (e.g., `figure1`, `figure2`), component rules (e.g., `celltype`, `t_cell`), and aggregation targets
- Depends on: Shell scripts (`run_container.sh`)
- Used by: Snakemake CLI invoked with targets like `snakemake -d . -s Snakefile figure1`

**Configuration Layer:**
- Purpose: Centralize project metadata, file paths, colors, and condition labels
- Location: `metadata/ggo_config.yml` (downloaded at runtime)
- Contains: Per-panel AnnData paths, color mappings, marker gene lists, clinical condition definitions
- Depends on: Zenodo/Box cloud storage (fallback URLs)
- Used by: All scripts via `load_config()` from `utils.py`

**Shared Utilities:**
- Purpose: Eliminate boilerplate for config, data loading, and figure saving
- Location: `scripts/utils.py`
- Contains: `load_config()`, `load_panels()`, `load_single_panel()`, `save_figure()`, `ensure_dir()`
- Depends on: scanpy, imc_analysis (external packages)
- Used by: All 9 figure-generation scripts

**Figure-Generation Scripts:**
- Purpose: Load data, apply analysis, generate and save publication-ready figures
- Location: `scripts/{celltype,roi_pca,t_cell,myeloid,epithelial,ue,patient}_*.py`
- Contains: Domain-specific analysis (clustering, density computation, statistical testing) and matplotlib plotting
- Depends on: Shared utilities, imc_analysis library, scanpy, pandas, numpy
- Used by: Snakemake rules via container runtime

**Utility Scripts:**
- Purpose: Data preparation and utility functions
- Location: `scripts/{download_yaml,concat_anndata,label_metadata,generate_manifests,benchmark_*}.py`
- Contains: Configuration download, data concatenation, metadata attachment, manifest generation, method benchmarking
- Depends on: imc_analysis, scanpy
- Used by: Manual execution or future pipeline phases

**Tests:**
- Purpose: Validate utility functions and pure functions
- Location: `tests/test_*.py`
- Contains: Unit tests for `utils.py` functions and `patient_group.py` conditional probability logic
- Depends on: pytest, pandas (test fixtures)
- Used by: Developers for regression testing

## Data Flow

**Main Pipeline (Figure 1-5 Generation):**

1. **Configuration Setup**
   - User invokes: `snakemake -d . -s Snakefile figure1`
   - Snakemake triggers `download` rule (if not cached)
   - `download_yaml.py` downloads `metadata/ggo_config.yml` from Box URL
   - Config is parsed as dict and injected into each script

2. **Panel Loading**
   - Scripts call `load_panels(metadata)` or `load_single_panel(metadata, 'PANEL_G')`
   - Utility checks `metadata['PANEL_*']['AnnData']['phenotyped_umap_name']` for local path
   - If missing, falls back to `metadata['PANEL_*']['AnnData']['backup_url']` (Zenodo)
   - AnnData objects returned as dict `{'PANEL_G': adata_g, 'PANEL_H': adata_h}`

3. **Analysis Execution**
   - Panel-specific analysis (e.g., `celltype_heatmap_info.py` clusters and plots)
   - Uses `imc_analysis` library for domain-specific tasks (heatmaps, density, statistical tests)
   - Scanpy/pandas for data manipulation
   - Matplotlib for visualization

4. **Output Generation**
   - All figures saved via `save_figure(path)` from utils
   - Path format: `figures/figure{N}/<component>/<panel>/<plot_type>.pdf`
   - Directory creation handled automatically via `os.makedirs(..., exist_ok=True)`
   - Snakemake creates `.done` sentinel files for dependency tracking

**Data Structure (AnnData):**
```
adata = AnnData(X=gene_expression_matrix, obs=cell_metadata, var=gene_metadata, obsm=embeddings, uns=colors_and_metadata)
├── .X              # Gene expression matrix (cells × markers)
├── .obs            # Cell annotations: celltype, celltype_broad, condition labels, spatial coords
├── .var            # Marker metadata: gene names, marker families
├── .obsm           # Embeddings: X_umap, X_pca (spatial projections)
└── .uns            # Metadata: celltype_colors, radio_colors, pathology_colors, marker groups
```

**State Management:**
- Immutable input: AnnData objects loaded from disk (copy-on-modify pattern used sparingly)
- Analysis-scoped state: temporary AnnData subsets within script lifetime
- Output state: figures written to disk, checksums tracked in `figures_checksums.md5`
- No persistent session state between scripts (idempotent execution)

## Key Abstractions

**Figure Target:**
- Purpose: Represent a publication figure as a collection of dependent analysis components
- Examples: `figure1` = `celltype` + `pca` rules
- Pattern: Snakemake aggregation rule with list of `.done` sentinel files

**Panel:**
- Purpose: Represent a distinct imaging mass cytometry panel with separate markers and cells
- Examples: `PANEL_G`, `PANEL_H`
- Pattern: String key in config; AnnData loaded per-panel to allow independent analysis

**Condition:**
- Purpose: Represent a clinical stratification axis (e.g., histology, radiology diagnosis)
- Examples: `'pathology'`, `'radio'` as obs columns
- Pattern: Categorical metadata in AnnData; used for grouping and statistical testing

**Microenvironment (uE):**
- Purpose: Represent spatial organization of cells (UTAG annotation)
- Examples: `'uE_broad'` in obs (used in Figure 4)
- Pattern: Categorical clustering in AnnData; accessed via obs groupby operations

## Entry Points

**Snakemake Orchestration:**
- Location: `Snakefile` + `run_container.sh`
- Triggers: `snakemake -d . -s Snakefile [target]`
- Responsibilities:
  - Resolves dependencies (celltype → pca → figure1)
  - Runs scripts inside container (Docker/Apptainer)
  - Tracks completion via `.done` sentinel files
  - Supports manual script execution with `SC_TOOLS_RUNTIME=none` env var

**Figure Generation Scripts:**
- Location: `scripts/celltype_heatmap_info.py`, `scripts/roi_pca_plot.py`, etc.
- Triggers: Invoked by Snakemake rules via `run_container("scripts/script_name.py")`
- Responsibilities:
  - Load config and panels
  - Perform domain-specific analysis
  - Generate and save figures
  - Exit with code 0 on success (Snakemake creates `.done` file)

**Direct Execution (Development):**
- Scripts can be run standalone: `python scripts/celltype_heatmap_info.py`
- Requires: Working directory set to repo root, conda environment activated
- Pattern: `if __name__ == '__main__':` guards all analysis code

## Error Handling

**Strategy:** Fail-fast with clear error messages

**Patterns:**

**Config Validation in utils.py:**
```python
# load_panels() validates presence of required keys before loading
if p not in metadata:
    raise ValueError(f"Panel '{p}' not found in config. Available keys: {sorted(metadata.keys())}")
if anndata_section is None:
    raise ValueError(f"Config section metadata['{p}']['AnnData'] is missing.")
```

**AnnData Backup (transparent fallback):**
```python
# scanpy.read() supports backup_url parameter for cloud fallback
adata = sc.read(path, backup_url=url)
```

**Container Exit Codes:**
- Snakemake interprets non-zero exit codes as rule failures
- `.done` sentinels only created on successful exit
- Failed rule prevents downstream rules from executing

## Cross-Cutting Concerns

**Logging:**
- Approach: `scanpy.logging` and `print()` statements for progress tracking
- Pattern: `sc.logging.info()` used in utils.py; scripts use `print()` for task descriptions

**Validation:**
- Approach: Inline assertions and explicit error messages
- Pattern: Config validation in `load_panels()`; figure output paths validated via `ensure_dir()`

**Authentication:**
- Approach: Cloud URLs hardcoded in config; no credentials in code
- Pattern: Box URL for config, Zenodo URLs for AnnData backups (public access)

**Reproducibility:**
- Approach: Deterministic random seeds and checksum verification
- Pattern: `random_state=0` in subsampling (Figure 1); `figures_checksums.md5` manifest for output verification

---

*Architecture analysis: 2026-03-19*
