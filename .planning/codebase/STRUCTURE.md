# Codebase Structure

**Analysis Date:** 2026-03-19

## Directory Layout

```
projects/imc/ggo-imc/
├── scripts/                    # Analysis scripts (pipeline + utilities)
│   ├── utils.py               # Shared utilities (config, loading, saving)
│   ├── celltype_heatmap_info.py     # Figure 1: cell type heatmap
│   ├── roi_pca_plot.py        # Figure 1: ROI PCA (all pathologies)
│   ├── celltype_differential_abundance.py  # Figures 2, 3: cell density plots
│   ├── t_cell_analysis.py     # Figure 2: T cell functional analysis
│   ├── myeloid_analysis.py    # Figure 2: myeloid polarization
│   ├── epithelial_characterization.py  # Figure 3: epithelial phenotypes, EMT
│   ├── ue_analysis.py         # Figure 4: microenvironments, TLS, interactions
│   ├── roi_pca_plot_group.py  # Figure 5: PCA by patient risk group
│   ├── patient_group.py       # Figure 5: patient risk stratification
│   ├── download_yaml.py       # Utility: download config from Zenodo/Box
│   ├── exploratory/           # Non-pipeline exploratory scripts
│   └── scripts_backup/        # Previous versions (archived)
├── tests/                      # Test suite (unit tests for utils, patient_group)
│   ├── __init__.py
│   ├── conftest.py            # pytest configuration
│   ├── test_units.py          # Unit tests for utility functions
│   └── test_imports.py        # Import validation
├── Snakefile                   # Workflow orchestration (12 rules)
├── config.yaml                 # Snakemake config (repo_root, project_rel)
├── metadata/                   # Runtime metadata and configuration
│   ├── ggo_config.yml         # Main config (downloaded at runtime)
│   ├── backup/                # Cached metadata
│   ├── interactions/          # Interaction analysis results
│   └── roi_pca_anndata/       # PCA density matrices
├── figures/                    # Generated output figures (main pipeline output)
│   ├── figure1/               # Celltype heatmap, UMAP, ROI PCA
│   ├── figure2/               # Immune dynamics (lymphoid, myeloid)
│   ├── figure3/               # Stromal expansion, epithelial, EMT
│   ├── figure4/               # Microenvironments, TLS, interactions
│   ├── figure5/               # Patient risk groups, diagnostic gap
│   ├── benchmark/             # Integration/segmentation benchmarks
│   ├── demographics/          # Cohort overview figures
│   ├── densities_*/           # Cell density outputs by condition
│   ├── QC/                    # Quality control plots
│   └── Figures/               # Manuscript composite figures
├── figures_baseline/          # Reference output (for reproducibility checks)
├── figures_backup/            # Previous generation outputs (archived)
├── processed/                 # Intermediate analysis results
│   ├── filtered/              # Filtered AnnData objects
│   ├── quantification/        # Quantified cell-type data
│   ├── PANEL_G/               # Panel-specific processed data
│   ├── PANEL_H/               # Panel-specific processed data
│   └── roi_pca_anndata/       # ROI-level density matrices
├── results/                   # Statistical results and outputs
│   ├── patient_csvs/          # Patient-level stratification results
│   ├── phenotyping/           # Cell phenotyping outputs
│   └── previous/              # Historical results
├── data/                      # Raw input data
│   └── IMAGES/                # Raw IMC image files
├── images/                    # Segmented/labeled images
│   ├── labeled/               # Cell segmentation masks
│   └── stack/                 # Image stacks
├── .planning/                 # GSD planning directory
│   └── codebase/              # This analysis (ARCHITECTURE.md, STRUCTURE.md)
├── .claude/                   # Claude agent profiles
├── .snakemake/                # Snakemake runtime metadata
├── documentation/             # Manuscript drafts, meeting notes
├── previous_code/             # Legacy code (archived)
├── README.md                  # Repository overview
├── Plan.md                    # Project plan and phase tracking
├── pyproject.toml             # Python project metadata
├── requirements.txt           # Package dependencies
├── Snakefile                  # Main workflow definition
├── figures_checksums.md5      # Reproducibility verification manifest
├── input_data_checksums.md5   # Input data integrity check
└── *.done                     # Sentinel files (rule completion markers)
```

## Directory Purposes

**`scripts/`:**
- Purpose: All Python analysis code for the pipeline
- Contains: 9 figure-generation scripts, utility functions, download script, archived versions
- Key files: `utils.py` (shared layer), `celltype_*.py`, `*_analysis.py`, `roi_pca_*.py`, `patient_*.py`
- Output: Figure files (PDF/PNG) saved to `figures/` by each script

**`tests/`:**
- Purpose: Unit test suite for utility functions and business logic
- Contains: pytest tests for `utils.py` and `patient_group.py`
- Key files: `test_units.py` (TestEnsureDir, TestLoadPanelsValidation, TestCondProb)
- Excludes: Integration tests, E2E tests (not present in codebase)

**`metadata/`:**
- Purpose: Store and manage configuration and intermediate analysis artifacts
- Contains: `ggo_config.yml` (main config, downloaded at runtime), interaction matrices, PCA density outputs
- Key files: `ggo_config.yml`, `roi_pca_anndata/`, `interactions/`
- Not committed: Raw config (downloaded on first run)

**`figures/`:**
- Purpose: Output directory for all generated figures
- Contains: Figure outputs organized by figure number (figure1–5), plus benchmarks, demographics, QC
- Key structure: `figures/figure{1-5}/` with subdirectories per analysis (e.g., `figure2/densities_immune/`, `figure4/uE_broad/`)
- Output by: Each analysis script via `save_figure()`

**`figures_baseline/`:**
- Purpose: Baseline reference figures for reproducibility verification
- Contains: Known-good outputs from verified pipeline run
- Usage: Compared against generated `figures/` using checksum manifest (`figures_checksums.md5`)

**`processed/`:**
- Purpose: Store intermediate analysis results and cached computations
- Contains: Filtered AnnData objects, quantification matrices, panel-specific data, ROI-level density matrices
- Examples: `processed/quantification/`, `processed/filtered/`, `processed/roi_pca_anndata/`
- Caching: Some scripts check for precomputed results (e.g., `myeloid_analysis.py` checks for cached UMAP before recomputation)

**`results/`:**
- Purpose: Final statistical results and patient-level stratification
- Contains: Patient-level CSV outputs, phenotyping results, historical results
- Examples: `results/patient_csvs/`, `results/phenotyping/`

**`data/` and `images/`:**
- Purpose: Input data (raw images, segmentation masks)
- Contains: Raw IMC image stacks, labeled segmentation masks, image quantification files
- Size: Large (raw images are multi-channel TIFF stacks)
- Committed: No (tracked separately, available on Zenodo)

**`.snakemake/`:**
- Purpose: Snakemake runtime metadata and caching
- Contains: Logs, metadata, shadows, conda environments, lock files
- Not committed: Generated at runtime

**`.planning/`:**
- Purpose: GSD (Claude orchestrator) planning documentation
- Contains: Phase tracking, execution plans, codebase analysis (this directory)
- Key files: `Plan.md` (updated by documentor), `codebase/ARCHITECTURE.md`, `codebase/STRUCTURE.md`

## Key File Locations

**Entry Points:**
- `Snakefile`: Workflow orchestration, defines all rules and targets
- `scripts/*.py`: Individual analysis scripts (9 pipeline scripts)
- `scripts/download_yaml.py`: First step in pipeline (downloads config)

**Configuration:**
- `config.yaml`: Snakemake configuration (repo_root, project_rel, container_sif)
- `metadata/ggo_config.yml`: Runtime metadata (downloaded from Box/Zenodo)

**Core Logic:**
- `scripts/utils.py`: Shared utilities (load_config, load_panels, save_figure, ensure_dir)
- `scripts/celltype_heatmap_info.py`: Cell type heatmap and UMAP for Figure 1
- `scripts/roi_pca_plot.py`: ROI PCA archetype plot (Figure 1)
- `scripts/celltype_differential_abundance.py`: Density comparison across conditions (Figures 2, 3)
- `scripts/t_cell_analysis.py`: Lymphocyte functional state analysis (Figure 2)
- `scripts/myeloid_analysis.py`: Myeloid cell polarization (Figure 2)
- `scripts/epithelial_characterization.py`: Epithelial subtypes and EMT (Figure 3)
- `scripts/ue_analysis.py`: Microenvironment composition and interactions (Figure 4)
- `scripts/roi_pca_plot_group.py`: PCA colored by patient risk group (Figure 5)
- `scripts/patient_group.py`: Patient risk stratification (Figure 5)

**Testing:**
- `tests/test_units.py`: Unit tests for utils.py and patient_group.py
- `tests/conftest.py`: pytest fixtures and configuration

**Verification:**
- `figures_checksums.md5`: MD5 hashes of baseline figures for reproducibility
- `input_data_checksums.md5`: MD5 hashes of input data files

## Naming Conventions

**Files:**
- Pipeline scripts: `{analysis}_{target}.py` (e.g., `celltype_heatmap_info.py`, `roi_pca_plot.py`)
- Utility scripts: `{purpose}.py` (e.g., `utils.py`, `download_yaml.py`)
- Test files: `test_*.py` (e.g., `test_units.py`)
- Config: YAML files in `metadata/` (e.g., `ggo_config.yml`)
- Output figures: Named by content (e.g., `celltype_heatmap.pdf`, `roi_pca_all_combined.pdf`)

**Directories:**
- Figure outputs: `figures/figure{1-5}/` (by figure number)
- Analysis outputs: `figures/{figure}/densities_{lineage}/{panel}/` (e.g., `figures/figure2/densities_immune/PANEL_G/`)
- Panels: Named `PANEL_G`, `PANEL_H` (marker sets)
- Conditions: `pathology`, `radio` (pathology and radiology), `Group` (patient risk group)

**Variables:**
- Panel iterator: `p` (e.g., `for p in ['PANEL_G', 'PANEL_H']`)
- Metadata dict: `metadata` (from `load_config()`)
- AnnData objects: `adata`, `adata_dict` (panel-keyed dict)
- Density matrices: `density`, `*_density` (output of `imc.tl.celltype_density()`)
- Cell lineage categories: `celltype`, `celltype_broad` (from config or AnnData.obs)

## Where to Add New Code

**New Figure/Analysis:**
- Primary code: `scripts/new_analysis.py` (import from `utils.py`, follow pattern of existing scripts)
- Add rule to `Snakefile` (define input/output, call `run_container()`)
- Add to appropriate figure aggregate rule (e.g., `rule figure6: input: ...`)
- Output figures to `figures/figure6/` directory
- Test: Add unit tests to `tests/test_units.py` if logic is extracted as a pure function

**New Utility Function:**
- Implementation: `scripts/utils.py` (add function with docstring)
- Pattern: Functions should be config-agnostic (take parameters explicitly) or use module-level constants (PANELS, CYTOKINE)
- Export: No explicit export list; all functions available to importing scripts
- Testing: Add test class to `tests/test_units.py` if function is testable without AnnData

**New Configuration Key:**
- Location: `metadata/ggo_config.yml` (YAML format)
- Access: Via `metadata[key]` or `metadata[panel]['AnnData'][key]`
- Pattern: Follow existing structure (panel-level dicts, AnnData section, color palettes)

**New Cell Type/Marker Set:**
- Location: Define in analysis script or `metadata/ggo_config.yml`
- Pattern: Existing scripts inline lists (e.g., `pro_inflammatory_markers`) or load from config
- Usage: Pass to `scanpy` functions or `imc_analysis` functions

**New Benchmark/Exploratory Script:**
- Location: `scripts/benchmark_*.py` (not included in pipeline) or `scripts/exploratory/`
- Pattern: Can import from `utils.py` but not required to follow pipeline discipline
- Execution: Run manually, not via Snakemake rule

## Special Directories

**`scripts_backup/`:**
- Purpose: Archive of previous script versions
- Generated: No (manually maintained)
- Committed: Yes (serves as historical record)

**`.snakemake/`:**
- Purpose: Snakemake internal state
- Generated: Yes (at runtime)
- Committed: No (in .gitignore)

**`figures_baseline/`:**
- Purpose: Reference baseline figures for reproducibility
- Generated: Yes (committed as baseline after verification)
- Committed: Yes (serves as gold standard)

**`.pytest_cache/`:**
- Purpose: pytest internal caching
- Generated: Yes (at runtime)
- Committed: No (in .gitignore)

**`.claude/`:**
- Purpose: Claude agent profiles and instructions
- Generated: No (manually maintained by orchestrator)
- Committed: Yes (shared context for agents)

---

*Structure analysis: 2026-03-19*
