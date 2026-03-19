# Codebase Structure

**Analysis Date:** 2026-03-19

## Directory Layout

```
ggo-imc/
├── Snakefile                           # Main pipeline orchestration
├── config.yaml                         # Snakemake configuration (repo_root, project_rel)
├── requirements.txt                    # Python dependencies
├── pyproject.toml                      # Project metadata
├── README.md                           # Documentation and reproducibility guide
├── figures_checksums.md5               # Output checksum manifest for verification
├── .claude/                            # Claude agent configuration
│   └── settings.local.json             # Local settings
├── metadata/                           # Configuration and metadata
│   └── ggo_config.yml                  # (Downloaded at runtime) Panel paths, colors, markers
├── scripts/                            # Analysis and utility scripts
│   ├── utils.py                        # Shared utilities: config/data loading, figure saving
│   ├── celltype_heatmap_info.py        # Figure 1: cell type heatmap and UMAP plots
│   ├── roi_pca_plot.py                 # Figure 1: ROI PCA archetype overlay
│   ├── celltype_differential_abundance.py # Figures 2, 3: differential density analysis
│   ├── t_cell_analysis.py              # Figure 2: T cell functional states
│   ├── myeloid_analysis.py             # Figure 2: myeloid/macrophage polarization
│   ├── epithelial_characterization.py  # Figure 3: epithelial phenotypes and EMT
│   ├── ue_analysis.py                  # Figure 4: UTAG microenvironment analysis
│   ├── roi_pca_plot_group.py           # Figure 5: patient group PCA clustering
│   ├── patient_group.py                # Figure 5: patient risk stratification and scores
│   ├── download_yaml.py                # Utility: download config from cloud
│   ├── concat_anndata.py               # Utility: concatenate per-sample AnnData
│   ├── label_metadata_anndata.py       # Utility: attach clinical metadata
│   ├── generate_manifests.py           # Utility: generate batch manifests
│   ├── benchmark_integration.py        # Utility: compare integration methods
│   ├── benchmark_segmentation.py       # Utility: compare segmentation approaches
│   ├── exploratory/                    # Exploratory scripts (not in main pipeline)
│   └── scripts_backup/                 # Archived scripts (legacy implementations)
├── tests/                              # Unit tests
│   ├── __init__.py                     # Test package marker
│   ├── conftest.py                     # Pytest configuration and fixtures
│   ├── test_imports.py                 # Verify module imports
│   └── test_units.py                   # Unit tests: utils.py, patient_group.py
├── figures/                            # Generated output figures (per-figure subdirectories)
│   ├── figure1/                        # Cell type heatmaps and ROI PCA
│   │   ├── celltype/                   # Per-panel heatmaps and UMAPs
│   │   │   ├── PANEL_G/
│   │   │   └── PANEL_H/
│   │   └── roi_pca_all_combined.pdf
│   ├── figure2/                        # Immune cell density and functional analysis
│   │   ├── densities_lymphoid/         # T cell and B cell density
│   │   ├── densities_myeloid/          # Macrophage and myeloid density
│   │   └── (t cell, myeloid plots)
│   ├── figure3/                        # Stromal and epithelial analysis
│   │   ├── densities_stromal/          # Fibroblast density
│   │   ├── densities_epithelial/       # Epithelial cell density
│   │   ├── epithelial_proportion_*.pdf # Epithelial phenotype proportions
│   │   └── EMT proportion.pdf          # PanCK vs Vimentin density
│   ├── figure4/                        # Microenvironment analysis
│   │   ├── PANEL_G_uE_broad_*/         # Per-panel microenvironment composition
│   │   ├── PANEL_H_uE_broad_*/
│   │   ├── PANEL_G/uE_area/            # Area-based statistical tests
│   │   └── uE_broad/pathology/         # Functional marker dotplots
│   └── figure5/                        # Patient stratification
│       ├── group_process_scores.pdf
│       ├── cond_prob_*.pdf             # Conditional probability plots
│       └── (group-level analysis plots)
├── figures_baseline/                   # Reference baseline figures (for testing)
├── figures_backup/                     # Historical figure versions
├── data/                               # Local data directory
│   └── IMAGES/                         # IMC image files (organized by panel)
├── documentation/                      # Papers, meeting notes, submission materials
│   └── 251229 Cancer Cell Submission/
│       └── Figures/                    # Manuscript figures
└── .snakemake/                         # Snakemake cache (auto-generated)
```

## Directory Purposes

**scripts/:**
- Purpose: All analysis code, organized by figure and function type
- Contains: Figure-generation scripts, utility scripts, benchmarks, exploratory code
- Key files: `utils.py` (core), figure-specific scripts (`celltype_heatmap_info.py`, `roi_pca_plot.py`, etc.)

**metadata/:**
- Purpose: Project-level configuration and annotations
- Contains: `ggo_config.yml` (downloaded at runtime from Box)
- Key files: Config YAML with per-panel AnnData paths, marker lists, color schemes, condition labels

**tests/:**
- Purpose: Validation of core utility functions and pure functions
- Contains: Pytest unit tests for `utils.py` and `patient_group.py`
- Key files: `test_units.py` (main tests), `conftest.py` (fixtures)

**figures/:**
- Purpose: All generated publication figures, organized by figure number
- Contains: PDFs, PNGs, SVGs organized by figure/component/panel
- Key files: Linked to Snakemake `.done` sentinels (via rule outputs)

**data/:**
- Purpose: Local IMC image and raw data storage
- Contains: IMAGES directory with per-panel and per-sample image files
- Key files: Referenced in config YAML metadata

**documentation/:**
- Purpose: Manuscript materials, meeting notes, submission documents
- Contains: Latex/PDF versions of figures, journal submission materials
- Key files: Not part of pipeline; for reference only

## Key File Locations

**Entry Points:**
- `Snakefile`: Main orchestration entry point; invoked as `snakemake -d . -s Snakefile [target]`
- `scripts/celltype_heatmap_info.py`: Figure 1 entry (loads config, generates heatmap)
- `scripts/roi_pca_plot.py`: Figure 1 entry (ROI PCA analysis)

**Configuration:**
- `config.yaml`: Snakemake config (repo paths)
- `metadata/ggo_config.yml`: Project config (AnnData paths, metadata, colors) — downloaded at runtime

**Core Logic:**
- `scripts/utils.py`: Shared loading/saving functions; imported by all figure scripts
- `scripts/celltype_differential_abundance.py`: Shared density analysis; called for Figures 2 and 3
- `imc_analysis` (external package): Domain-specific tasks (heatmaps, density, statistical tests)

**Testing:**
- `tests/test_units.py`: Unit tests for `utils.ensure_dir()`, `utils.load_panels()`, `patient_group.cond_prob()`
- `tests/conftest.py`: Pytest fixtures and import configuration

## Naming Conventions

**Files:**
- Pipeline scripts: `{cell_type}_{analysis_type}.py` (e.g., `t_cell_analysis.py`, `epithelial_characterization.py`)
- Utility scripts: Descriptive verb-noun (e.g., `download_yaml.py`, `concat_anndata.py`)
- Test files: `test_{module}.py` (e.g., `test_units.py`)
- Config files: All lowercase with underscores (e.g., `ggo_config.yml`)

**Directories:**
- Figure output: `figures/figure{N}/` for each manuscript figure
- Component subdivision: `figures/figure{N}/{component_name}/{panel}/`
- Test directory: `tests/` at project root

**Functions in utils.py:**
- camelCase with `_` separators: `load_config()`, `load_panels()`, `load_single_panel()`, `save_figure()`, `ensure_dir()`
- Constants: SCREAMING_SNAKE_CASE: `CONFIG_PATH`, `PANELS`, `CYTOKINE`

**Variables in scripts:**
- Abbreviations: `adata` (AnnData object), `p` (panel), `cond` (condition), `fig` (matplotlib figure)
- Boolean prefixes: `is_`, e.g., no usage observed but pattern would be standard
- Iterate vars: Single letter for inner loops (e.g., `for i, ax in enumerate(axes)`)

## Where to Add New Code

**New Figure Analysis:**
1. Create script: `scripts/{domain}_{analysis}.py`
2. Import from utils: `from utils import load_config, load_panels, save_figure, ensure_dir`
3. Follow pattern:
   ```python
   if __name__ == '__main__':
       metadata = load_config()
       adata_dict = load_panels(metadata)
       # ... analysis ...
       save_figure('figures/figure{N}/{name}.pdf')
   ```
4. Create Snakemake rule in `Snakefile`:
   ```snakemake
   rule my_analysis:
       output: touch("my_analysis.done")
       input: "scripts/my_analysis.py"
       shell: run_container("scripts/my_analysis.py") + " && touch my_analysis.done"
   ```
5. Add to aggregation rule: `rule figure{N}: input: "my_analysis.done", ...`

**New Utility Function:**
1. Add to `scripts/utils.py` with docstring
2. Follow existing patterns (see `load_panels()` signature and validation)
3. Add unit test to `tests/test_units.py`
4. Import in figure scripts via `from utils import new_function`

**New Test:**
1. Add test class to `tests/test_units.py`
2. Use pytest fixtures from `conftest.py` (e.g., `tmp_path`)
3. Follow AAA pattern (Arrange, Act, Assert)
4. Run with: `pytest tests/test_units.py::TestClassName::test_method`

## Special Directories

**figures/:**
- Purpose: Generated manuscript outputs
- Generated: Yes (created by figure scripts via `save_figure()`)
- Committed: No (*.pdf, *.png ignored via .gitignore, except checksums file)
- Cleanup: `snakemake -d . -s Snakefile clean` removes `.done` sentinels (not figures themselves)

**data/:**
- Purpose: Local IMC raw data and images
- Generated: No (user-downloaded)
- Committed: No (large binary files)
- Structure: Referenced in `ggo_config.yml` under per-panel paths

**.snakemake/:**
- Purpose: Snakemake runtime cache and logs
- Generated: Yes (auto-created by Snakemake)
- Committed: No (.gitignore)
- Content: Rule dependency graphs, execution logs

**scripts_backup/:**
- Purpose: Legacy script versions (pre-refactor)
- Generated: No (manually archived)
- Committed: Yes (for historical reference)
- Status: Not used in current pipeline; can be deleted

**metadata/:**
- Purpose: Configuration storage (ggo_config.yml downloaded at runtime)
- Generated: Yes (downloaded from Box in `download` rule)
- Committed: No (ggo_config.yml has dynamic paths; .gitignore)
- Note: Directory itself committed; YAML file is generated

---

*Structure analysis: 2026-03-19*
