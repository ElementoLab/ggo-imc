# Codebase Concerns

**Analysis Date:** 2026-03-19

## Missing Critical Files

**Figure 5 R script not committed:**
- Issue: The Snakefile references `scripts/asd.R` (lines 73-78), but this file is not present in the repository
- Files: `Snakefile`, `README.md` (line 63)
- Impact: The `patient` rule cannot run without this file; Figure 5 cannot be fully reproduced from the repo. The `all` target explicitly excludes `patient.done` (line 111)
- Fix approach: Either commit `scripts/asd.R` with R/bioconductor dependencies documented, or remove the `patient` rule and update documentation to clarify reproducibility limitations

## Large Backup Directories

**Deprecated code and figures not cleaned up:**
- Issue: `scripts_backup/` (29M) and `figures_backup/` (6.3G) contain old versions taking significant disk space
- Files: `scripts_backup/`, `figures_backup/`
- Impact: Repository size bloated; clutter when browsing codebase; no clear purpose for these backups
- Fix approach: Move to external archive or remove entirely. `.gitignore` (lines 3-4) already excludes them from new commits, but historical commits still contain them. Consider `git filter-branch` or `git-filter-repo` if backups should be removed from history

## Incomplete Test Coverage

**Pipeline scripts untested:**
- Issue: Only `utils.py` and `patient_group.py` have unit tests. Nine pipeline scripts (celltype_heatmap_info.py, roi_pca_plot.py, t_cell_analysis.py, myeloid_analysis.py, epithelial_characterization.py, ue_analysis.py, roi_pca_plot_group.py, celltype_differential_abundance.py, benchmark_*.py) have no automated tests
- Files: `tests/test_units.py` (81 lines total; only tests utils and patient_group), all figure-generating scripts in `scripts/`
- Impact: Figure-generating scripts can silently fail or produce incorrect figures. Regressions in shared utilities (via utils.py) could propagate across all figures. Manual verification via checksum comparison (README.md lines 83-97) is the only safety net
- Fix approach: Add integration tests that verify pipeline outputs match baseline checksums. Add unit tests for functions in figure scripts (e.g., plotting logic in ue_analysis.py, t_cell_analysis.py). Use pytest markers to distinguish unit vs. integration tests

## Hardcoded Metadata Keys

**Fragile coupling to config structure:**
- Issue: Scripts reference specific nested dictionary keys from `metadata` dict (e.g., `metadata['PANEL_H']['AnnData']['lymphocytes_url']`) without validation. Different scripts use different key names and structures
- Files:
  - `scripts/t_cell_analysis.py` (lines 23, 46): uses `lymphocytes_url`
  - `scripts/myeloid_analysis.py` (lines 46, 74): uses `myeloids_url`
  - `scripts/patient_group.py` (line 51): uses `patient_group_url`
  - `scripts/roi_pca_plot.py` (line 18): uses `PCA_URL`
  - `scripts/roi_pca_plot_group.py` (lines 21, 32): uses `PCA_URL`, `patient_group_url`
  - `scripts/ue_analysis.py` (line 19): uses custom `utag_url` key
- Impact: Config schema is implicit and undocumented. Adding a new panel or changing metadata structure requires manual updates to multiple scripts. Easy to miss and cause runtime KeyError exceptions
- Fix approach: Create a schema validation module (e.g., `scripts/config_schema.py`) that defines required keys with defaults. Use `utils.load_config()` to validate at load time. Document the config structure in README.md or a separate SCHEMA.md

## Bare Exception Handling

**Generic exception catching masks real errors:**
- Issue: Multiple scripts catch generic `Exception` without logging useful context
- Files:
  - `scripts/benchmark_segmentation.py` (lines 107-108, 130-131, 202-203): `except Exception: logger.exception()` but proceeds silently
  - `scripts/benchmark_integration.py`: assumes backup_url fallback exists without validation
- Impact: If segmentation fails for a reason other than missing imports, the error is logged but execution continues with empty results. May produce incomplete or misleading benchmark reports
- Fix approach: Catch specific exceptions (e.g., `ImportError`, `RuntimeError`) and re-raise or fail gracefully. For benchmark.py, validate that comparison results are complete before aggregating (check line 207-209 which only logs "No results" but doesn't exit early enough)

## Sentinel Files for Workflow

**Fragile state tracking:**
- Issue: Snakefile uses `.done` files as outputs (e.g., `celltype.done`, `pca.done`) instead of real figure outputs
- Files: `Snakefile` (lines 19, 25, etc.), all `.done` files in root
- Impact: If a figure is partially generated (e.g., PDF written but metadata incomplete), the `.done` file still marks it complete. Running `snakemake clean` removes all sentinel files, forcing full rerun even if some figures are valid. No checksum validation during pipeline
- Fix approach: Use actual figure outputs as Snakemake targets (e.g., `figures/figure1/celltype_heatmap.pdf` instead of `celltype.done`). Optionally compute and validate checksums after each figure generation

## Missing Error Handling in Data Loading

**Silent failures on missing backup URLs:**
- Issue: `utils.py` `load_panels()` function validates presence of metadata keys (lines 96-112) but doesn't validate that backup URLs are reachable or well-formed
- Files: `scripts/utils.py` (lines 62-121)
- Impact: If a backup URL is malformed or the service is down, `sc.read(..., backup_url=url)` (line 119) will attempt to download and fail at read time with a cryptic error. The validation only checks key existence, not value validity
- Fix approach: Add optional URL validation in `load_config()` or before calling `sc.read()`. Document expected URL format and error recovery strategy

## Unsupported Python/Dependency Versions

**Broad version constraints:**
- Issue: `pyproject.toml` (line 9) specifies `requires-python = ">=3.10"` but README.md (line 26) shows Python 3.9 in conda setup. Dependencies in `requirements.txt` have no upper bounds
- Files: `pyproject.toml`, `requirements.txt`, `README.md`
- Impact: Unclear which Python version is actually supported. Dependency updates could introduce breaking changes (e.g., scimap or squidpy API changes). IPython pinned to 8.18.1 (requirements.txt line 4) for unclear reason
- Fix approach: Test and document actual supported Python versions (3.9 vs 3.10+). Pin major versions of key dependencies (scanpy, anndata, scipy) with upper bounds or use a lock file (e.g., uv.lock or pip-compile). Document why IPython 8.18.1 is required

## Incomplete Data Checksums

**Baseline checksums may be stale:**
- Issue: `figures_checksums.md5` contains checksums for verified outputs, but only covers PDF/PNG/SVG files. Binary outputs (.h5ad, intermediate data) are not verified
- Files: `figures_checksums.md5`, README.md (lines 81-97)
- Impact: Intermediate data corruption (e.g., in `results/` or processed files) could silently propagate to figures that still match checksums because they're derived from corrupted input
- Fix approach: Extend checksum manifest to include key data files. Compute checksums at each pipeline stage and verify during the Snakefile run

## Unused or Archived Code Not Cleaned

**Exploratory scripts mixed with pipeline:**
- Issue: `scripts/exploratory/` contains 20+ exploratory/development scripts (not called by Snakefile) alongside the 9 main pipeline scripts
- Files: `scripts/exploratory/` directory (celltype_interaction.py, cell_proportion_across_condition.py, etc.)
- Impact: When new developers review code, it's unclear which scripts are "canonical" vs. experimental. Maintenance burden if exploratory code has bugs or outdated dependencies
- Fix approach: Move exploratory scripts to a separate `scripts_exploratory/` directory (outside main pipeline). Document in README which scripts are experimental. Alternatively, remove if no longer needed

## Config Download May Fail Silently

**Network request with no retry logic:**
- Issue: `utils.load_config()` calls `imc.utils.download()` (line 54) which may fail if network is unavailable. No retry logic or timeout specification
- Files: `scripts/utils.py` (lines 33-55)
- Impact: If Zenodo or Box is temporarily unavailable, the entire pipeline halts. No fallback to cached config or graceful degradation
- Fix approach: Add retry logic with exponential backoff. Cache the config locally after first successful download. Allow override via environment variable (e.g., `GGO_CONFIG_PATH`)

## Type Hints Missing in Key Functions

**Runtime type validation gaps:**
- Issue: Most figure-generating scripts lack type hints or parameter validation
- Files: `scripts/t_cell_analysis.py`, `scripts/epithelial_characterization.py`, `scripts/myeloid_analysis.py`, `scripts/ue_analysis.py` (all lack comprehensive type hints)
- Impact: Passing wrong AnnData subset or forgetting to load a panel doesn't fail until deep in the plotting code. Static type checkers (mypy) would catch this before runtime
- Fix approach: Add type hints to all script functions (especially those taking adata, metadata, fig, ax parameters). Run `mypy` in CI

## Hardcoded Colors and Plot Parameters

**Magic values scattered across scripts:**
- Issue: Plot colors, sizes, and figure dimensions are hardcoded in multiple scripts
- Files:
  - `scripts/patient_group.py` (line 83): `adata.uns['Group_colors'] = ['#99B898', '#FECEAB', '#E84A5F', '#2A363B']`
  - `scripts/t_cell_analysis.py`, `scripts/epithelial_characterization.py`: figure dims like `figsize=(12, 5)` scattered throughout
- Impact: Updating a color scheme or figure size requires searching and editing multiple files. Inconsistent styling across figures
- Fix approach: Centralize plot config in `metadata/ggo_config.yml` or a new `scripts/plot_config.py` module with color palettes, fonts, sizes

## No Logging Configuration in Pipeline Scripts

**Difficult to debug pipeline runs:**
- Issue: Most pipeline scripts have no logging setup. Only benchmark scripts call `logging.basicConfig()` (benchmark_segmentation.py line 32)
- Files: `scripts/celltype_heatmap_info.py`, `scripts/roi_pca_plot.py`, `scripts/t_cell_analysis.py`, `scripts/epithelial_characterization.py`, `scripts/ue_analysis.py`
- Impact: Progress is silent until figures are written. If a script hangs, no indication of where. Debugging requires adding print statements manually
- Fix approach: Add logging module setup to utils.py. Provide `--verbose` flag to scripts for optional debug logging

## Memory Management in Benchmark Scripts

**Potential memory leaks with large image files:**
- Issue: `benchmark_segmentation.py` loads large TIFF images and deletes references (lines 106, 129, 205) but doesn't explicitly call `gc.collect()` between methods until line 112
- Files: `scripts/benchmark_segmentation.py` (lines 94-256)
- Impact: On systems with limited RAM, processing 20 ROIs with multiple segmentation methods could OOM before garbage collection. Cellpose and StarDist both allocate GPU/CPU memory
- Fix approach: Add `gc.collect()` and explicit memory cleanup between major pipeline stages. Profile memory usage during benchmark

## Documentation Drift Risk

**README outdated elements:**
- Issue: README.md (line 62) mentions `snakemake figure5` doesn't include patient rule, and notes that R script is "not yet committed" (line 63). This is accurate but highlights fragility
- Files: `README.md`
- Impact: As codebase evolves, the README can fall out of sync. Currently there's no automation to verify the example commands work
- Fix approach: Add a CI check that runs the Snakemake examples in the README (at least with `--dryrun`). Update README to document the missing asd.R and when it will be available

---

*Concerns audit: 2026-03-19*
