# Codebase Concerns

**Analysis Date:** 2026-03-19

## Tech Debt

**Deprecated scripts not removed:**
- Issue: Large collection of backup scripts in `scripts_backup/` (40+ files) creating maintenance burden and confusion about which code is active
- Files: `scripts_backup/` directory (e.g., `scripts_backup/celltype_differential_abundance_backup.py`, `scripts_backup/patient_group.py`, `scripts_backup/t_cell_analysis.py`)
- Impact: Developers must maintain both old and new versions; unclear which version is the source of truth; increases repository size and cognitive load during debugging
- Fix approach: Audit backup scripts to confirm all functionality exists in active versions, then delete the entire `scripts_backup/` directory. Use git history if restoration is needed.

**Exploratory scripts not integrated into pipeline:**
- Issue: 27+ exploratory scripts in `scripts/exploratory/` directory that are not part of the main pipeline and lack documentation on their purpose or status
- Files: `scripts/exploratory/` (e.g., `rapids_scanpy_funcs.py` - 924 lines of RAPIDS GPU code, `cell_proportion_across_condition.py` - 1046 lines, `patient_feature_selection.py` - 676 lines)
- Impact: Unclear which exploratory work is experimental vs. production-ready; large unmaintained code paths; potential for stale dependencies (e.g., RAPIDS/CUDA libraries); confusion about which analysis should be run
- Fix approach: Move exploratory scripts to a separate `archive/` directory and document their purpose. Document which exploratory findings fed into the main pipeline figures.

**Snakemake directory in active repo:**
- Issue: `.snakemake/` directory is uncommitted (appears in git status as `??`) suggesting pipeline management is not fully integrated
- Files: `.snakemake/` directory
- Impact: Pipeline reproducibility unclear; intermediate files and cache state not under version control; difficult to share pipeline state across environments
- Fix approach: Add `.snakemake/` to `.gitignore` if it's a build artifact. Document the pipeline execution path (manual scripts vs. Snakemake-managed).

**Baseline figures directory uncommitted:**
- Issue: `figures_baseline/` directory present but not tracked in git (appears in git status as `??`)
- Files: `figures_baseline/` directory
- Impact: Unclear if this is reference output for regression testing; difficult to audit changes to figures
- Fix approach: Clarify purpose (baseline for comparison? previous results?). If it's for testing, add to `.gitignore` and document how to regenerate. If it's reference data, commit it.

## Missing Error Handling & Validation

**Unvalidated config keys in scripts:**
- Issue: Multiple scripts assume config keys exist without checking. If `load_config()` or YAML parsing fails, runtime errors occur downstream
- Files:
  - `scripts/utils.py:97-112` - Validates config structure but only for loaded panels
  - `scripts/t_cell_analysis.py:22-24` - Assumes `metadata['PANEL_H']['AnnData']['lymphocytes']` exists
  - `scripts/myeloid_analysis.py:14,43,71` - Assumes multiple metadata keys without validation
- Impact: Script crashes with unhelpful KeyError if YAML is malformed or missing sections; difficult to debug config issues in large projects
- Fix approach: Add a `validate_config()` function to `utils.py` that checks all required keys at startup. Return clear error messages listing missing keys.

**No validation of data shape/content in analysis scripts:**
- Issue: Scripts call `.isin()`, `.loc[]`, and column access without checking if required columns exist
- Files:
  - `scripts/epithelial_characterization.py:26-38` - Assumes `celltype` column exists; will fail silently if empty
  - `scripts/epithelial_characterization.py:79` - Accesses `pg_epi.obs[feature]` but doesn't check if `feature` is categorical
  - `scripts/myeloid_analysis.py:30-32` - Assumes specific celltype values exist
  - `scripts/patient_group.py:76` - Calls `.str.replace()` on `Group` obs column without checking if it exists
- Impact: Cryptic errors when input AnnData objects lack expected columns; difficult to debug data preparation issues
- Fix approach: Add assertions or try-except blocks with informative error messages checking for required obs/var columns at script start.

**Hardcoded marker lists not validated:**
- Issue: Marker/cytokine lists hardcoded and assumed to exist in all datasets
- Files:
  - `scripts/utils.py:26` - `CYTOKINE` list hardcoded
  - `scripts/t_cell_analysis.py:50,74` - `pro_inflammatory_markers`, `anti_inflammatory_markers` hardcoded
  - `scripts/myeloid_analysis.py:24-28, 60-64` - `h_myeloid_markers`, `g_myeloid_markers` hardcoded
- Impact: Script fails if a marker is missing from input data; no graceful degradation
- Fix approach: Check that all markers exist before subsetting: `missing = set(markers) - set(adata.var_names)`. Log warning and skip missing markers or fail with clear message.

## Known Bugs

**Potential edge case in epithelial_characterization.py:**
- Symptoms: Line 79 accesses `pg_epi.obs[feature]` which returns a categorical. If a category has no data, the loop on line 77 still runs, potentially creating empty plots
- Files: `scripts/epithelial_characterization.py:75-85`
- Trigger: Run script with pathology or radio condition that has no representation in `pg_epi` subset
- Workaround: Filter out empty categories before loop: `for ax in axes.flatten()` loop checks `pg_tmp = pg_epi[pg_epi.obs[feature] == rad].copy()` but doesn't validate `len(pg_tmp) > 0`

**Embedding alignment issue in benchmark_integration.py:**
- Symptoms: When multiple methods succeed but produce different cell orderings, embeddings are misaligned using index intersection (lines 256-269)
- Files: `scripts/benchmark_integration.py:256-269`
- Trigger: If `run_method_a/b/c` reorder cells during processing (e.g., via `anndata.concat`), the common index alignment is incorrect
- Workaround: Current code uses `get_indexer()` to align but relies on index being consistent across all method results. If raw data is modified differently per method, alignment fails silently

**Macrophage polarization loss in myeloid_analysis.py:**
- Symptoms: Line 102 computes polarization proportions but doesn't store result, so it's computed but discarded
- Files: `scripts/myeloid_analysis.py:99-102`
- Trigger: Always (computation is never used)
- Impact: Dead code; removed in refactored version but indicates incomplete implementation

## Test Coverage Gaps

**Minimal test coverage:**
- What's not tested: None of the actual pipeline figure generation is tested; only imports are checked
- Files:
  - `tests/test_imports.py` - Only imports scripts (no assertions)
  - `tests/test_units.py` - Tests a few individual utility functions
  - `tests/conftest.py` - Mocks all external dependencies so scripts can be imported
- Risk: Figure generation logic untested; assumes output will be correct; errors discovered only in visual inspection
- Priority: **High** — figures are the deliverable

**Mock-heavy tests obscure real failures:**
- What's not tested: Integration with actual data; config parsing with real YAML; actual AnnData I/O
- Files: `tests/conftest.py:206-305` - All tests use mocked imc_analysis, scanpy, seaborn
- Risk: Script may import successfully but fail on real data due to missing columns, incompatible versions, or YAML parsing errors
- Priority: **High** — should have smoke tests with fixture data

## Performance Bottlenecks

**Large exploratory scripts never executed:**
- Problem: `scripts/exploratory/cell_proportion_across_condition.py` (1046 lines) and `rapids_scanpy_funcs.py` (924 lines) require RAPIDS/cuDF/CUDA
- Files: `scripts/exploratory/rapids_scanpy_funcs.py`, `scripts/exploratory/cell_proportion_across_condition.py`
- Cause: GPU-specific code path requires NVIDIA CUDA libraries; not portable to CPU-only environments
- Improvement path: Separate RAPIDS code into optional submodule with graceful fallback to scanpy. Document GPU requirements.

**Inefficient string operations in plots:**
- Problem: Each figure generation script independently loads full AnnData, subsets by celltype, then regenerates embeddings (PCA/UMAP) for each figure
- Files: `scripts/myeloid_analysis.py:49-52`, `scripts/ue_analysis.py:16-57`
- Cause: No shared preprocessing cache; O(n) AnnData I/O repeated across scripts
- Improvement path: Pre-compute and cache embeddings in the input AnnData. Or run all figures in a single script with shared data.

## Fragile Areas

**Script module-level side effects:**
- Files: All scripts in `scripts/*.py` (e.g., `epithelial_characterization.py:15-85`)
- Why fragile: Code runs at import time in `if __name__ == '__main__'` block, but relies on global state (config, file paths, mocked data). If imports change, side effects change. Makes scripts hard to test and reuse.
- Safe modification: Extract all figure-generation logic into functions; keep `if __name__ == '__main__'` as thin wrapper calling functions. This allows `test_imports.py` to import without executing.

**Config structure assumptions throughout codebase:**
- Files: Multiple scripts assume nested dict structure: `metadata['PANEL_G']['AnnData']['myeloids']`
- Why fragile: If config schema changes, all scripts must be updated. No schema validation.
- Safe modification: Add a Pydantic model or dataclass for config schema in `utils.py`. Validate at startup.

**Celltype name coupling:**
- Files: `scripts/epithelial_characterization.py:26-32`, `scripts/myeloid_analysis.py:30-32`, `scripts/patient_group.py:65-72`
- Why fragile: Cell type names hardcoded in multiple places (e.g., `'Tumor-like (RAGE+)'`, `'Mac. (CD163+)'`). If cell typing logic changes, must update all scripts
- Safe modification: Store celltype mappings in config YAML or a constants module. Import once in utils.

**File path assumptions:**
- Files: All scripts assume relative paths work (e.g., `'figures/figure1/'`, `metadata['phenotyped_umap_name']`)
- Why fragile: If scripts are run from different directories, paths break. No normalization of paths.
- Safe modification: Use `Path(__file__).resolve().parent` to find script directory, then resolve all paths relative to project root stored in config.

## Missing Critical Features

**No logging infrastructure:**
- Problem: Scripts use `print()` and `sc.logging.info()` inconsistently; no structured logging
- Impact: Cannot audit which steps ran, how long they took, or capture warnings/errors systematically
- Blocks: Future CI/CD, batch job monitoring, reproducibility audits

**No progress tracking or checkpointing:**
- Problem: All scripts recompute from scratch on each run; no intermediate caching
- Impact: 50k cell subsampling in benchmark scripts takes time; rerunning any figure regenerates all inputs
- Blocks: Iterative development; quick re-renders of figures

**No input/output contracts:**
- Problem: No clear documentation of required AnnData columns, expected config keys, or output file structure
- Impact: Difficult to validate inputs; scripts fail with cryptic errors; hard to parallelize or reuse components
- Blocks: Scaling pipeline; automation; collaboration

## Security Considerations

**Config downloads over HTTP (not HTTPS):**
- Risk: Config file downloaded from `CONFIG_URL = 'https://wcm.box.com/...'` but URL uses Box shared link (not direct HTTPS with cert pinning)
- Files: `scripts/utils.py:22`
- Current mitigation: Downloaded over HTTPS; Box is trusted provider
- Recommendations: Verify Box cert chain; consider storing config in git or a secure artifact store; do not embed credentials in config URL

**No input sanitization:**
- Risk: Config keys and column names used directly in file paths (e.g., `f'figures/{figure_number}/'`); could allow directory traversal if config is compromised
- Files: `scripts/celltype_differential_abundance.py:61`
- Current mitigation: Config is trusted internal file; column names from AnnData (internal data)
- Recommendations: Validate file paths before writing; use `pathlib.Path.resolve()` to prevent `..` traversal

## Scaling Limits

**Benchmark scripts designed for ~50k cells:**
- Current capacity: Designed to subsample to 50k cells; full panel may be 100k+
- Limit: Integration benchmark may become slow for larger panels or if multiple panels run in parallel
- Scaling path: Profile benchmark scripts; consider streaming benchmarks or distributed comparison metrics

**All figures regenerated on every run:**
- Current capacity: ~10-15 minutes for all scripts to run sequentially
- Limit: If more panels added or scripts parallelized, need caching strategy
- Scaling path: Implement Snakemake rules or DVC for caching; memoize expensive steps (PCA, UMAP)

## Dependencies at Risk

**Exploratory RAPIDS code may become obsolete:**
- Risk: `scripts/exploratory/rapids_scanpy_funcs.py` requires cuML, cuDF, cuPy (NVIDIA RAPIDS ecosystem); these are GPU-specific and may not be maintained
- Impact: Cannot run on CPU-only clusters; breaks if CUDA is unavailable
- Migration plan: Rewrite GPU-heavy operations using scanpy + sklearn + numpy; extract into optional `scripts/gpu/` module with clear documentation of when to use

**imc_analysis dependency is not vendored:**
- Risk: Project depends on external `imc_analysis` library (likely in parent project); if that changes API, scripts break
- Impact: Cannot reproduce results if imc_analysis version changes
- Migration plan: Pin imc_analysis version in requirements. Consider vendoring key functions into utils.py if API is unstable.

**Seaborn version dependency:**
- Risk: Scripts use seaborn KDE plots and `despine()` without version pinning; may break if seaborn changes plot defaults
- Impact: Figures may look different between runs with different seaborn versions
- Migration plan: Pin matplotlib and seaborn versions in requirements.txt

## Dependencies Versions

**No requirements.txt or environment specification:**
- Problem: No `requirements.txt`, `environment.yml`, or `pyproject.toml` with pinned versions
- Impact: Reproducing results with different versions of scanpy, anndata, matplotlib will produce different figures
- Fix approach: Create `requirements.txt` with `pip freeze` or use Poetry/conda for environment management

**Inconsistent import patterns:**
- Problem: Some scripts import from `imc_analysis as imc`, others use `import imc_analysis as imc`; some use relative imports (from utils)
- Impact: Harder to track dependencies; makes it harder to refactor
- Fix approach: Standardize to `import imc_analysis as imc` and `from scripts.utils` in exploratory notebooks

---

*Concerns audit: 2026-03-19*
