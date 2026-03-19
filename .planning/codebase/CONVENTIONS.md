# Coding Conventions

**Analysis Date:** 2026-03-19

## Naming Patterns

**Files:**
- Snake case: `epithelial_characterization.py`, `t_cell_analysis.py`, `celltype_heatmap_info.py`
- Figure scripts named by content: `figure1/roi_pca_plot.py`, `figure2/t_cell_analysis.py`, `figure5/patient_group.py`
- Utility modules: `utils.py` (shared functions), `conftest.py` (test configuration)
- Exploratory scripts in `scripts/exploratory/` subdirectory

**Functions:**
- Snake case: `load_config()`, `load_panels()`, `load_single_panel()`, `save_figure()`, `ensure_dir()`, `cond_prob()`
- Private helpers use underscore prefix: `_build_obs()`, `_make_fake_adata()`, `_make_imc_mock()`, `_make_sc_mock()`
- Test fixtures prefixed with underscore: `_install_global_mocks()`, `_install_sc_mock()`

**Variables:**
- Snake case for local variables: `h_myeloid`, `g_lymphocytes`, `celltype_map`, `phenotype`
- Constants in UPPER_CASE: `CONFIG_PATH`, `CONFIG_URL`, `PANELS`, `CYTOKINE`
- Dictionary keys use appropriate case for context: `'pathology'`, `'radio'`, `'celltype'`, `'celltype_broad'`

**Types:**
- Type hints in docstrings and function signatures: `dict`, `list[str]`, `str | None`, `AnnData`
- Modern PEP 604 union syntax in type hints: `str | None` instead of `Optional[str]`

## Code Style

**Formatting:**
- No explicit formatter configured (black/autopep8)
- Consistent 4-space indentation
- Docstrings use triple-quoted format
- Long function signatures break across lines with keyword-only arguments after `*`

**Linting:**
- No explicit linter configured (flake8/ruff/pylint)
- Follows PEP 8 conventions informally
- Imports include `# noqa: F401` for intentionally unused imports in tests

## Import Organization

**Order:**
1. Standard library: `sys`, `os`, `logging`, `from pathlib import Path`
2. Third-party core: `numpy as np`, `pandas as pd`
3. Bioinformatics libraries: `scanpy as sc`, `anndata`
4. Visualization: `matplotlib.pyplot as plt`, `seaborn as sns`
5. Project-specific: `imc_analysis as imc`, `from utils import ...`

**Path Aliases:**
- No path aliases configured in pyproject.toml
- Relative imports from `scripts/` work because `sys.path` is prepended in test setup
- Bare `from utils import load_config` works in pipeline scripts (set in `pytool.pytest.ini_options.pythonpath`)

**Example from `myeloid_analysis.py`:**
```python
import os

import pandas as pd
import numpy as np

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import imc_analysis as imc

from utils import load_config, load_panels, save_figure, CYTOKINE as cytokine
```

## Error Handling

**Validation patterns:**
- Explicit validation in utility functions with descriptive `ValueError` messages
- Example from `utils.py`:
```python
def load_panels(metadata: dict, panels: list[str] = PANELS, ...):
    for p in panels:
        if p not in metadata:
            raise ValueError(
                f"Panel '{p}' not found in config. "
                f"Available keys: {sorted(metadata.keys())}"
            )
        anndata_section = metadata[p].get('AnnData')
        if anndata_section is None:
            raise ValueError(
                f"Config section metadata['{p}']['AnnData'] is missing."
            )
```

**Error handling in scripts:**
- Conditional checks using `os.path.exists()` to decide between loading cached and computing
- Example from `myeloid_analysis.py`:
```python
if os.path.exists(metadata['PANEL_H']['AnnData']['myeloids']):
    h_myeloid = sc.read(...)
else:
    sc.pp.scale(h_myeloid)
    # ... processing ...
    h_myeloid.write(metadata['PANEL_H']['AnnData']['myeloids'])
```

**Exception handling:**
- Explicit try-except in benchmarking/intensive scripts
- Example from `benchmark_integration.py`:
```python
try:
    # process
except ZeroDivisionError:
    logger.warning("ZeroDivisionError: skipping ROI %s", roi_name)
except Exception:
    logger.exception("Method A failed")
```

## Logging

**Framework:**
- Primarily `scanpy` logging: `sc.logging.info(f"Reading {path}...")`
- Standard `logging` module in benchmarking scripts: `logger = logging.getLogger(__name__)`

**Patterns:**
- Scanpy logging for file I/O and data loading: `sc.logging.info(f"Reading {path}...")`
- Print statements for optional command-line arguments: `print(sys.argv)`, `print(f'Subsetting samples for cell lineage: "{x}"')`
- Standard logging module with INFO/WARNING/ERROR levels in computation-heavy scripts
- Setup: `logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")`

## Comments

**When to Comment:**
- Section headers using block comments with dashes:
```python
# ------------------------------------------------------------------
# PANEL_H — epithelial subtype cell proportion
# ------------------------------------------------------------------
```
- Module docstrings at top of every script describing purpose
- Complex data transformations and business logic

**JSDoc/TSDoc:**
- Uses NumPy-style docstrings (not Google or Sphinx)
- Example from `utils.py`:
```python
def load_config(
    config_path: str = CONFIG_PATH,
    *,
    download_url: str | None = CONFIG_URL,
) -> dict:
    """Load the project YAML config, optionally downloading it first.

    Parameters
    ----------
    config_path:
        Path to ``ggo_config.yml`` (relative to the working directory).
    download_url:
        If provided and the file does not already exist, download from this
        URL before parsing.  Pass ``None`` to skip the download step.

    Returns
    -------
    dict
        Parsed YAML content.
    """
```

## Function Design

**Size:**
- Utility functions typically 10-20 lines
- Script bodies often 50+ lines (acceptable for exploratory/pipeline code)
- Keep repeated logic in utils; scripts call utilities

**Parameters:**
- Use keyword-only arguments after `*` for optional configuration:
```python
def load_config(
    config_path: str = CONFIG_PATH,
    *,
    download_url: str | None = CONFIG_URL,
) -> dict:
```
- Default parameters preferred over None checks
- Positional arguments for required inputs, keyword-only for options

**Return Values:**
- Explicit return types in signature: `-> dict`, `-> AnnData`, `-> str`, `-> None`
- Single return value per function (no tuple returns in utils)
- Functions return data, not status codes or error codes

## Module Design

**Exports:**
- `utils.py` exports core utilities used by all scripts:
  - Configuration: `load_config()`, `CONFIG_PATH`, `CONFIG_URL`
  - Data loading: `load_panels()`, `load_single_panel()`, `PANELS`
  - Figure saving: `save_figure()`, `ensure_dir()`
  - Constants: `CYTOKINE` list

**Barrel Files:**
- No barrel files (index.py) used
- Direct function imports from modules

## Python Version

**Target:** Python >=3.10 (from `pyproject.toml`)
- Modern type hint syntax allowed: `str | None`, `list[str]`, `dict[str, Any]`
- f-strings with expressions

---

*Convention analysis: 2026-03-19*
