# Coding Conventions

**Analysis Date:** 2026-03-19

## Naming Patterns

**Files:**
- Snake case: `epithelial_characterization.py`, `roi_pca_plot.py`, `celltype_heatmap_info.py`
- Descriptive names indicating figure/analysis type or utility purpose
- Shared utilities module: `utils.py`

**Functions:**
- Snake case: `load_config()`, `load_panels()`, `save_figure()`, `ensure_dir()`, `cond_prob()`
- Verb-noun pattern for action functions: `load_*`, `save_*`, `ensure_*`

**Variables:**
- Snake case for local variables: `adata_dict`, `h_myeloid`, `cytokine_expression`, `pathology_colors`
- Concise abbreviations for AnnData objects: `adata`, `pg_epi`, `h_myeloid`, `g_myeloid`
- Dictionary keys match data structure names from config: `'PANEL_G'`, `'PANEL_H'`, `'AnnData'`

**Types:**
- Modern type hints used throughout: `list[str]`, `dict`, `str | None`
- From `__future__ import annotations` enables PEP 563 postponed evaluation

## Code Style

**Formatting:**
- Python 3.10+ (requires-python = ">=3.10")
- No explicit formatter configured (no .black, .isort, .flake8 config files)
- Lines are generally kept under 100 characters
- 4-space indentation (standard Python)

**Linting:**
- No active linting config detected
- No explicit pre-commit hooks or CI linting rules in place

**Imports:**
- `from __future__ import annotations` at top of most files for forward compatibility
- Standard library imports first (`os`, `sys`, `logging`, `argparse`)
- Third-party imports second (`pandas`, `numpy`, `scanpy`, `matplotlib`, `imc_analysis`)
- Local imports last (`from utils import ...`)
- Imports typically grouped by category but not explicitly organized with blank lines

**Path Aliases:**
- No path aliases configured (PYTHONPATH includes `scripts/` via pytest.ini)
- Relative imports from same directory used when needed: `from utils import load_config`

## Error Handling

**Patterns:**
- Explicit validation with `ValueError` for configuration errors: `raise ValueError(f"Panel '{p}' not found in config...")`
- Validation occurs early in loading functions (e.g., `load_panels()` validates all panels and keys before loading)
- File existence checks with `os.path.exists()` before attempting file operations
- Informative error messages include context and available options: `f"Available keys: {sorted(metadata.keys())}"`
- No try-except blocks in main pipeline scripts (assume upstream data is valid)

**Logging:**
- Uses scanpy logging: `sc.logging.info(f"Reading {path}...")`
- Informative messages for file loading operations
- No explicit logging configuration

## Comments

**When to Comment:**
- Module docstrings describe script purpose and major sections: `"""Figure 3: epithelial phenotype proportions and EMT (PanCK/Vimentin) density plots."""`
- Inline comments explain non-obvious grouping/filtering logic
- Section headers use comment blocks (dashed lines) to separate major phases

**Example:**
```python
# ------------------------------------------------------------------
# PANEL_H — epithelial subtype cell proportion
# ------------------------------------------------------------------
```

**JSDoc/TSDoc:**
- NumPy-style docstrings used for functions with parameters and returns
- Format: description, Parameters section, Returns section, Raises section (if applicable)
- Example from `utils.py`:
```python
def load_panels(
    metadata: dict,
    panels: list[str] = PANELS,
    *,
    anndata_key: str = 'phenotyped_umap_name',
    backup_key: str = 'backup_url',
) -> dict:
    """Load per-panel AnnData objects from the config.

    Reads ``metadata[panel]['AnnData'][anndata_key]`` for each panel and
    falls back to ``metadata[panel]['AnnData'][backup_key]`` if the local
    file is absent.

    Parameters
    ----------
    metadata:
        Parsed config dict (from :func:`load_config`).
    panels:
        Which panels to load, e.g. ``['PANEL_G', 'PANEL_H']``.
    anndata_key:
        Key inside ``metadata[panel]['AnnData']`` for the local file path.
    backup_key:
        Key inside ``metadata[panel]['AnnData']`` for the remote backup URL.

    Returns
    -------
    dict
        Mapping ``panel -> AnnData``.

    Raises
    ------
    ValueError
        If a panel name or required key is missing from *metadata*.
    """
```

## Function Design

**Size:**
- Pipeline scripts use `if __name__ == '__main__':` block for execution
- Main blocks typically span 50-100 lines for complex analyses
- Shared utilities extracted to `utils.py` for reuse across scripts
- Single-purpose functions preferred: `load_config()`, `save_figure()`, `ensure_dir()`

**Parameters:**
- Keyword-only parameters after `*` to enforce explicit calling: `def save_figure(path: str, *, fig: plt.Figure | None = None, ...)`
- Default values provided for optional parameters
- Type hints on all parameters and return values

**Return Values:**
- Functions return data directly (e.g., `dict` from `load_config()`)
- Side effects (file I/O) clearly documented in docstrings
- No implicit None returns; functions either return data or explicitly state side-effect-only behavior

## Module Design

**Exports:**
- `utils.py` exports constants: `CONFIG_PATH`, `CONFIG_URL`, `PANELS`, `CYTOKINE`
- All functions in `utils.py` are public (no underscore prefix)
- Pipeline scripts define one main execution block

**Constants:**
- Module-level constants defined near imports: `CONFIG_PATH = 'metadata/ggo_config.yml'`
- Listed alphabetically or by logical grouping
- All caps for global constants: `PANELS`, `CYTOKINE`

**Data Structures:**
- Configuration loaded as plain `dict` from YAML
- AnnData objects for biological data
- Pandas DataFrames for tabular results
- Matplotlib figures for plotting

---

*Convention analysis: 2026-03-19*
