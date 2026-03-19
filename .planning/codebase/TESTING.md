# Testing Patterns

**Analysis Date:** 2026-03-19

## Test Framework

**Runner:**
- pytest (configured in `pyproject.toml` with `[tool.pytest.ini_options]`)
- Config: `pyproject.toml` (lines 15-16)

**Assertion Library:**
- pytest's built-in assertions and pytest.raises context manager

**Run Commands:**
```bash
pytest                      # Run all tests
pytest -v                   # Verbose output
pytest tests/test_imports.py -v  # Run specific test file
pytest tests/test_units.py::TestEnsureDir -v  # Run specific test class
```

## Test File Organization

**Location:**
- Co-located in `tests/` directory (separate from source)
- Directory structure: `tests/conftest.py`, `tests/__init__.py`, `tests/test_imports.py`, `tests/test_units.py`

**Naming:**
- Test files prefixed with `test_`: `test_imports.py`, `test_units.py`
- Test functions prefixed with `test_`: `test_import_utils()`, `test_creates_directory()`
- Test classes use `Test` prefix: `TestEnsureDir`, `TestLoadPanelsValidation`, `TestCondProb`

**Structure:**
```
tests/
├── __init__.py              # Package marker
├── conftest.py              # Shared pytest fixtures and mocks
├── test_imports.py          # Smoke tests for script imports
└── test_units.py            # Unit tests for utility functions
```

## Test Structure

**Suite Organization:**
Test modules organize tests by function or component. Example from `test_units.py`:

```python
class TestEnsureDir:
    def test_creates_directory(self, tmp_path):
        from utils import ensure_dir
        target = str(tmp_path / "a" / "b" / "c")
        ensure_dir(target)
        assert (tmp_path / "a" / "b" / "c").is_dir()

    def test_returns_path(self, tmp_path):
        from utils import ensure_dir
        target = str(tmp_path / "x")
        result = ensure_dir(target)
        assert result == target

    def test_idempotent(self, tmp_path):
        from utils import ensure_dir
        target = str(tmp_path / "d")
        ensure_dir(target)
        ensure_dir(target)  # second call must not raise


class TestLoadPanelsValidation:
    def test_missing_panel_raises(self):
        from utils import load_panels
        with pytest.raises(ValueError, match="not found in config"):
            load_panels({"PANEL_G": {"AnnData": {}}}, panels=["PANEL_X"])
```

**Patterns:**
- One class per function or cohesive feature being tested
- Each test method tests a single aspect (arrange-act-assert)
- Imports deferred to test methods to work with conftest fixtures
- Direct imports from `utils` and pipeline modules by name (no relative paths) due to PYTHONPATH setup in conftest

## Mocking

**Framework:**
- `unittest.mock` (Python standard library)
- `MagicMock` for creating mock objects with realistic behavior

**Patterns:**
From `conftest.py`, mocks are created as factory functions:

```python
def _make_imc_mock() -> MagicMock:
    imc = MagicMock()
    imc.utils.parse_yaml.return_value = FAKE_CONFIG
    imc.utils.download.return_value = None
    # celltype_density returns an object whose .to_df() gives a usable DataFrame
    _ldens = pd.DataFrame({
        "CD8 T": np.random.rand(5),
        "CD4 T": np.random.rand(5),
        "T reg": np.random.rand(5),
    })
    _area_obs = pd.DataFrame(
        {"pathology": ["AAH"] * 4, "radio": ["pGGO"] * 4},
        index=[f"roi{i + 1}" for i in range(4)],
    )
    _area_var = pd.DataFrame(index=["A", "B", "C"])
    imc.tl.celltype_density.return_value.to_df.return_value = _ldens
    imc.tl.celltype_density.return_value.obs = _area_obs
    imc.tl.celltype_density.return_value.var = _area_var
    return imc
```

**What to Mock:**
- External dependencies: `imc_analysis`, `scanpy`, `seaborn`
- File I/O operations that would require real data files
- Heavy computations and UMAP/PCA calculations
- Network operations (YAML downloads)

**What NOT to Mock:**
- Core utility functions in `utils.py` (test these directly)
- Pandas operations (test with real DataFrames)
- AnnData object construction and basic operations
- Pytest fixtures and built-in context managers

**Fixture Management:**
- Session-scoped setup for persistent mocks (reused across all tests)
- Per-test cleanup via autouse fixtures
- Per-test side effects for specific test variations (e.g., `roi_pca_plot_group` requires different mock behavior)

## Fixtures and Factories

**Test Data:**
The conftest provides factory functions that construct realistic fake AnnData objects matching the pipeline requirements:

```python
def _build_obs(include_group: bool = True) -> pd.DataFrame:
    """Build the obs DataFrame used in fake AnnData objects."""
    ct_all = list(_CELLTYPES[:_N])
    cb_all = list(_CELLTYPE_BROAD[:_N])
    # Override the last 4 rows' celltype_broad to support epithelial_characterization L69
    cb_all[20] = "Epithelial-like"
    cb_all[21] = "Epithelial-like (Ki67+)"
    cb_all[22] = "Mesenchymal-like"
    cb_all[23] = "Epithelial-like"

    d = {
        "radio": pd.Categorical(...),
        "pathology": pd.Categorical(...),
        "celltype": pd.Categorical(...),
        # ... more obs columns
    }
    if include_group:
        d["Group"] = pd.Categorical(...)
    return pd.DataFrame(d)


def _make_fake_adata(include_group: bool = True) -> anndata.AnnData:
    obs = _build_obs(include_group=include_group)
    X = np.random.rand(_N, len(_ALL_MARKERS)).astype(np.float32)
    var = pd.DataFrame(index=list(_ALL_MARKERS))
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    # ... add obsm, uns
    return adata
```

**Location:**
- Factories defined in `conftest.py` (lines 142-200)
- Shared constants at module level: `_CELLTYPES`, `_CELLTYPE_BROAD`, `_ALL_MARKERS`, `_GGO_IDS`

## Test Types

### Smoke Tests (test_imports.py)

**Scope and Approach:**
- Verify all pipeline scripts can be imported without error
- One test per script: `test_import_utils()`, `test_import_epithelial_characterization()`
- Tests that scripts work with mocked dependencies (no real data required)
- Catches module-level import errors and syntax issues early

**Example:**
```python
def test_import_epithelial_characterization():
    import epithelial_characterization  # noqa: F401
```

### Unit Tests (test_units.py)

**Scope and Approach:**
- Test pure functions from `utils.py` and `patient_group.py` in isolation
- Use pytest fixtures (`tmp_path` for filesystem testing)
- Test both success cases and error conditions
- Organize by class (one class per function being tested)

**Example - Error Validation:**
```python
class TestLoadPanelsValidation:
    def test_missing_panel_raises(self):
        from utils import load_panels
        with pytest.raises(ValueError, match="not found in config"):
            load_panels({"PANEL_G": {"AnnData": {}}}, panels=["PANEL_X"])

    def test_missing_anndata_section_raises(self):
        from utils import load_panels
        with pytest.raises(ValueError, match="AnnData.*missing"):
            load_panels({"PANEL_G": {}}, panels=["PANEL_G"])
```

**Example - Data Transformation:**
```python
class TestCondProb:
    def test_output_columns(self, sample_obs):
        from patient_group import cond_prob
        result = cond_prob(sample_obs, y="pathology", x="radio")
        assert "radio" in result.columns
        pathologies = set(sample_obs["pathology"].unique())
        assert pathologies.issubset(set(result.columns))

    def test_rows_sum_to_100(self, sample_obs):
        from patient_group import cond_prob
        result = cond_prob(sample_obs, y="pathology", x="radio")
        numeric = result.drop(columns="radio").fillna(0)
        row_sums = numeric.sum(axis=1)
        assert (row_sums.round(4) == 100.0).all()
```

### Integration Tests

**Status:** Not detected in codebase. Smoke tests serve as minimal integration verification (scripts run without error).

## Coverage

**Requirements:** Not enforced

**View Coverage:**
```bash
pytest --cov=scripts --cov-report=term-only tests/
```

(Coverage tracking not configured but can be added with pytest-cov plugin)

## Common Patterns

### Async Testing

Not applicable (synchronous Python scripts, no async/await)

### Error Testing

Tests use `pytest.raises()` context manager with optional `match` parameter for regex matching on error messages:

```python
with pytest.raises(ValueError, match="not found in config"):
    load_panels({"PANEL_G": {"AnnData": {}}}, panels=["PANEL_X"])
```

### Fixture Patterns

**Autouse Fixtures:**
Two session-wide autouse fixtures in `conftest.py` (lines 247-306):

```python
@pytest.fixture(autouse=True)
def _install_global_mocks():
    """Install imc_analysis and seaborn mocks for every test."""
    imc_mock = _make_imc_mock()
    sns_mock = MagicMock()

    _originals = {}
    for name, obj in [("imc_analysis", imc_mock), ("seaborn", sns_mock)]:
        _originals[name] = sys.modules.get(name)
        sys.modules[name] = obj

    yield imc_mock, sns_mock

    # Restore originals
    for name, orig in _originals.items():
        if orig is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = orig
    # Clear pipeline script modules so next test gets a fresh import
    for script in PIPELINE_SCRIPTS:
        sys.modules.pop(script, None)
```

**Per-Test Fixtures:**
`_install_sc_mock` handles test-specific scanpy mock behavior (e.g., different return values for multiple `sc.read()` calls in `roi_pca_plot_group`)

### Cleanup Strategy

- After each test, pipeline script modules are cleared from `sys.modules` (line 267-268)
- Allows fresh import of scripts with fresh mocks for next test
- Prevents stale module state from interfering with subsequent tests

---

*Testing analysis: 2026-03-19*
