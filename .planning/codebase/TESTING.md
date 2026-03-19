# Testing Patterns

**Analysis Date:** 2026-03-19

## Test Framework

**Runner:**
- pytest 7.x+ (indicated by `.pytest_cache` directory and `conftest.py` structure)
- Config: `pyproject.toml` with `[tool.pytest.ini_options]`

**Assertion Library:**
- Standard `assert` statements
- pytest built-in assertions with match patterns: `pytest.raises(ValueError, match="not found in config")`

**Run Commands:**
```bash
pytest                     # Run all tests
pytest tests/test_imports.py          # Run smoke tests only
pytest tests/test_units.py            # Run unit tests only
pytest -v                  # Verbose output
pytest --tb=short          # Short traceback format
```

Note: Tests are run from the project root; pytest.ini_options configures `pythonpath = ["scripts"]` to make pipeline scripts importable.

## Test File Organization

**Location:**
- `tests/` directory in project root

**Naming:**
- `test_imports.py` - smoke tests for all 10 pipeline scripts
- `test_units.py` - unit tests for pure functions
- `conftest.py` - pytest configuration and fixtures

**Structure:**
```
tests/
├── __init__.py
├── conftest.py          # Global fixtures, mocks setup, fake data builders
├── test_imports.py      # Smoke tests (one test per pipeline script)
└── test_units.py        # Unit tests grouped in test classes
```

## Test Structure

**Suite Organization:**
- Smoke tests: 10 simple tests, one per pipeline script to verify imports
- Unit tests: Grouped by class using pytest test classes

**Smoke Test Pattern:**
```python
def test_import_epithelial_characterization():
    import epithelial_characterization  # noqa: F401
```

**Unit Test Structure:**
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
```

**Patterns:**
- Test classes group related tests: `TestEnsureDir`, `TestLoadPanelsValidation`, `TestCondProb`
- Each test method tests one behavior
- Fixture functions provide test data via pytest's `@pytest.fixture` decorator
- Test-specific imports inside test methods allow fresh mocks between tests

## Mocking

**Framework:**
- Python's built-in `unittest.mock.MagicMock`
- Patching via `sys.modules` replacement (not `@patch` decorator)

**Mock Strategy:**
The test setup (`conftest.py`) creates lightweight mocks that allow pipeline scripts to import without loading real data. Mocks persist across tests but are reset between tests via `sys.modules` cleanup.

**Key Mocks:**
1. `imc_analysis` (imc_mock) - Returns FAKE_CONFIG from `parse_yaml()`
2. `scanpy` (sc_mock) - Returns fake AnnData objects with all required obs columns/var names
3. `seaborn` (sns_mock) - Pure mock to prevent rendering attempts with fake data

**Patterns:**
```python
def _make_imc_mock() -> MagicMock:
    imc = MagicMock()
    imc.utils.parse_yaml.return_value = FAKE_CONFIG
    imc.utils.download.return_value = None
    imc.tl.celltype_density.return_value.to_df.return_value = _ldens
    imc.tl.celltype_density.return_value.obs = _area_obs
    imc.tl.celltype_density.return_value.var = _area_var
    return imc

def _make_sc_mock(read_side_effect=None) -> MagicMock:
    import scanpy as _real_sc
    sc = MagicMock()
    sc.logging = _real_sc.logging  # Use real scanpy logging for console output
    sc.settings = _real_sc.settings
    sc.pp.subsample = MagicMock(return_value=_make_fake_adata())
    sc.tl.score_genes = MagicMock(return_value=None)
    sc.read = MagicMock(return_value=_make_fake_adata())
    return sc
```

**What to Mock:**
- External libraries that scripts import at module level
- Data-loading functions that would hit disk/network
- Plotting functions to avoid rendering

**What NOT to Mock:**
- Core Python functionality (os, sys, pandas, numpy)
- Functions under test (utils.py functions are real)
- Real scanpy logging (use real `sc.logging` to test logging output)

## Fixtures and Factories

**Test Data:**
Test data is built dynamically in conftest.py using factory functions:

```python
def _build_obs(include_group: bool = True) -> pd.DataFrame:
    """Build the obs DataFrame used in fake AnnData objects."""
    d = {
        "radio": pd.Categorical(
            [_RADIO_CATS[i % 4] for i in range(_N)],
            categories=_RADIO_CATS,
        ),
        "pathology": pd.Categorical(
            [_PATHO_CATS[i % len(_PATHO_CATS)] for i in range(_N)],
            categories=_PATHO_CATS,
        ),
        "celltype": pd.Categorical(ct_all, categories=sorted(set(ct_all))),
        # ... more columns required by scripts
    }
    if include_group:
        d["Group"] = pd.Categorical(
            [f"Group{i % 4 + 1}" for i in range(_N)],
            categories=["Group1", "Group2", "Group3", "Group4"],
        )
    return pd.DataFrame(d)

def _make_fake_adata(include_group: bool = True) -> anndata.AnnData:
    obs = _build_obs(include_group=include_group)
    X = np.random.rand(_N, len(_ALL_MARKERS)).astype(np.float32)
    var = pd.DataFrame(index=list(_ALL_MARKERS))
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_pca"] = np.random.rand(_N, 10)
    # ...
    return adata
```

**Key Design Detail:** Test AnnData includes every marker and obs column that ANY pipeline script requires. The `_ALL_MARKERS` set is comprehensive, and obs columns are built to satisfy all script constraints (e.g., specific celltype categories for epithelial_characterization).

**Location:**
- Factory functions in `tests/conftest.py` (lines 142-199)
- `_build_obs()`, `_make_fake_adata()`, `_make_imc_mock()`, `_make_sc_mock()`

**Fixture Functions:**
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

@pytest.fixture(autouse=True)
def _install_sc_mock(request):
    """Install a scanpy mock for every test."""
    test_name = request.node.name

    if test_name == "test_import_roi_pca_plot_group":
        # Special handling: first sc.read returns no-Group, second returns with-Group
        _call_count = [0]
        def _roi_pca_read(*args, **kwargs):
            idx = _call_count[0]
            _call_count[0] += 1
            return _no_group if idx == 0 else _with_group
        sc_mock = _make_sc_mock(read_side_effect=_roi_pca_read)
    else:
        sc_mock = _make_sc_mock()

    orig = sys.modules.get("scanpy")
    sys.modules["scanpy"] = sc_mock

    yield sc_mock

    if orig is None:
        sys.modules.pop("scanpy", None)
    else:
        sys.modules["scanpy"] = orig
```

**Custom Fixtures:**
```python
@pytest.fixture
def sample_obs():
    return pd.DataFrame({
        "radio": ["pGGO", "pGGO", "mGGO", "mGGO", "Solid", "Solid"],
        "pathology": ["AAH", "AIS", "MIA", "LUAD", "MIA", "LUAD"],
    })
```

## Coverage

**Requirements:** No coverage threshold enforced in `pyproject.toml`

**View Coverage:**
```bash
pytest --cov=scripts --cov-report=html
pytest --cov=scripts --cov-report=term
```

Note: Coverage measurement is not automated; pytest-cov can be installed separately if needed.

## Test Types

**Smoke Tests:**
- Type: Import verification
- Scope: Each of 10 pipeline scripts (utils, celltype_heatmap_info, roi_pca_plot, etc.)
- Approach: Import script at module level, verify no errors with mocked dependencies
- Location: `tests/test_imports.py` (lines 4-42)
- Purpose: Catch module-level syntax errors and import failures before execution

**Unit Tests:**
- Type: Pure function testing
- Scope: `ensure_dir()`, `load_panels()`, `cond_prob()`
- Approach: Test with real pandas/numpy data; no mocking required
- Location: `tests/test_units.py` (lines 12-81)
- Purpose: Verify data validation and transformation logic

**Example Unit Test:**
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

**Integration Tests:** Not present. Scripts perform integration at runtime; no integration test suite.

**E2E Tests:** Not applicable (scripts are exploratory/analytical, not user-facing APIs).

## Common Patterns

**Async Testing:**
Not used (Python 3.10 synchronous code only).

**Error Testing:**
```python
def test_missing_panel_raises(self):
    from utils import load_panels
    with pytest.raises(ValueError, match="not found in config"):
        load_panels({"PANEL_G": {"AnnData": {}}}, panels=["PANEL_X"])
```

**Fixture Setup/Teardown:**
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

    # Cleanup/restoration logic
    for name, orig in _originals.items():
        if orig is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = orig
    for script in PIPELINE_SCRIPTS:
        sys.modules.pop(script, None)
```

## Special Test Design: roi_pca_plot_group

The test for `roi_pca_plot_group.py` has special mocking because the script calls `sc.read()` twice with different requirements:

```python
if test_name == "test_import_roi_pca_plot_group":
    _call_count = [0]
    _no_group = _make_fake_adata(include_group=False)
    _with_group = _make_fake_adata(include_group=True)

    def _roi_pca_read(*args, **kwargs):
        idx = _call_count[0]
        _call_count[0] += 1
        return _no_group if idx == 0 else _with_group

    sc_mock = _make_sc_mock(read_side_effect=_roi_pca_read)
```

This ensures the merge operation inside `roi_pca_plot_group.py` (line 30-31) succeeds by returning appropriate AnnData structures at each call.

## Constraints Encoded in Test Data

The fake AnnData is carefully constructed to satisfy all constraints across 9 pipeline scripts:
- `_N = 24` (divisible by 4 pathologies and 3 radio categories)
- `_CELLTYPES` includes all cell types referenced by epithelial_characterization, myeloid_analysis, etc.
- `_GGO_IDS` are unique to prevent merge errors in roi_pca_plot_group
- Categorical columns have matching categories across all scripts (e.g., radio: ["pGGO", "mGGO", "Solid", "pGGO (partial)"])

---

*Testing analysis: 2026-03-19*
