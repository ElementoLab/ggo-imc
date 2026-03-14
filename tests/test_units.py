"""Unit tests for pure functions in utils.py and patient_group.py."""
from __future__ import annotations

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# utils.py tests
# ---------------------------------------------------------------------------

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

    def test_missing_anndata_section_raises(self):
        from utils import load_panels
        with pytest.raises(ValueError, match="AnnData.*missing"):
            load_panels({"PANEL_G": {}}, panels=["PANEL_G"])

    def test_missing_anndata_key_raises(self):
        from utils import load_panels
        metadata = {"PANEL_G": {"AnnData": {"backup_url": "http://example.com"}}}
        with pytest.raises(ValueError, match="phenotyped_umap_name"):
            load_panels(metadata, panels=["PANEL_G"])


# ---------------------------------------------------------------------------
# patient_group.py tests
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_obs():
    return pd.DataFrame({
        "radio": ["pGGO", "pGGO", "mGGO", "mGGO", "Solid", "Solid"],
        "pathology": ["AAH", "AIS", "MIA", "LUAD", "MIA", "LUAD"],
    })


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
        assert (row_sums.round(4) == 100.0).all(), f"Row sums: {row_sums.tolist()}"

    def test_one_row_per_radio_group(self, sample_obs):
        from patient_group import cond_prob
        result = cond_prob(sample_obs, y="pathology", x="radio")
        assert len(result) == sample_obs["radio"].nunique()
