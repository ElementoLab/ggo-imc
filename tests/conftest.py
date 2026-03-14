"""
conftest.py — pytest configuration for ggo-imc smoke tests.

Sets up sys.path and installs lightweight mocks for imc_analysis and scanpy
so that pipeline scripts (which have module-level data-loading code) can be
imported without touching real data files.

Design
------
* ``scripts/`` is prepended to sys.path so bare ``import utils`` etc. work.
* ``imc_analysis`` is replaced with a MagicMock that returns a realistic
  config dict from ``utils.parse_yaml``.
* ``scanpy`` is replaced with a MagicMock; ``sc.read`` returns a pre-built
  AnnData that has every obs column / var name required by any pipeline script.
* ``seaborn`` is replaced with a MagicMock to avoid seaborn trying to render
  figures with MagicMock data.
* Each test clears its target script from ``sys.modules`` so a fresh import
  runs against fresh mocks (important for scripts that call sc.read multiple
  times with different expectations).
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata
import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# sys.path: make scripts/ importable
# ---------------------------------------------------------------------------
SCRIPTS_DIR = str(Path(__file__).parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

# ---------------------------------------------------------------------------
# Config dict returned by imc.utils.parse_yaml in all scripts
# ---------------------------------------------------------------------------
_PATHOLOGY_COLORS = ["#000000", "#111111", "#222222", "#333333", "#444444"]

_PANEL_ANNDATA = {
    "phenotyped_umap_name": "/dev/null",
    "backup_url": "http://example.com",
    "lymphocytes": "/dev/null",
    "lymphocytes_url": "http://example.com",
    "utag_labeled_name": "/dev/null",
    "utag_url": "http://example.com",
    # myeloids paths — set to non-existent so os.path.exists returns False
    "myeloids": "/tmp/_ggo_imc_test_myeloids_nonexistent.h5ad",
    "myeloids_url": "http://example.com",
    "var_celltype_groups": [],
}

FAKE_CONFIG: dict = {
    "pathology_color": _PATHOLOGY_COLORS,
    "Radiology_color": ["#aaaaaa", "#bbbbbb", "#cccccc"],
    "CONDITIONS": ["radio", "pathology"],
    "PCA_ANNDATA_NAME": "/dev/null",
    "PCA_URL": "http://example.com",
    "patient_celltype_broad_clustered": "/dev/null",
    "patient_group_url": "http://example.com",
    # var_celltype_groups is accessed as metadata[p]['var_celltype_groups'] (not under AnnData)
    "PANEL_G": {"AnnData": dict(_PANEL_ANNDATA), "var_celltype_groups": []},
    "PANEL_H": {"AnnData": dict(_PANEL_ANNDATA), "var_celltype_groups": []},
}

# ---------------------------------------------------------------------------
# Pipeline script names (used to clean sys.modules between tests)
# ---------------------------------------------------------------------------
PIPELINE_SCRIPTS = [
    "utils",
    "celltype_heatmap_info",
    "roi_pca_plot",
    "celltype_differential_abundance",
    "t_cell_analysis",
    "myeloid_analysis",
    "epithelial_characterization",
    "ue_analysis",
    "roi_pca_plot_group",
    "patient_group",
]

# ---------------------------------------------------------------------------
# Fake AnnData factory
# ---------------------------------------------------------------------------
_ALL_MARKERS = sorted({
    "CCR4", "CCR7", "GranzymeB", "CD27", "CD28", "ICOS", "Tbet", "41BB",
    "FoxP3", "CTLA4", "PD1", "LAG3", "TIM3",
    "CD56", "NKp44", "CD15", "CD66b", "VISTA", "CD117",
    "CD14", "CD16", "CD68", "CD163", "HLADR", "MMP7", "CD33",
    "IL1alpha", "IL1beta", "IL1R1", "IL12p40", "IL17A", "IL23p19", "IL23R", "IFNg",
    "CD206", "CD11c", "CD11b", "CD57", "PDL1",
    "PanCK", "Vimentin", "SMA", "ColI", "Ecad",
    "RAGE", "SFTPC", "CD31", "FOLR1", "CD3", "CD4", "CD8a", "CD20",
    "Ki67", "CD45RA", "CD45RO", "CD25",
})

_SCORE_COLS = [
    "Epithelial score", "Mesenchymal score", "Fibrosis score",
    "Angiogenesis score", "B&T score", "Pan-Immune score",
]

# celltypes must satisfy:
#   epithelial_characterization L31: .isin(['Tumor-like (RAGE+)','Tumor-like (SFTPC+)','Tumor-like','Epi.-like (RAGE+)','Epi.-like (SFTPC+)'])
#   epithelial_characterization L69: pg_epi filters celltype_broad.isin(['Epithelial-like','Epithelial-like (Ki67+)','Mesenchymal-like'])
#   myeloid_analysis: celltype.isin(['Mac. (CD68+, CD163-)','Mac. (CD163+)','Mast','NK','PMN-MDSC','Neut.'])
#   myeloid_analysis: celltype_broad == 'Macrophage'
_N = 24  # divisible by 4 (pathology) and 3 (radio)
_CELLTYPES = (
    ["Tumor-like (RAGE+)", "Tumor-like (SFTPC+)", "Tumor-like",
     "Epi.-like (RAGE+)", "Epi.-like (SFTPC+)",
     "Mac. (CD68+, CD163-)", "Mac. (CD163+)", "Mast", "NK", "PMN-MDSC", "Neut.",
     "Mac. (CD163+, CD206+)", "CD4 T", "CD8 T", "T reg", "B", "Fib.",
     "Tumor-like (RAGE+)", "Epi.-like (RAGE+)", "Mac. (CD163+)",
     "Mac. (CD68+, CD163-)", "Mast", "NK", "PMN-MDSC"]
)
_CELLTYPE_BROAD = (
    ["Tumor-like", "Tumor-like", "Tumor-like",
     "Epi.-like", "Epi.-like",
     "Macrophage", "Macrophage", "Mast", "NK", "PMN-MDSC", "Neutrophil",
     "Macrophage", "CD4 T", "CD8 T", "T reg", "B", "Fib.",
     "Tumor-like", "Epi.-like", "Macrophage",
     "Macrophage", "Mast", "NK", "PMN-MDSC"]
)
# celltype_broad for epithelial_characterization L69
_CELLTYPE_BROAD_EPI_VALUES = ["Epithelial-like", "Epithelial-like (Ki67+)", "Mesenchymal-like"]

# GGO IDs must be UNIQUE so that density_all[tmp.index] works in roi_pca_plot_group.
# The merge uses left_on='GGO ID', right_index=True; duplicate GGO IDs create a
# many-to-many join which produces a non-unique index, causing AnnData to reject it.
_GGO_IDS = [f"ggo{i + 1}" for i in range(_N)]
# 'N' pathology is required by roi_pca_plot_group.py (line 36: tmp.loc[tmp['pathology']=='N','Group']='N')
# Without 'N' rows, the 'N' group would be empty, causing np.linalg.eig to fail on a NaN cov matrix.
_PATHO_CATS = ["N", "AAH", "AIS", "MIA", "LUAD"]
_RADIO_CATS = ["pGGO", "mGGO", "Solid", "pGGO (partial)"]  # 4 categories required by epithelial_characterization


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
        "radio": pd.Categorical(
            [_RADIO_CATS[i % 4] for i in range(_N)],
            categories=_RADIO_CATS,
        ),
        # Cycle through all 5 pathology categories (i%5) so that
        # the 'N' rows (at i%5==0) do not perfectly align with any single
        # group (i%4==0), avoiding an empty group after roi_pca_plot_group
        # reassigns pathology-N rows to Group='N'.
        "pathology": pd.Categorical(
            [_PATHO_CATS[i % len(_PATHO_CATS)] for i in range(_N)],
            categories=_PATHO_CATS,
        ),
        "celltype": pd.Categorical(ct_all, categories=sorted(set(ct_all))),
        "celltype_broad": pd.Categorical(
            cb_all, categories=sorted(set(cb_all + _CELLTYPE_BROAD_EPI_VALUES))
        ),
        "uE_broad": pd.Categorical(
            [chr(ord("A") + i % 3) for i in range(_N)], categories=["A", "B", "C"]
        ),
        "sample": [f"s{i % 5 + 1}" for i in range(_N)],
        "roi": [f"roi{i % 4 + 1}" for i in range(_N)],
        "area": np.random.rand(_N) * 1000.0,
        "GGO ID": _GGO_IDS,
        **{s: np.random.rand(_N) for s in _SCORE_COLS},
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
    n_ct = obs["celltype_broad"].nunique()
    adata.uns["celltype_broad_colors"] = ["#1f77b4"] * n_ct
    adata.uns["pathology_colors"] = _PATHOLOGY_COLORS
    adata.uns["radio_colors"] = ["#ff7f0e"] * 4
    adata.uns["Group_colors"] = ["#99B898", "#FECEAB", "#E84A5F", "#2A363B"]
    # Index = GGO IDs so that roi_pca_plot_group merge works
    adata.obs.index = pd.Index(_GGO_IDS)
    return adata


# ---------------------------------------------------------------------------
# Mock factories
# ---------------------------------------------------------------------------

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


def _make_sc_mock(read_side_effect=None) -> MagicMock:
    import scanpy as _real_sc

    sc = MagicMock()
    sc.logging = _real_sc.logging
    sc.settings = _real_sc.settings
    sc.pp.subsample = MagicMock(return_value=_make_fake_adata())
    sc.tl.score_genes = MagicMock(return_value=None)
    if read_side_effect is not None:
        sc.read = MagicMock(side_effect=read_side_effect)
    else:
        sc.read = MagicMock(return_value=_make_fake_adata())
    return sc


# ---------------------------------------------------------------------------
# Session-scoped base setup: install mocks that persist for the whole session
# (individual tests may re-mock for script-specific needs)
# ---------------------------------------------------------------------------

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
    """Install a scanpy mock for every test.

    For roi_pca_plot_group the mock uses a side_effect that returns
    a no-Group AnnData on the first sc.read call and a with-Group AnnData
    on the second call, so that the obs.merge inside that script succeeds.
    """
    test_name = request.node.name

    if test_name == "test_import_roi_pca_plot_group":
        # First sc.read call: density_all (needs NO Group column to avoid merge collision)
        # Second sc.read call: a (needs Group column)
        _call_count = [0]
        _no_group = _make_fake_adata(include_group=False)
        _with_group = _make_fake_adata(include_group=True)

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
