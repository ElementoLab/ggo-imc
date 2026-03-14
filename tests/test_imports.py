"""Smoke tests: verify all pipeline scripts import without error."""


def test_import_utils():
    import utils  # noqa: F401


def test_import_celltype_heatmap_info():
    import celltype_heatmap_info  # noqa: F401


def test_import_roi_pca_plot():
    import roi_pca_plot  # noqa: F401


def test_import_celltype_differential_abundance():
    import celltype_differential_abundance  # noqa: F401


def test_import_t_cell_analysis():
    import t_cell_analysis  # noqa: F401


def test_import_myeloid_analysis():
    import myeloid_analysis  # noqa: F401


def test_import_epithelial_characterization():
    import epithelial_characterization  # noqa: F401


def test_import_ue_analysis():
    import ue_analysis  # noqa: F401


def test_import_roi_pca_plot_group():
    import roi_pca_plot_group  # noqa: F401


def test_import_patient_group():
    import patient_group  # noqa: F401
