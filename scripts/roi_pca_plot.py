"""Figure 1: ROI PCA archetype plot (all pathologies combined)."""

import imc_analysis as imc
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)

from utils import load_config, save_figure


if __name__ == '__main__':
    metadata = load_config()

    density_all = sc.read(
        metadata['PCA_ANNDATA_NAME'],
        backup_url=metadata['PCA_URL'],
    )

    # Perform PCA — flip x so that normal tissue is on the negative side
    sc.pp.scale(density_all)
    sc.pp.pca(density_all)
    density_all.obsm['X_pca'][:, 0] = -density_all.obsm['X_pca'][:, 0]

    fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=300)
    sc.pl.pca(
        density_all,
        color='pathology',
        size=200,
        add_outline=True,
        outline_width=(0.15, 0.),
        annotate_var_explained=True,
        title='',
        ax=ax,
        show=False,
    )
    sns.despine()
    plt.tight_layout()
    save_figure('figures/figure1/roi_pca_all_combined.pdf')
