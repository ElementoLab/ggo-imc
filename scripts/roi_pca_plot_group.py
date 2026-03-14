"""Figure 5: ROI PCA coloured by patient risk group with covariance ellipses."""

import imc_analysis as imc
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse

sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)

from utils import load_config, save_figure, ensure_dir


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

    # Load patient group annotations and merge onto the PCA AnnData
    a = sc.read(
        metadata['patient_celltype_broad_clustered'],
        backup_url=metadata['patient_group_url'],
    )

    tmp = density_all.obs.merge(a.obs[['Group']], left_on='GGO ID', right_index=True)
    tmp['Group'] = tmp['Group'].astype(str)
    tmp.loc[tmp['pathology'] == 'N', 'Group'] = 'N'
    tmp['Group'] = pd.Categorical(
        tmp['Group'],
        categories=['N', 'Group1', 'Group2', 'Group3', 'Group4'],
        ordered=True,
    )

    da = density_all[tmp.index]
    da.obs = tmp
    da.uns['Group_colors'] = ['#2AC674', '#8DAC91', '#F9C8AB', '#DB4962', '#252F32']

    # PCA scatter
    fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)
    sc.pl.pca(
        da,
        color='Group',
        size=200,
        add_outline=True,
        outline_width=(0.15, 0.),
        annotate_var_explained=True,
        title='',
        legend_loc=None,
        ax=ax,
        show=False,
    )

    # Overlay mean + covariance ellipse per group
    for i, group in enumerate(da.obs['Group'].cat.categories):
        da_group = da[da.obs['Group'] == group]
        mean = da_group.obsm['X_pca'][:, :2].mean(axis=0)
        cov = np.cov(da_group.obsm['X_pca'][:, :2].T)

        eigenvalues, eigenvectors = np.linalg.eig(cov)
        angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
        width, height = 2 * np.sqrt(eigenvalues)
        ellipse = Ellipse(
            xy=mean, width=width, height=height, angle=angle,
            facecolor=da.uns['Group_colors'][i], edgecolor=None, alpha=0.3,
        )
        plt.gca().add_patch(ellipse)

    ensure_dir('figures/figure5/')
    sns.despine()
    plt.tight_layout()
    save_figure('figures/figure5/roi_pca_group.pdf')
