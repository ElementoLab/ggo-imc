"""Figure 1: cell type heatmap and UMAP plots for PANEL_G and PANEL_H."""

import imc_analysis as imc
import scanpy as sc
import matplotlib.pyplot as plt

from utils import load_config, load_panels, save_figure, ensure_dir


if __name__ == '__main__':
    metadata = load_config()
    adata_dict = load_panels(metadata)

    ensure_dir('figures/figure1/')

    print('Plot celltype heatmap')
    for p in adata_dict:
        imc.pl.celltype_heatmap(
            adata_dict[p],
            cluster_ids=['celltype_broad'],
            var_names=metadata[p]['var_celltype_groups'],
            out_dir=f'figures/figure1/celltype/{p}/'
        )

    print('Plot UMAP')
    for p in adata_dict:
        adata = sc.pp.subsample(
            adata_dict[p],
            copy=True,
            n_obs=100000,
            random_state=0,
        )

        celltype = 'celltype_broad'

        fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)
        sc.pl.umap(adata, color=celltype, frameon=False, show=False)

        ensure_dir(f'figures/figure1/celltype/{p}/')
        plt.tight_layout()
        save_figure(f'figures/figure1/celltype/{p}/umap.pdf')
