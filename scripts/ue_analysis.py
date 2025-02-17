import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import os

import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc
from tqdm import tqdm

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

adata_dict = dict()

for p in ['PANEL_G', 'PANEL_H']:
    adata_dict[p] = sc.read(
        metadata[p]['AnnData']['utag_labeled_name'],
        backup_url = metadata[p]['AnnData']['utag_url']
    )

# Plot Celltype Composition
import seaborn as sns
celltype = 'celltype_broad'
uE = 'uE_broad'
for p in ['PANEL_G', 'PANEL_H']:
    count = adata_dict[p].obs.groupby([celltype, uE]).count()
    count = count / 1000
    g = sns.FacetGrid(
        data=count["sample"].reset_index(),
        hue=celltype,
        col=uE,
        col_wrap=2,
        aspect=1.3,
        palette=adata_dict[p].uns['celltype_broad_colors'],
        sharex=False,
        height = 3,
    )

    g.map(sns.barplot, "sample",  celltype)

    # iterate over axes of FacetGrid
    for i, ax in enumerate(g.axes.flat):
        labels = ax.get_yticklabels()  # get x labels
        #ax.set_xticklabels(labels, rotation=90, fontsize=10)  # set new labels
        titles = ax.get_title()
        titles = titles.replace(f"{uE} = ", "")
        ax.set_title(titles, fontsize=12)
        ax.set_ylabel("", fontsize=12)
        ax.set_xlabel("")
    plt.tight_layout()
    plt.xlabel('x1000 cells')
    plt.savefig(f"figures/figure4/{p}_{uE}_{celltype}_composition.pdf")


import seaborn as sns

p = 'PANEL_G'
uE = 'uE_broad'
density = imc.tl.celltype_density(
    adata_dict[p],
    celltype = uE,
    condition_keys = ['pathology', 'radio'])

# for cond in ['pathology', 'radio']:
#     imc.tl.grouped_mwu_test(
#         density,
#         condition_keys = [cond]
#     )

#     for pval_form in ['star', 'sci_notation']:
#         imc.pl.plot_mwu(
#             density,
#             save_dir=f'figures/figure4/{p}/{uE}_density/',
#             kind = 'box-line',
#             pval_form=pval_form
#         )

# Area based comparison
area = adata_dict[p].obs.groupby(['roi',uE])[['area']].sum().pivot_table(index = 'roi', columns = uE)
area = area / 1e6
area.columns = area.columns.droplevel(0)

area_adata = anndata.AnnData(X = area, obs = density.obs, var = density.var)
imc.tl.grouped_mwu_test(
    area_adata,
    condition_keys = ['radio', 'pathology']
)

area_adata.uns['radio_colors'] = adata_dict[p].uns['radio_colors']
area_adata.uns['pathology_colors'] = adata_dict[p].uns['pathology_colors']

for pval_form in ['star', 'sci_notation']:
    imc.pl.plot_mwu(
        area_adata,
        save_dir=f'figures/figure4/{p}/{uE}_area/',
        pval_form=pval_form,
        kind = 'violin',
    )


functional_markers = ['HLADR', 'GranzymeB', 'Tbet', 'CD163',
'ICOS', 'CCR4', 'LAG3', 'CD4',
'CCR7', 'CD28', 'CD11b', 
'CD45RA', 'CD206', 'Ki67', 'CD57', 'CD27', '41BB', 'CD45RO',
'TIM3', 'CD25', 'CTLA4']

#for p in adata_dict:
celltype = 'uE_broad'
p = 'PANEL_G'
# p = 'PANEL_H'
for feat in ['pathology', 'radio']:
    for ct in tqdm(adata_dict[p].obs[celltype].unique()):

        a = adata_dict[p][adata_dict[p].obs[celltype] == ct, functional_markers].copy()
        sc.tl.rank_genes_groups(a, groupby = feat, method = 'wilcoxon', use_raw = False)
        sc.pl.rank_genes_groups_dotplot(a, n_genes = 5, values_to_plot = 'logfoldchanges', min_logfoldchange = 0.2,cmap = 'bwr', show = False, title = ct, vmax = 1, vmin = -1, dendrogram = False, use_raw = False)
        # sc.pl.rank_genes_groups_dotplot(
        #     a, n_genes = 4, values_to_plot = 'scores',
        #     cmap = 'bwr', show = False, title = ct,
        #     dendrogram = False, use_raw = False)
        
        cts = ct.replace('/','')
        path = f'figures/{p}/{celltype}/{feat}/'

        os.makedirs(path, exist_ok = True)
        plt.savefig(path + f'dotplot_log_foldchange_{p}_{feat}_{cts}.pdf', bbox_inches = 'tight')
        plt.close()

        sc.pl.rank_genes_groups_dotplot(
            a, n_genes = 4, values_to_plot = 'scores',
            cmap = 'bwr', show = False, title = ct,
            dendrogram = False, use_raw = False)
        
        cts = ct.replace('/','')
        path = f'figures/{p}/{celltype}/{feat}/'

        os.makedirs(path, exist_ok = True)
        plt.savefig(path + f'dotplot_score_{p}_{feat}_{cts}.pdf', bbox_inches = 'tight')
        plt.close()
