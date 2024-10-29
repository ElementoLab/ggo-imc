
import os
from glob import glob
from pathlib import Path
import tifffile
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

# from utag import utag
import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc

from scripts.load_yaml import load
import anndata
import warnings
warnings.simplefilter("ignore", UserWarning)

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
#matplotlib.use('Agg')

import yaml

metadata = load('metadata/ggo_config.yml')

adata_dict = dict()

for p in ['PANEL_G', 'PANEL_H']:
    adata_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])


import seaborn as sns
for p in ['PANEL_G', 'PANEL_H']:
    for celltype in ['celltype', 'celltype_broad']:
        for uE in ['uE', 'uE_broad']:
            count = adata_dict[p].obs.groupby([celltype, uE]).count()
            count = count / 1000
            g = sns.FacetGrid(
                data=count["sample"].reset_index(),
                hue=celltype,
                col=uE,
                col_wrap=2,
                aspect=1.3,
                palette="colorblind",
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
            plt.savefig(f"figures/ue/{p}_{uE}_{celltype}_composition.pdf")


import seaborn as sns

p = 'PANEL_G'
p = 'PANEL_H'
p = 'PANEL_G'
celltype = 'celltype_broad'
uE = 'uE_broad'

adata = adata_dict[p].copy()
radio = adata_dict[p].obs.groupby(['radio', uE]).count()[['sample']].reset_index().pivot_table(index = 'radio', columns = uE)
radio / np.array(radio).sum(axis=1)[:,None]

pathology = adata_dict[p].obs.groupby(['pathology', uE]).count()[['sample']].reset_index().pivot_table(index = 'pathology', columns = uE)
pathology / np.array(pathology).sum(axis=1)[:,None]

for uE in ['uE', 'uE_broad']:
    del adata_dict[p].obs['ROI_area']
    density = imc.tl.celltype_density(
        adata_dict[p],
        celltype = uE,
        condition_keys = ['pathology', 'radio'])

    # density.write(metadata[p]['AnnData'][f'roi_{uE}_name'])

    for cond in ['pathology', 'radio']:
        imc.tl.grouped_mwu_test(
            density,
            condition_keys = [cond]
        )
        if f'{cond}_color' in metadata:
            palette = metadata[f'{cond}_color']
        else:
            palette = None #'tab10'

        for pval_form in ['star', 'sci_notation']:
            imc.pl.plot_mwu(
                density,
                save_dir=f'figures/{p}/differential_{uE}_density/',
                palette=palette,
                pval_form=pval_form
            )





p = 'PANEL_G'
uE = 'uE_broad'
cond = 'pathology'
# density comparison frequency
density = imc.tl.celltype_density(
    adata_dict[p],
    celltype = uE,
    condition_keys = ['pathology'])

area = adata_dict[p].obs.groupby(['roi',uE])[['area']].sum().pivot_table(index = 'roi', columns = uE)
area = area / 1e6
area.columns = area.columns.droplevel(0)

area_adata = anndata.AnnData(X = area, obs = density.obs, var = density.var)
imc.tl.grouped_mwu_test(
    area_adata,
    condition_keys = [cond]
)
if f'{cond}_color' in metadata:
    palette = metadata[f'{cond}_color']
else:
    palette = None #'tab10'

for pval_form in ['star', 'sci_notation']:
    imc.pl.plot_mwu(
        area_adata,
        save_dir=f'figures/{p}/{uE}_area/',
        palette=palette,
        pval_form=pval_form,
        kind = 'violin',
    )




tumor_density = density[density.obs['pathology']!='N'].copy()
density_df = tumor_density.to_df()

# per ROI count lymphonet connectivity (network degree)










'''
adata_dict['PANEL_G'].var.index
['HLADR', 'GranzymeB', 'Vimentin', 'CD14', 'Tbet', 'CD163',
'ICOS', 'CCR4', 'PDL1', 'LAG3', 'CD11c', 'FoxP3', 'CD4',
'CCR7', 'CD28', 'CD11b', 'PD1',
'CD45RA', 'CD206', 'Ki67', 'CD57', 'CD3', 'CD27', '41BB', 'CD45RO',
'TIM3', 'CD25', 'CTLA4']'PD1', 'PDL1', 
'''
functional_markers = ['HLADR', 'GranzymeB', 'Tbet', 'CD163',
'ICOS', 'CCR4', 'LAG3', 'CD4',
'CCR7', 'CD28', 'CD11b', 
'CD45RA', 'CD206', 'Ki67', 'CD57', 'CD27', '41BB', 'CD45RO',
'TIM3', 'CD25', 'CTLA4']

#for p in adata_dict:
celltype = 'uE_broad'
p = 'PANEL_G'
p = 'PANEL_H'
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
