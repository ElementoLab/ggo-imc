import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

adata_dict = dict()

for p in ['PANEL_G', 'PANEL_H']:
    adata_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])
    adata_dict[p].obs['Sample ID'] = adata_dict[p].obs['GGO ID'].str[:6]


pdf = pd.read_csv('metadata/mutation_piv.csv', index_col = 0)
genes_ = ['EGFR', 'KRAS', 'TP53']
pdf = pdf.rename(columns = {g: f'WES_{g}' for g in genes_})


for p in ['PANEL_G', 'PANEL_H']:
    adata = adata_dict[p]
    obs = adata.obs.merge(pdf, left_on = 'Sample ID', right_index = True, how="left")
    idx = adata.obs.index
    adata.obs = obs
    adata.obs.index = idx

    # Molecular Prediction
    adata.obs.loc[adata.obs['Molecular'].str.lower() == 'not performed', 'Molecular'] = np.nan
    adata.obs.loc[adata.obs['Molecular'].str.lower() == 'nor performed', 'Molecular'] = np.nan

    # mutational status
    adata.obs['Molecular_EGFR'] = adata.obs['Molecular'].str.upper().str.contains('EGFR')
    adata.obs['Molecular_KRAS'] = adata.obs['Molecular'].str.upper().str.contains('KRAS')
    adata.obs['Molecular_TP53'] = adata.obs['Molecular'].str.upper().str.contains('TP53')

    # idx = adata.obs[['EGFR', 'KRAS', 'TP53']].isna().any(axis = 1)


conditions = [f'Molecular_{x}' for x in ['EGFR', 'KRAS', 'TP53']] + [f'WES_{x}' for x in ['EGFR', 'KRAS', 'TP53']]

celltype = 'celltype_broad'
for p in ['PANEL_G', 'PANEL_H']:
    tumor = adata_dict[p][adata_dict[p].obs['radio']!='N']
    # computing celltype density
    density = imc.tl.celltype_density(
        tumor,
        celltype = celltype,
        condition_keys = conditions)

    # per condition
    for cond in conditions:
        
        density.obs[cond] = density.obs[cond].replace({False:'WT', True:'MUT'})
        density.obs[cond] = pd.Categorical(density.obs[cond], categories = ['WT', 'MUT'])

        # statistical testing
        imc.tl.grouped_mwu_test(
            density,
            condition_keys = [cond]
        )

        # produce figures
        for pval_form in ['star', 'sci_notation']:
            imc.pl.plot_mwu(
                density,
                kind = 'violin',
                save_dir=f'figures/{p}/mutation_density/',
                pval_form=pval_form
            )

