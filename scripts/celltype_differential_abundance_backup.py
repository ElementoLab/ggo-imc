import imc_analysis as imc
import os
from pathlib import Path
from tqdm import tqdm

import scanpy as sc
import anndata

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=20)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

from load_yaml import load
metadata = load('metadata/ggo_config.yml')

adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])
    adata_dict[p].obs['PID'] = adata_dict[p].obs['description'].str[:6]
    
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
    # adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)
    adata_dict[p].obs['pathology'] = adata_dict[p].obs['pathology'].astype(str).str.replace('Normal', 'N')
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['N', 'AIS', 'MIA', 'IAC'], ordered = True)
    
    adata_dict[p].uns['radio_color'] = metadata['radio_color']
    adata_dict[p].uns['pathology_color'] = metadata['pathology_color']
    
    del adata_dict[p].obs['ROI_area']

conditions = metadata['CONDITIONS']
density_dict = dict()
# overall celltype
for p in ['PANEL_G', 'PANEL_H']:
    
    adata = adata_dict[p]
    adata.obs['all_cells'] = 'all_cells'
    
    density = imc.tl.celltype_density(
        adata,
        celltype = 'all_cells',#'celltype_broad',
        condition_keys = conditions)
    density_dict[p] = density
    for cond in ['pathology', 'radio']:
        imc.tl.grouped_mwu_test(
            density,
            condition_keys = [cond]
        )

        for pval_form in ['star', 'sci_notation']:
            imc.pl.plot_mwu(
                density,
                save_dir=f'figures/{p}/all_cell_density/',
                pval_form=pval_form
            )


# conditions = ['radio']
for p in ['PANEL_G', 'PANEL_H']:
    for celltype in metadata['CELLTYPES']:
        
        density = imc.tl.celltype_density(
            adata_dict[p],
            celltype = celltype,
            condition_keys = conditions)

        density.write(metadata[p]['AnnData'][f'roi_{celltype}_name'])

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
                plot_mwu(
                    density,
                    save_dir=f'figures/{p}/differential_{celltype}_density/',
                    palette=palette,
                    pval_form=pval_form
                )


for p in ['PANEL_G', 'PANEL_H']:
    for celltype in metadata['CELLTYPES']:
        adata = adata_dict[p].copy()
        adata.obs['TN'] = 'T'
        adata.obs.loc[adata.obs['radio'] == 'N', 'TN'] = 'N'
        adata.obs['P_ID'] = adata.obs['PID'].astype(str) + '_' + adata.obs['TN']
        
        density = imc.tl.patient_density(
            adata,
            celltype_key = celltype,
            condition_keys = conditions,
            patient_key = 'P_ID')

        density.write(metadata[p]['AnnData'][f'patient_{celltype}_name'])

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
                plot_mwu(
                    density,
                    save_dir=f'figures/{p}/patient/differential_{celltype}_density/',
                    palette=palette,
                    pval_form=pval_form
                )





# Exome Comparison
metadata = load('metadata/ggo_config.yml')

# adata = adata[:,metadata['celltype_broad_orders']]
# for cond in ['pathology', 'radio', 'EGFR', 'KRAS', 'RBM10', 'TP53', 'Smoking Status', 'Gender', 'Race', 'Group']:
#     imc.tl.grouped_mwu_test(
#         adata,
#         condition_keys = [cond]
#     )
#     if f'{cond}_color' in metadata:
#         palette = metadata[f'{cond}_color']
#     else:
#         palette = None #'tab10'

#     for pval_form in ['star', 'sci_notation']:
#         plot_mwu(
#             adata,
#             save_dir=f'figures/patient/density/',
#             palette=palette,
#             pval_form=pval_form
#         )

print(f"Reading {metadata['patient_celltype_broad_clustered']}...")
exome = sc.read(metadata['patient_celltype_broad_clustered'])
df = exome.obs[['EGFR', 'KRAS', 'RBM10', 'TP53']]
conditions = ['EGFR', 'KRAS', 'RBM10', 'TP53', 'Smoking Status', 'Gender', 'Race']
for p in ['PANEL_G', 'PANEL_H']:
    for celltype in metadata['CELLTYPES']:
        adata = adata_dict[p].copy()
        adata = adata[adata.obs['pathology'] != 'N']
        joint = adata.obs.merge(df.reset_index(), left_on='PID', right_on = 'index', how = 'left')
        for feature in ['EGFR', 'KRAS', 'RBM10', 'TP53']:
            adata.obs[feature] = pd.Categorical(joint[feature].tolist(), categories = ['WT', 'MUT'])

        density = imc.tl.celltype_density(
            adata,
            celltype = celltype,
            condition_keys = conditions)

        for cond in conditions:
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
                    save_dir=f'figures/{p}/differential_{celltype}_density/',
                    palette=palette,
                    pval_form=pval_form
                )

df = exome.obs[['EGFR', 'KRAS', 'RBM10', 'TP53']]
conditions = ['EGFR', 'KRAS', 'RBM10', 'TP53', 'Smoking Status', 'Gender', 'Race']
for p in ['PANEL_G', 'PANEL_H']:
    for celltype in metadata['CELLTYPES']:
        adata = adata_dict[p].copy()
        adata = adata[adata.obs['pathology'] == 'N']
        joint = adata.obs.merge(df.reset_index(), left_on='PID', right_on = 'index', how = 'left')
        for feature in ['EGFR', 'KRAS', 'RBM10', 'TP53']:
            adata.obs[feature] = pd.Categorical(joint[feature].tolist(), categories = ['WT', 'MUT'])

        density = imc.tl.celltype_density(
            adata,
            celltype = celltype,
            condition_keys = conditions)

        for cond in conditions:
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
                    save_dir=f'figures/{p}/differential_{celltype}_density/Normal/',
                    palette=palette,
                    pval_form=pval_form
                )


'''
repeat the same, but for samples that are N, PNS-AIS, PS-MIA, S-IAC
'''
# for p in ['PANEL_G', 'PANEL_H']:
#     for celltype in metadata['CELLTYPES']:
        
#         density = imc.tl.celltype_density(
#             adata_dict[p],
#             celltype = celltype,
#             condition_keys = metadata['CONDITIONS'])

#         filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]
        
#         density = anndata.concat([
#             density[(density.obs['pathology'] == s['pathology']) & (density.obs['radio'] == s['radio'])] for s in filt
#         ])
#         density.obs['radio'] = pd.Categorical(density.obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered = True)
#         density.obs['pathology'] = pd.Categorical(density.obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)

#         density = density[:,~density.var.index.isin(['Low Expr.', 'Other'])]
        
#         for cond in conditions:
#             imc.tl.grouped_mwu_test(
#                 density,
#                 condition_keys = [cond]
#             )
#             for pval_form in ['star', 'sci_notation']:
#                 if f'{cond}_color' in metadata:
#                     palette = metadata[f'{cond}_color']
#                 else:
#                     palette = None #'tab10'

#                 imc.pl.plot_mwu(
#                     density,
#                     save_dir=f'figures/{p}/differential_{celltype}_density_filtered/',
#                     palette=palette,
#                     pval_form=pval_form
#                 )

# conditions = metadata['CONDITIONS'] + ['PID','description']

'''
Perform three comparisons at a patient level

1. Between tumor samples
2. Between normal and tumor samples combined. Here we take an approach saying that 
    normal samples and tumor samples are completely independent identity, and not from
    the same patient
3. Tumor - normal differentials per patient
'''

# logger = imc.logging.logger
# logger.setLevel(logging.DEBUG)  # Set the log level

# conditions = metadata['CONDITIONS'] + ['Smoking Status', 'Race']
# for p in ['PANEL_G', 'PANEL_H']:
#     for celltype in metadata['CELLTYPES']:
#         adata = adata_dict[p].copy()
        
#         adata.obs['TN'] = 'T'
#         adata.obs.loc[adata.obs['radio'] == 'N', 'TN'] = 'N'
#         adata.obs['P_ID'] = adata.obs['PID'].astype(str) + '_' + adata.obs['TN']

#         p_density = imc.tl.patient_density(
#             adata, # d[k],
#             patient_key = 'P_ID',
#             celltype_key = celltype,
#             condition_keys = conditions
#         )

#         p_density.write(metadata[p]['AnnData'][f'patient_{celltype}_name'])

#         n2 = p_density[p_density.obs['radio'] == 'N'].copy()
#         t2 = p_density[p_density.obs['radio'] != 'N'].copy()
#         n2.obs.index = n2.obs.index.str[:6]
#         t2.obs.index = t2.obs.index.str[:6]
#         intersection = n2.obs.index.intersection(t2.obs.index)
#         n2 = n2[intersection]
#         t2 = t2[intersection]
#         t3 = t2.copy()
#         t3.X = t3.X - n2.X

#         comps = {
#             'tumor_only': p_density[p_density.obs['radio'] != 'N'],
#             'combined': p_density,
#             'differential': t3, # need to do differential separately with shared index without N and T
#         }

#         for c in comps:
#             p_density_ = comps[c]
#             for cond in conditions:
#                 p_density_.obs[cond] = pd.Categorical(p_density_.obs[cond])
#                 imc.tl.grouped_mwu_test(
#                     p_density_,
#                     condition_keys = [cond]
#                 )
#                 for pval_form in ['star', 'sci_notation']:
#                     if f'{cond}_color' in metadata:
#                         palette = metadata[f'{cond}_color']
#                         if cond == 'differential':
#                             palette = palette[1:]
#                     else:
#                         palette = None #'tab10'

#                     imc.pl.plot_mwu(
#                         p_density_,
#                         palette = palette,
#                         pval_form = pval_form,
#                         save_dir = f'figures/{p}/patient/{c}/{celltype}/'
#                     )


# EMT extent
p = 'PANEL_G'
celltype = 'celltype_broad'
adata = adata_dict[p].copy()
adata.obs['TN'] = 'T'
adata.obs.loc[adata.obs['radio'] == 'N', 'TN'] = 'N'
adata.obs['P_ID'] = adata.obs['PID'].astype(str) + '_' + adata.obs['TN']

p_density = imc.tl.patient_density(
    adata,
    patient_key = 'P_ID',
    celltype_key = celltype,
    condition_keys = metadata['CONDITIONS']
)
emt = anndata.AnnData(
    X = np.concatenate([
            p_density[:,'Mesen.-like'].X / p_density[:,'Epi.-like'].X,
            p_density[:,'Mesen.-like'].X + p_density[:,'Epi.-like'].X,
        ], axis = 1),
    obs = p_density.obs)
emt.var.index = ['M/E ratio', 'M+E count']

for cond in metadata['CONDITIONS']:
    emt.obs[cond] = pd.Categorical(emt.obs[cond])
    imc.tl.grouped_mwu_test(
        emt,
        condition_keys = [cond]
    )
    for pval_form in ['star', 'sci_notation']:
        if f'{cond}_color' in metadata:
            palette = metadata[f'{cond}_color']
            if cond == 'differential':
                palette = palette[1:]
        else:
            palette = None #'tab10'

        imc.pl.plot_mwu(
            emt,
            palette = palette,
            pval_form = pval_form,
            save_dir = f'figures/{p}/patient/emt/'
        )



# # patient feature selection
# # measure mean intensity on patient level

# metadata = load('metadata/ggo_config.yml')

# adata_dict = dict()
# for p in ['PANEL_G', 'PANEL_H']:
#     print(f"Reading {metadata[p]['AnnData']['phenotyped_file_name']}...")
#     adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])


# conditions = metadata['CONDITIONS'] + ['PID', 'Race', 'Smoking', 'Packyrs', 'Age', 'Molecular']
# for p in ['PANEL_G', 'PANEL_H']:

#     adata = adata_dict[p]
#     adata.obs['PID'] = adata.obs['description'].str[:6]

#     n = adata[adata.obs['radio']=='N']
#     t = adata[adata.obs['radio']!='N']

#     n_ = grouped_obs_mean(n, group_key = 'PID')
#     n_obs = n.obs[conditions].drop_duplicates().sort_values('PID').set_index('PID')
#     n_obs.index.name = None
#     n_data = anndata.AnnData(X = n_, obs = n_obs)
#     n_data.obs['group'] = 'N'
#     n_data = n_data[:,~n_data.var.index.str.contains('EMPTY')]
#     n_data.obs.index = n_data.obs.index + '_N'
    
#     t_ = grouped_obs_mean(t, group_key = 'PID')
#     t_obs = t.obs[conditions].drop_duplicates().sort_values('PID').set_index('PID')
#     t_obs.index.name = None
#     t_data = anndata.AnnData(X = t_, obs = t_obs)
#     t_data.obs['group'] = 'T'
#     t_data = t_data[:,~t_data.var.index.str.contains('EMPTY')]
#     t_data.obs.index = t_data.obs.index + '_T'

#     n_data.write(metadata[p]['AnnData']['patient_normal_expression_name'])
#     t_data.write(metadata[p]['AnnData']['patient_tumor_expression_name'])



# normal.var['group'] = 'N'
# tumor.var['group'] = 'T'
# adata = anndata.concat([normal, tumor], axis = 1, join = 'outer')
# adata.obs['pathology'] = tumor.obs['pathology']
# adata.obs['radio'] = tumor.obs['radio']
# adata.obs['ROI_area_T'] = tumor.obs['ROI_area']
# adata.obs['ROI_area_N'] = normal.obs['ROI_area']

# adata.write('results/Patient_cellcount_all.h5ad')

