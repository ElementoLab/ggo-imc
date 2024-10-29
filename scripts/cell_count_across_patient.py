import scanpy as sc
from glob import glob
import tifffile
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)

import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False


def get_feature_prop(
    adata,
    sample_key = 'sample_y',
    feature_key = 'celltype',
    proportion = False,
    per_thousand = True,
    drop_unknown = False,
):
    count = adata.obs.groupby([sample_key, feature_key]).count()['sample'].reset_index()
    count = count.pivot(index = sample_key, columns = feature_key, values= 'sample')
    if proportion:
        prop = count / count.sum(axis = 1)[:, np.newaxis]
    elif per_thousand:
        prop = count / 1000
    prop = prop.merge(adata.obs[[sample_key, 'ggo', 'sample status']].drop_duplicates().set_index(sample_key), left_index = True, right_index = True)
    
    if drop_unknown:
        other_key = prop.columns[~(prop.columns.str.contains('Other') | prop.columns.str.contains('other'))]
        prop = prop[other_key]
        
    columns = prop.columns
    epi_columns = columns.str.contains('Tumor') | columns.str.contains('Alveolar') | columns.str.contains('AT')
    stromal_columns = columns.str.contains('Endothelial') | columns.str.contains('Fibroblast') | columns.str.contains('Mesenchymal')
    immune_columns = ~(columns.isin(columns[epi_columns]) | columns.isin(columns[stromal_columns]) | columns.str.contains('Other'))

    prop_ = prop.copy()
    prop_['epi'] = prop[prop.columns[epi_columns]].sum(axis = 1)
    prop_['stro'] = prop[prop.columns[stromal_columns]].sum(axis = 1)
    prop_['imm'] = prop[prop.columns[immune_columns]].sum(axis = 1)

    prop_ = prop_.sort_values(['sample status', 'ggo', 'epi', 'stro', 'imm'])
    order = prop.columns[epi_columns].tolist() +  prop.columns[stromal_columns].tolist() + prop.columns[immune_columns].tolist()[:-2] + columns[columns.str.contains('Other')].tolist() + prop.columns[immune_columns].tolist()[-2:]

    del prop_['epi']
    del prop_['stro']
    del prop_['imm']

    pg_prop = prop_[order]
    return pg_prop

def get_feature_prop(
    adata,
    sample_key = 'sample_y',
    feature_key = 'celltype',
    proportion = False,
    per_thousand = True,
    drop_unknown = False,
):
    count = adata.obs.groupby([sample_key, feature_key]).count()['sample'].reset_index()
    count = count.pivot(index = sample_key, columns = feature_key, values= 'sample')
    if proportion:
        prop = count / count.sum(axis = 1)[:, np.newaxis]
    elif per_thousand:
        prop = count / 1000
    prop = prop.merge(adata.obs[[sample_key, 'ggo', 'sample status']].drop_duplicates().set_index(sample_key), left_index = True, right_index = True)
    
    if drop_unknown:
        other_key = prop.columns[~(prop.columns.str.contains('Other') | prop.columns.str.contains('other'))]
        prop = prop[other_key]
        
    columns = prop.columns
    epi_columns = columns.str.contains('Tumor') | columns.str.contains('Alveolar') | columns.str.contains('AT')
    stromal_columns = columns.str.contains('Endothelial') | columns.str.contains('Fibroblast') | columns.str.contains('Mesenchymal')
    immune_columns = ~(columns.isin(columns[epi_columns]) | columns.isin(columns[stromal_columns]) | columns.str.contains('Other'))

    prop_ = prop.copy()
    prop_['epi'] = prop[prop.columns[epi_columns]].sum(axis = 1)
    prop_['stro'] = prop[prop.columns[stromal_columns]].sum(axis = 1)
    prop_['imm'] = prop[prop.columns[immune_columns]].sum(axis = 1)

    prop_ = prop_.sort_values(['sample status', 'ggo', 'epi', 'stro', 'imm'])
    order = prop.columns[epi_columns].tolist() +  prop.columns[stromal_columns].tolist() + prop.columns[immune_columns].tolist()[:-2] + columns[columns.str.contains('Other')].tolist() + prop.columns[immune_columns].tolist()[-2:]

    del prop_['epi']
    del prop_['stro']
    del prop_['imm']

    pg_prop = prop_[order]
    return pg_prop

adata_dict = dict()
adata_dict['PANEL_G'] = sc.read(metadata['panel_g_phenotyped_file_name'])
adata_dict['PANEL_H'] = sc.read(metadata['panel_h_phenotyped_file_name'])

for panel in adata_dict:
	adata = adata_dict[panel]
	prop = get_feature_prop(adata, drop_unknown=False)

	fig, ax = plt.subplots(2,1,figsize = (11,6), dpi = 300)

	colors = {celltype: color for celltype, color in zip(adata.obs['celltype'].cat.categories, adata.uns['celltype_colors'])}
	colors = [colors[col] for col in prop.columns[:-2]]

	prop[prop['Radiology'] == 'Normal'].plot(kind='bar', stacked=True, color=colors, ax = ax[0], width = 0.9)
	ax[0].legend().remove()
	ax[0].set_xticks([])
	ax[0].set_xlabel('')
	ax[0].set_ylabel('Normal Cell Counts (x1,000)')

	prop[prop['Radiology'] != 'Normal'].plot(kind='bar', stacked=True, color=colors, ax = ax[1], width = 0.9)
	ax[1].legend().remove()
	ax[1].set_xticks([])
	ax[1].set_xlabel('')
	ax[1].set_ylabel('Tumor Cell Counts (x1,000)')
	plt.tight_layout()
	plt.savefig(f'figures/celltype_proportion_{PANEL}.pdf')
	plt.show()