
import os
from glob import glob
from pathlib import Path
import tifffile
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt

import anndata
import warnings
warnings.simplefilter("ignore", UserWarning)

import imc_analysis as imc
from scripts.load_yaml import load

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

metadata = load('metadata/ggo_config.yml')

p = 'PANEL_H'
print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
adata = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])
adata.obs['PID'] = adata.obs['description'].str[:6]
adata.obs['radio'] = pd.Categorical(adata.obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
# adata.obs['pathology'] = pd.Categorical(adata.obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)
adata.obs['pathology'] = adata.obs['pathology'].astype(str).str.replace('Normal', 'N')
adata.obs['pathology'] = pd.Categorical(adata.obs['pathology'], categories = ['N', 'AIS', 'MIA', 'IAC'], ordered = True)

#epithelial_celltypes = ['Tumor-like (RAGE+)', 'Tumor-like (SFTPC+)', 'Tumor-like', 'Epi.-like (RAGE+)', 'Epi.-like (SFTPC+)']
celltype_map = {
	'Tumor-like (RAGE+)':'RAGE+',
	'Tumor-like (SFTPC+)':'SFTPC+',
	'Tumor-like':'PanCK+',
	'Epi.-like (RAGE+)':'RAGE+',
	'Epi.-like (SFTPC+)':'SFTPC+',
}
adata = adata[adata.obs['celltype'].isin(celltype_map.keys())]
adata.obs['phenotype'] = adata.obs['celltype'].replace(celltype_map)
adata.obs['phenotype'] = pd.Categorical(adata.obs['phenotype'], categories = ['PanCK+', 'RAGE+', 'SFTPC+'], ordered = True)
path_phenotype = adata.obs.groupby(['phenotype', 'pathology']).count()['sample'].unstack()
path_phenotype = path_phenotype / np.array(path_phenotype.sum(axis = 0))[None,:]

fig, ax = plt.subplots(1,1,figsize = (3,2), dpi = 300)
path_phenotype.T.reset_index().plot(
	x='pathology',
	kind='bar',
	stacked=True,
	title='',
	ax = ax)

df_rel = path_phenotype.T
df_pos = df_rel.cumsum(1)-df_rel/2

for i, x in enumerate(df_rel.index):
	for y in df_rel.columns:
		# print(x, df_pos.loc[x,y], str(np.round(df_rel.loc[x,y] *100, 1)))
		plt.text(i, df_pos.loc[x,y], str(np.round(df_rel.loc[x,y] *100, 1)) + '%', va = 'center', ha = 'center')
sns.despine()
ax.legend().remove()
ax.set_xlabel('')
ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
plt.tight_layout()
plt.savefig('figures/epithelial_proportion_pathology.pdf')
plt.close()

radio_phenotype = adata.obs.groupby(['phenotype', 'radio']).count()['sample'].unstack()
radio_phenotype = radio_phenotype / np.array(radio_phenotype.sum(axis = 0))[None,:]

fig, ax = plt.subplots(1,1,figsize = (3,2), dpi = 300)
radio_phenotype.T.reset_index().plot(
	x='radio',
	kind='bar',
	stacked=True,
	title='',
	ax = ax)

df_rel = radio_phenotype.T
df_pos = df_rel.cumsum(1)-df_rel/2

for i, x in enumerate(df_rel.index):
	for y in df_rel.columns:
		# print(x, df_pos.loc[x,y], str(np.round(df_rel.loc[x,y] *100, 1)))
		plt.text(i, df_pos.loc[x,y], str(np.round(df_rel.loc[x,y] *100, 1)) + '%', va = 'center', ha = 'center')
sns.despine()
ax.legend().remove()
ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
ax.set_xlabel('')
plt.tight_layout()
plt.savefig('figures/epithelial_proportion_radiology.pdf')
plt.close()

