
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
# matplotlib.use('Agg')

metadata = load('metadata/ggo_config.yml')

print(f"Reading {metadata['patient_celltype_broad_clustered']}...")
adata = sc.read(metadata['patient_celltype_broad_clustered'])
adata = adata[:,metadata['celltype_broad_orders']]
df = adata.to_df() / 1000
prop = df / np.array(df.sum(axis = 1))[:,None] * 100

# cell type density
fig, ax = plt.subplots(1,1,figsize = (8,2.5), dpi = 300)
df.sort_values(['Epi.-like', 'Fib.', 'Mac.', 'CD8 T'], ascending = False).plot(kind='bar', stacked=True, color=metadata['celltype_broad_colors'], ax = ax, width = 0.9)
ax.legend().remove()
ax.set_xticks([])
ax.set_xlabel('')
ax.set_ylabel(r'Cell Density (x1000 cells / mm$^2$)')
plt.yticks(rotation = 90)
sns.despine()
plt.savefig('figures/celltype_density.pdf')
plt.close()

# cell type density normalized to percentage
fig, ax = plt.subplots(1,1,figsize = (8,2.5), dpi = 300)
prop.sort_values(['Epi.-like', 'Fib.', 'Mac.', 'CD8 T'], ascending = False).plot(kind='bar', stacked=True, color=metadata['celltype_broad_colors'], ax = ax, width = 0.9)
ax.legend().remove()
ax.set_xticks([])
ax.set_xlabel('')
plt.yticks(rotation = 90)
sns.despine()
ax.set_ylabel('Cell Proportion')
plt.savefig('figures/celltype_proportion.pdf')
plt.close()

heatmap_data = adata[prop.sort_values(['Epi.-like', 'Fib.', 'Mac.', 'CD8 T'], ascending = False).index].obs[['pathology', 'radio', 'Group']]
heatmap_data.to_csv('results/heatmap_data.csv')


from sklearn.metrics import silhouette_score, silhouette_samples, calinski_harabasz_score
score_data = adata.copy()
score_data.X = prop

# score_dict = dict()
# for g in ['pathology', 'radio', 'Group']:
#     score_data.obs[f'silhouette_score_{g}'] = silhouette_samples(score_data.X, score_data.obs[g].astype(str).fillna(''))
#     score_dict[g] = silhouette_score(score_data.X, score_data.obs[g].astype(str).fillna(''))
# print(score_dict)
# print(score_data.obs.describe())


score_dict = dict()
for g in ['pathology', 'radio', 'Group']:
    score_dict[g] = calinski_harabasz_score(score_data.X, score_data.obs[g].astype(str).fillna(''))    
vrc = pd.DataFrame(score_dict.values(), index = score_dict.keys(), columns = ['Variance Ratio Criterion'])
vrc = vrc.loc[['Group','pathology', 'radio']]
fig, ax = plt.subplots(1,1,figsize = (3,2), dpi = 300)
vrc.round(2).plot(kind='barh', ax = ax, width = 0.8, color = 'lightgray')
for container in ax.containers:
    ax.bar_label(container)
ax.legend().remove()
ax.set_xlabel('Variance Ratio Criterion')
ax.set_ylabel('Stratification Method')
sns.despine()
plt.tight_layout()
plt.savefig('figures/stratification_vrc.pdf')
plt.close()
