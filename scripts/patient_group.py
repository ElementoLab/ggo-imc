import pandas as pd
import numpy as np
import os
import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt

import imc_analysis as imc

# load data
metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')
os.makedirs('figures/figure5', exist_ok=True)

# Read in data
adata = sc.read(
    metadata['patient_celltype_broad_clustered'],
    backup_url = metadata['patient_group_url']
)

# clustering
sc.pp.scale(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

group = 'Group'
features = metadata['CONDITIONS']

# Measure Scores for processes based on cellular abundance
processes = {
    'Epithelial score': ['Epi.','Epi. Prol.','Epi.-like','Tumor-like'],
    'Mesenchymal score': ['Mesen.-like'],
    'Fibrosis score': ['Fib.'],
    'Angiogenesis score': ['Endo.'],
    'B&T score': ['B', 'CD4 T','CD8 T',],
    'Pan-Immune score': ['CD4 T', 'CD8 T', 'Mac.', 'NK', 'T reg', 'Mast', 'Mono.', 'Neut.', 'PMN-MDSC']
}

adata.obs[group] = adata.obs[group].str.replace('roup','')
for s in processes:
    sc.tl.score_genes(adata, gene_list = processes[s],score_name=s, ctrl_size = 20, n_bins = 10)

adata.obs['EMT score'] = adata.obs['Mesenchymal score'] - adata.obs['Epithelial score']
scores = ['EMT score', 'Epithelial score','Mesenchymal score','Fibrosis score', 'Angiogenesis score', 'B&T score', 'Pan-Immune score']
# scores = ['EMT score', 'Fibrosis score', 'Angiogenesis score', 'B&T score', 'Pan-Immune score']
adata.uns['Group_colors'] = ['#99B898', '#FECEAB', '#E84A5F', '#2A363B']

from scipy.stats import zscore

fig, axes = plt.subplots(1,len(scores), figsize = (12,2.5), dpi = 300)
for i, ax in enumerate(axes.flatten()):
    adata.obs[scores[i]] = zscore(adata.obs[scores[i]])
    sc.pl.violin(adata, scores[i], groupby=group, stripplot=False, inner = 'box', show = False, ax = ax)
    #ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    ax.set_xlabel('')
sns.despine()
plt.tight_layout()
plt.savefig('figures/figure5/group_process_scores.pdf')
plt.close()

# Conditional Probability
def cond_prob(
    data,y,x,
):
    return (data
    .groupby(x)[y]
    .value_counts(normalize=True)
    .mul(100)
    .rename('percent')
    .reset_index()
    .pivot_table(index = x, columns = y, values='percent')
    .reset_index())


# P(histology | radiology)
p_hist_radio = cond_prob(adata.obs, x = 'radio', y = 'pathology')

# P(histology | radiology, group)
p_hist_radio_group = cond_prob(adata.obs, x = ['radio', group],y = 'pathology')
groups = p_hist_radio_group.groupby(group)

fig, axes = plt.subplots(1,len(groups) + 1, sharey=True, figsize=(7,2.5), dpi = 300)

p_hist_radio.plot(x = 'radio',kind='bar',stacked = True, ax=axes[0], color = metadata['pathology_color'][1:])
axes[0].set_xlabel('All Groups')
axes[0].set_ylabel(r'P($Histology|Radiology$)')
axes[0].legend().remove()

for ax, (group, data) in zip(axes[1:], groups):
    data.set_index("radio").plot(kind="bar", stacked = True, ax=ax, color = metadata['pathology_color'][1:])
    ax.tick_params(axis='both', which='both', length=0)
    ax.set_xlabel(f'{group}')
    ax.set_ylabel(r'P($Histology|Radiology$)')
    ax.legend().remove()
sns.despine()
plt.tight_layout()
plt.savefig('figures/figure5/group_conditional_probability_hist_radio.pdf')
plt.close()
# plt.show()


# Information Gained
# I(histology | radiology, group) = H(histology | radiology) - H(histology | radiology, group)
print(p_hist_radio)
print(p_hist_radio_group)
