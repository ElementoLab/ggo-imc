import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt

import imc_analysis as imc
import os
from utils import CYTOKINE as cytokine

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(
        metadata[p]['AnnData']['phenotyped_umap_name'],
        backup_url = metadata[p]['AnnData']['backup_url']
    )


'''
PANEL_H
Myeloid Analysis
'''
h_myeloid_markers = [
    'CD56', 'NKp44', 'IFNg', 'CD15', 'CD66b', 'VISTA', 'CD117',
    'CD14', 'CD16', 'CD68', 'CD163', 'HLADR', 'MMP7', 'CD33',
    'IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R'
]


h_myeloid_idx = adata_dict['PANEL_H'].obs['celltype'].isin(['Mac. (CD68+, CD163-)', 'Mac. (CD163+)', 'Mast', 'NK','PMN-MDSC', 'Neut.'])
# myeloid_var = adata_dict['PANEL_H'].var.index.isin(['aSMA','PanCK', 'Vimentin', 'RAGE', 'SFTPC', 'CD31', 'FOLR1', 'CD3', 'CD4', 'CD8a', 'CD20'])
h_myeloid = adata_dict['PANEL_H'][h_myeloid_idx,h_myeloid_markers].copy()


h_myeloid.obs['celltype_broad'] = h_myeloid.obs['celltype_broad'].astype(str).replace({'Neutrophil':'Monocyte'})
h_myeloid.obs['celltype_broad'] = pd.Categorical(h_myeloid.obs['celltype_broad'], categories = ['Monocyte', 'Macrophage', 'Mast', 'NK', 'PMN-MDSC'])
# myeloid.write('results/myeloid.h5ad')


if os.path.exists(metadata['PANEL_H']['AnnData']['myeloids']):
    h_myeloid = sc.read(metadata['PANEL_H']['AnnData']['myeloids'],
        backup_url = metadata['PANEL_H']['AnnData']['myeloids_url'])
else:
    sc.pp.scale(h_myeloid)
    sc.tl.pca(h_myeloid)
    sc.pp.neighbors(h_myeloid)
    sc.tl.umap(h_myeloid)

    h_myeloid.write(metadata['PANEL_H']['AnnData']['myeloids'])

cytokine_expression = h_myeloid[:,cytokine].X.mean(axis=1)

'''
PANEL_G
Myeloid Analysis
'''
g_myeloid_markers = ['HLADR', 'CD14', 'CD16', 'CD68', 'CD163', 'CD206', 'CD11c', 'CD11b', 'CD56', 'CD57',
'GranzymeB', 'Tbet', 'ICOS', 'CD25', 'CD27', 'CD28','41BB', 'CCR4', 'CCR7', 'PD1', 'PDL1', 'TIM3', 'CTLA4', 'LAG3']

functional_markers = ['GranzymeB', 'Tbet', 'ICOS', 'CD25', 'CD27', 'CD28','41BB', 'CCR4', 'CCR7', 'PD1', 'PDL1', 'TIM3', 'CTLA4', 'LAG3']
pro_inflammatory_markers = ['GranzymeB', 'Tbet', 'ICOS', 'CD27', 'CD28','41BB', 'CCR4', 'CCR7']
anti_inflammatory_markers = ['PD1', 'PDL1', 'TIM3', 'CTLA4', 'CD25', 'LAG3']

g_myeloid_idx = adata_dict['PANEL_G'].obs['celltype'].isin(['Mac. (CD163+)', 'Mac. (CD163+, CD206+)', 'NK'])
# g_myeloid_var = adata_dict['PANEL_H'].var.index.isin(g_myeloid_markers)
g_myeloid = adata_dict['PANEL_G'][g_myeloid_idx,g_myeloid_markers].copy()

# g_myeloid.write('results/g_myeloid.h5ad')

if os.path.exists(metadata['PANEL_G']['AnnData']['myeloids']):
    g_myeloid = sc.read(metadata['PANEL_G']['AnnData']['myeloids'],
        backup_url = metadata['PANEL_G']['AnnData']['myeloids_url'])
else:
    sc.pp.scale(g_myeloid)
    sc.tl.pca(g_myeloid)
    sc.pp.neighbors(g_myeloid)
    sc.tl.umap(g_myeloid)

    g_myeloid.write(metadata['PANEL_G']['AnnData']['myeloids'])


# M2/M1 (non-M2) ratio
macrophages = h_myeloid[h_myeloid.obs['celltype_broad'] == 'Macrophage'].copy()
macrophages.obs['macrophage_subtype'] = (macrophages[:,'CD163'].X > 0).flatten().tolist()
macrophages.obs['macrophage_subtype'] = macrophages.obs['macrophage_subtype'].replace({True: 'M2', False: 'M1'})
# tmp = macrophages[:,['CD14', 'CD16', 'CD68', 'CD163']]

polar = macrophages[:,['CD68', 'CD163']].to_df()
polar['subtype'] = macrophages.obs['macrophage_subtype']
sns.kdeplot(
    data=polar, x="CD68", y="CD163", fill = True, hue= 'subtype',
)
sns.despine()
plt.savefig('figures/figure2/myeloid_polarization.pdf')
plt.close()

counts = macrophages.obs.groupby(['pathology','macrophage_subtype']).count()[['sample']].pivot_table(index = 'pathology', columns = 'macrophage_subtype')
counts / np.array(counts.sum(axis=1))[:,None]


macrophage_polarization = imc.tl.celltype_density(macrophages, celltype = 'macrophage_subtype', condition_keys = ['pathology'])
imc.tl.grouped_mwu_test(macrophage_polarization, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        macrophage_polarization,
        save_dir=f'figures/figure2/macrophage_polarization_density/',
        kind = 'bar',
        pval_form=pval,
        y_max = 500,
    )





# g_myeloid.uns['pathology_colors'] = adata_dict['PANEL_G'].uns['pathology_colors']
# h_myeloid.uns['pathology_colors'] = adata_dict['PANEL_G'].uns['pathology_colors']
# g_myeloid.uns['radio_colors'] = adata_dict['PANEL_G'].uns['radio_colors']
# h_myeloid.uns['radio_colors'] = adata_dict['PANEL_G'].uns['radio_colors']
# g_myeloid.write(metadata['PANEL_G']['AnnData']['myeloids'])
# h_myeloid.write(metadata['PANEL_H']['AnnData']['myeloids'])
