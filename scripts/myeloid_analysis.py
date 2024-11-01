import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt

import imc_analysis as imc
import os

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])

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

# sc.pl.umap(h_myeloid, color = 'celltype_broad', save = 'panel_h_myeloid_umap', frameon = False, title='')
# sc.pl.umap(h_myeloid, color = h_myeloid.var.index, save = 'panel_h_myeloid_umap_var', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)


cytokine = ['IFNg','IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R']
cytokine_expression = h_myeloid[:,cytokine].X.mean(axis=1)

# barplot cytokine positive cell density
h_myeloid.obs['mean cytokine expression'] = cytokine_expression.tolist()
# sc.pl.umap(h_myeloid, color = 'mean cytokine expression', save = 'panel_h_myeloid_umap_cytokine', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
h_myeloid.obs['cytokine+'] = (h_myeloid.obs['mean cytokine expression'] > 0).ravel().tolist()

cytokine_density = imc.tl.celltype_density(h_myeloid, celltype = 'cytokine+', condition_keys = ['pathology', 'radio'])
cytokine_density = cytokine_density[:,True]
cytokine_density.var.index = ['Cytokine+']
imc.tl.grouped_mwu_test(cytokine_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        cytokine_density,
        save_dir=f'figures/panel_h/cytokine+_myeloid_density/',
        kind = 'bar',
        pval_form=pval,
        y_max = 500,
    )

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

# sc.pl.umap(g_myeloid, color = 'pathology', save = 'panel_g_myeloid_umap_pathology', frameon = False, title = '')
# sc.pl.umap(g_myeloid, color = 'celltype_broad', save = 'panel_g_myeloid_umap', frameon = False, title = '')
# sc.pl.umap(g_myeloid, color = g_myeloid.var.index, save = 'panel_g_myeloid_umap_var', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)


# inflammatory cell density
pro_inflammatory_marker_expression = g_myeloid[:,pro_inflammatory_markers].X.mean(axis=1)
g_myeloid.obs['mean inflammatory marker expression'] = pro_inflammatory_marker_expression.tolist()
# sc.pl.umap(g_myeloid, color = 'mean inflammatory marker expression', save = 'panel_h_lymphocyte_umap_pro_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_myeloid.obs['pro_inflammatory_marker+'] = (g_myeloid.obs['mean inflammatory marker expression'] > 0).ravel().tolist()

pro_inflammatory_marker_density = imc.tl.celltype_density(g_myeloid, celltype = 'pro_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
pro_inflammatory_marker_density = pro_inflammatory_marker_density[:,True]
pro_inflammatory_marker_density.var.index = ['pro-inflammatory']
imc.tl.grouped_mwu_test(pro_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        pro_inflammatory_marker_density,
        save_dir=f'figures/panel_g/pro_inflammatory_myeloid_cell_density/',
        kind = 'bar',
        pval_form=pval,
        y_max = 500,
    )


# regulatory cell density
anti_inflammatory_marker_expression = g_myeloid[:,anti_inflammatory_markers].X.mean(axis=1)
g_myeloid.obs['mean regulatory marker expression'] = anti_inflammatory_marker_expression.tolist()
# sc.pl.umap(g_myeloid, color = 'mean regulatory marker expression', save = 'panel_h_lymphocyte_umap_anti_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_myeloid.obs['anti_inflammatory_marker+'] = (g_myeloid.obs['mean regulatory marker expression'] > 0).ravel().tolist()

anti_inflammatory_marker_density = imc.tl.celltype_density(g_myeloid, celltype = 'anti_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
anti_inflammatory_marker_density = anti_inflammatory_marker_density[:,True]
anti_inflammatory_marker_density.var.index = ['anti-inflammatory']

imc.tl.grouped_mwu_test(anti_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        anti_inflammatory_marker_density,
        save_dir=f'figures/panel_g/anti_inflammatory_myeloid_cell_density/',
        kind = 'bar',
        pval_form=pval,
        y_max = 500,
    )


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
plt.savefig('figures/myeloid_polarization.pdf')
plt.close()

counts = macrophages.obs.groupby(['pathology','macrophage_subtype']).count()[['sample']].pivot_table(index = 'pathology', columns = 'macrophage_subtype')
counts / np.array(counts.sum(axis=1))[:,None]


macrophage_polarization = imc.tl.celltype_density(macrophages, celltype = 'macrophage_subtype', condition_keys = ['pathology'])
imc.tl.grouped_mwu_test(macrophage_polarization, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        macrophage_polarization,
        save_dir=f'figures/panel_h/macrophage_polarization_density/',
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
