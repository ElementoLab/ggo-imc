import imc_analysis as imc
import scanpy as sc
import anndata

import matplotlib

matplotlib.use('Agg')
metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])


'''
PANEL_H
Lymphocyte Analysis
'''
h_lymphocyte_markers = ['CD20', 'CD3', 'CD4', 'CD8a', 'CD103', 'HLAABC', 'VISTA', 'MMP7', 'IFNg','IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R']
h_lymphocyte_idx = adata_dict['PANEL_H'].obs['celltype'].isin(['T reg', 'CD4 T', 'CD8 T', 'B'])
h_lymphocytes = adata_dict['PANEL_H'][h_lymphocyte_idx,h_lymphocyte_markers].copy()

# sc.pp.scale(h_lymphocytes)
# sc.tl.pca(h_lymphocytes)
# sc.tl.umap(h_lymphocytes)
# h_lymphocytes.write(metadata['PANEL_H']['AnnData']['lymphocytes'])

h_lymphocytes = sc.read(
    metadata['PANEL_H']['AnnData']['lymphocytes'],
    backup_url = metadata['PANEL_H']['AnnData']['lymphocytes_url']
)

# % of cells with proinflammatory signals
cytokine = ['IFNg','IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R']
cytokine_expression = h_lymphocytes[:,cytokine].X.mean(axis=1)

# barplot cytokine positive cell density
h_lymphocytes.obs['mean cytokine expression'] = cytokine_expression.tolist()
h_lymphocytes.obs['cytokine+'] = (h_lymphocytes.obs['mean cytokine expression'] > 0).ravel().tolist()

cytokine_density = imc.tl.celltype_density(h_lymphocytes, celltype = 'cytokine+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(cytokine_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        cytokine_density,
        kind = 'bar',
        save_dir=f'figures/panel_h/cytokine+_density/',
        pval_form=pval
    )


'''
PANEL_G
Lymphocyte Analysis
'''
g_lymphocyte_markers = ['CD20', 'CD3', 'CD4', 'FoxP3', 'CD8a', 'CD45', 'CD45RA', 'CD45RO', 'GranzymeB', 'Tbet', 'ICOS', 'CD25', 'CD27', 'CD28', '41BB', 'CCR4', 'CCR7', 'PD1', 'TIM3', 'CTLA4', 'LAG3']
g_lymphocyte_idx = adata_dict['PANEL_G'].obs['celltype'].isin(['T reg', 'CD4 T', 'CD8 T', 'B'])
g_lymphocytes = adata_dict['PANEL_G'][g_lymphocyte_idx,g_lymphocyte_markers].copy()

# sc.pp.scale(g_lymphocytes)
# sc.tl.pca(g_lymphocytes)
# sc.tl.umap(g_lymphocytes)
# g_lymphocytes.write(metadata['PANEL_G']['AnnData']['lymphocytes'])

g_lymphocytes = sc.read(
    metadata['PANEL_G']['AnnData']['lymphocytes'],
    backup_url = metadata['PANEL_G']['AnnData']['lymphocytes_url']
)

pro_inflammatory_markers = ['CCR4', 'CCR7', 'GranzymeB', 'CD27', 'CD28', 'ICOS', 'Tbet', '41BB']

# inflammatory cell density
pro_inflammatory_marker_expression = g_lymphocytes[:,pro_inflammatory_markers].X.mean(axis=1)
g_lymphocytes.obs['mean inflammatory marker expression'] = pro_inflammatory_marker_expression.tolist()
sc.pl.umap(g_lymphocytes, color = 'mean inflammatory marker expression', save = 'panel_h_lymphocyte_umap_pro_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_lymphocytes.obs['pro_inflammatory_marker+'] = (g_lymphocytes.obs['mean inflammatory marker expression'] > 0).ravel().tolist()

pro_inflammatory_marker_density = imc.tl.celltype_density(g_lymphocytes, celltype = 'pro_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
pro_inflammatory_marker_density = pro_inflammatory_marker_density[:,True]
pro_inflammatory_marker_density.var.index = ['pro-inflammatory']
imc.tl.grouped_mwu_test(pro_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        pro_inflammatory_marker_density,
        kind = 'bar',
        y_max = 500,
        save_dir=f'figures/panel_g/pro_inflammatory_cell_density/',
        pval_form=pval
    )


anti_inflammatory_markers = ['FoxP3', 'CTLA4', 'PD1', 'LAG3', 'TIM3']
# regulatory cell density
anti_inflammatory_marker_expression = g_lymphocytes[:,anti_inflammatory_markers].X.mean(axis=1)
g_lymphocytes.obs['mean regulatory marker expression'] = anti_inflammatory_marker_expression.tolist()
sc.pl.umap(g_lymphocytes, color = 'mean regulatory marker expression', save = 'panel_h_lymphocyte_umap_anti_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_lymphocytes.obs['anti_inflammatory_marker+'] = (g_lymphocytes.obs['mean regulatory marker expression'] > 0).ravel().tolist()

anti_inflammatory_marker_density = imc.tl.celltype_density(g_lymphocytes, celltype = 'anti_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
anti_inflammatory_marker_density = anti_inflammatory_marker_density[:,True]
anti_inflammatory_marker_density.var.index = ['anti-inflammatory']

imc.tl.grouped_mwu_test(anti_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        anti_inflammatory_marker_density,
        kind = 'bar',
        y_max = 500,
        save_dir=f'figures/panel_g/anti_inflammatory_cell_density/',
        pval_form=pval
    )




import scipy as sp
import matplotlib.pyplot as plt
import seaborn as sns

# FoxP3+ cells with CD4 and CD8 T cells
lymphocyte_density = imc.tl.celltype_density(g_lymphocytes, celltype = 'celltype', condition_keys = ['pathology', 'radio'])
lymphocyte_density_df = lymphocyte_density.to_df()
fig, ax = plt.subplots(1,1,figsize = (3,3), dpi = 300)
sns.regplot(x="CD8 T", y="T reg", data=lymphocyte_density_df, ax = ax, color='#2b7bba', truncate = False, line_kws=dict(color="#2b7bba", linewidth = 0.5, alpha = 0.5), scatter_kws={'s':0.5})
sns.regplot(x="CD4 T", y="T reg", data=lymphocyte_density_df, ax = ax, color='#3030fd', truncate = False, line_kws=dict(color="#3030fd", linewidth = 0.5, alpha = 0.5), scatter_kws={'s':0.5})

ax.set_xlabel('')
ax.set_ylabel('')
r, p_ = sp.stats.spearmanr(lymphocyte_density_df['CD8 T'], lymphocyte_density_df['T reg'])
ax.text(.05, .9, 'r={:.2f}, p={:.2g}'.format(r, p_), color = '#2b7bba',
        transform=ax.transAxes)
r, p_ = sp.stats.spearmanr(lymphocyte_density_df['CD4 T'], lymphocyte_density_df['T reg'])
ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p_), color = '#3030fd',
        transform=ax.transAxes)
ax.set_xlim(0, 1000)
ax.set_ylim(0, 80)
sns.despine()
plt.savefig('figures/T reg relative abundance.pdf', bbox_inches = 'tight')
plt.close()