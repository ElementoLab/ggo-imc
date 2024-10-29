
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

import logging
logger = logging.getLogger(__name__)


def plot_mwu_bar(
    adata: anndata.AnnData,
    line_width:float = 0.5,
    save_dir: Path = 'figures/',
    palette: list = [],
    pval_form: str = 'star',
    show: bool = False,
    nrow = 3,
    ncol = 6,
    figsize = (12,6)
):

    if 'mwu' not in adata.uns:
        logger.error("'mwu' not found in adata.uns layer. Run imc_analysis.tl.grouped_mwu_test()")
        return
    
    density = adata.to_df()
    mwu_ = adata.uns['mwu']
    # for each conditions
    for cond in tqdm(mwu_['condition'].unique()):
        
        logger.info(f'Producing figure for condition: {cond}')
        
        fig, axes = plt.subplots(nrow,ncol, dpi = 300, figsize = figsize)
        
        # for each celltype
        for i, ax in enumerate(axes.flatten()):
        
            if i >= density.shape[1]:
                ax.axis('off')
            else:
                
                ct = density.columns[i]
                # create swarmboxen plot
                sns.barplot(y = density[ct], x = adata.obs[cond], hue = adata.obs[cond], ax = ax, palette = palette, dodge = False, width = 0.7)
                # sns.swarmplot(y = density[ct], x = adata.obs[cond], color = 'black', ax = ax, s = 2, alpha=0.5)
                # sns.lineplot(y = density[ct], x = adata.obs[cond], ax = ax, color = 'black', estimator='mean', markers=True, dashes=False, linewidth = 2, marker='o', markeredgecolor='black')
                ax.set_title(ct)
                ax.set_ylabel('')
                ax.set_xlabel('')
                ax.set_ylim(0,500)

                n98 = np.percentile(density[ct], 98)
                n98 = 500
                ax.set_ylim(0 - n98*0.05, n98 * 1.2)

                y1 = n98 * 1.05
                r = y1 * 0.03
                l = adata.obs[cond].cat.categories.tolist()

                sig_mwu_ = mwu_[(mwu_['celltype'] == ct) & (mwu_['condition'] == cond) & (mwu_['adj. p-val'] < 0.05)]

                sig_n = 0
                for i, idx in enumerate(sig_mwu_.index):
                    pair = sig_mwu_.loc[idx, 'pair']
                    
                    if pval_form == 'star':
                        pval = sig_mwu_.loc[idx, 'significance']
                    else:
                        pval = sig_mwu_.loc[idx, 'adj. p-val sci']
                    
                    if len(sig_mwu_) == 1:

                        ax.plot([pair[0], pair[1]], [y1 + r*i, y1 + r*i], lw=line_width, c='black')
                        ax.text(s = pval, x = 0.5, y = y1, fontsize = 8, va = 'bottom', ha = 'center')

                    else:
                        p = sig_mwu_.iloc[i]['adj. p-val']
                        
                        ax.plot([pair[0], pair[1]], [y1 + r*sig_n, y1 + r*sig_n], lw=line_width, c='black')
                        ax.text(s = pval, x = pair[1], y = y1 + r*sig_n, fontsize = 8, va = 'top', ha = 'left')
                        sig_n += 1
                ax.legend().remove()
                sns.despine()
        plt.tight_layout()

        dir_path = f'{save_dir}/{cond}_{pval_form}.pdf'
        # check if directory exists
        if not os.path.exists(save_dir):
            # create directory if it doesn't exist
            os.makedirs(save_dir)
            print(f"Directory '{save_dir}' created.")
        plt.savefig(dir_path, bbox_inches = 'tight')
        
        if show:
            plt.show()
        
        plt.close()

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

    del adata_dict[p].obs['ROI_area']


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


for p in ['PANEL_G', 'PANEL_H']:
    sc.pl.umap(adata_dict[p], color = adata_dict[p].var.index, save = f'{p}_umap_var.png', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)


'''
PANEL_H
Lymphocyte Analysis
'''
['Mac. (CD68+, CD163-)', 'Mac. (CD163+)', 'Mast', 'NK','PMN-MDSC', 'Neut.', 'B', 'CD4 T', 'CD8 T']
['Mac. (CD163+)', 'Mac. (CD163+, CD206+)', 'NK', 'T reg', 'CD4 T', 'CD8 T', 'B']
['aSMA','PanCK', 'Vimentin', 'RAGE', 'SFTPC', 'CD31', 'FOLR1', 'NKp44', 'CD14', 'CD33', 'CD163', 'CD16', 'CD66b']

h_lymphocyte_markers = ['CD20', 'CD3', 'CD4', 'CD8a', 'CD103', 'HLAABC', 'VISTA', 'MMP7', 'IFNg','IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R']

h_lymphocyte_idx = adata_dict['PANEL_H'].obs['celltype'].isin(['T reg', 'CD4 T', 'CD8 T', 'B'])
# h_lymphocyte_var = adata_dict['PANEL_H'].var.index.isin(h_lymphocyte_markers)
h_lymphocytes = adata_dict['PANEL_H'][h_lymphocyte_idx,h_lymphocyte_markers].copy()

# h_lymphocytes.write('results/h_lymphocytes.h5ad')

sc.pp.scale(h_lymphocytes)
sc.tl.pca(h_lymphocytes)
sc.pp.neighbors(h_lymphocytes)
sc.tl.umap(h_lymphocytes)

h_lymphocytes.write('results/h_lymphocytes.h5ad')

sc.pl.umap(h_lymphocytes, color = 'celltype_broad', save = 'panel_h_lymphocyte_umap', frameon = False, title = '')
sc.pl.umap(h_lymphocytes, color = h_lymphocytes.var.index, save = 'panel_h_lymphocyte_umap_var', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)

# % of cells with proinflammatory signals
# inflammatory_score = np.exp()/sum(np.exp())
cytokine = ['IFNg','IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R']
cytokine_expression = h_lymphocytes[:,cytokine].X.mean(axis=1)

# barplot cytokine positive cell density
h_lymphocytes.obs['mean cytokine expression'] = cytokine_expression.tolist()
sc.pl.umap(h_lymphocytes, color = 'mean cytokine expression', save = 'panel_h_lymphocyte_umap_cytokine', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
h_lymphocytes.obs['cytokine+'] = (h_lymphocytes.obs['mean cytokine expression'] > 0).ravel().tolist()

cytokine_density = imc.tl.celltype_density(h_lymphocytes, celltype = 'cytokine+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(cytokine_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        cytokine_density,
        save_dir=f'figures/panel_h/cytokine+_density/',
        palette=metadata['pathology_color'],
        pval_form=pval
    )

# upsetplot cytokine positive cell
h_lymphocytes_cytokine = h_lymphocytes[:,cytokine].copy()
for cyto in cytokine:
    h_lymphocytes_cytokine.obs[f'{cyto}+'] = (h_lymphocytes[:,cyto].X > 0).ravel().tolist()
cytokine_pos_cellcount = h_lymphocytes_cytokine.obs.groupby(['IFNg+', 'IL1alpha+', 'IL1beta+', 'IL1R1+', 'IL12p40+', 'IL17A+', 'IL23p19+', 'IL23R+']).count()['sample']

import upsetplot
upsetplot.plot(cytokine_pos_cellcount, sort_by='cardinality')
plt.savefig('figures/cytokine_upsetplot.pdf')
plt.close()


# for p in ['PANEL_G', 'PANEL_H']:
#     adata = adata_dict[p].copy()
#     adata


'''
PANEL_G
Lymphocyte Analysis
'''



for idx in g_t:
    if os.path.exists(f'results/g_{idx}.h5ad'):
        g_t[idx] = sc.read(f'results/g_{idx}.h5ad')
    else:
        sc.pp.scale(g_t[idx])
        sc.tl.pca(g_t[idx])
        sc.pp.neighbors(g_t[idx])
        sc.tl.umap(g_t[idx])

        g_t[idx].write(f'results/g_{idx}.h5ad')

for idx in g_t:
    sc.pl.umap(g_t[idx], color = 'pathology', save = f'PANEL_G/{idx}_umap_pathology', frameon = False, title = '')
    sc.pl.umap(g_t[idx], color = 'celltype_broad', save = f'PANEL_G/{idx}_umap', frameon = False, title = '')
    sc.pl.umap(g_t[idx], color = g_t[idx].var.index, save = f'PANEL_G/{idx}_umap_var', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)


# inflammatory cell upset plot
g_lymphocytes_functional_marker = g_lymphocytes[:,pro_inflammatory_markers].copy()
for cyto in pro_inflammatory_markers:
    g_lymphocytes_functional_marker.obs[f'{cyto}+'] = (g_lymphocytes[:,cyto].X > 0).ravel().tolist()
functional_marker_pos_cellcount = g_lymphocytes_functional_marker.obs.groupby([f'{x}+' for x in pro_inflammatory_markers]).count()['sample']

upsetplot.plot(functional_marker_pos_cellcount, sort_by='cardinality')
plt.savefig('figures/inflammatory+_marker_upsetplot.pdf')
plt.close()

# regulatory cell upset plot
g_lymphocytes_functional_marker = g_lymphocytes[:,anti_inflammatory_markers].copy()
for cyto in anti_inflammatory_markers:
    g_lymphocytes_functional_marker.obs[f'{cyto}+'] = (g_lymphocytes[:,cyto].X > 0).ravel().tolist()
functional_marker_pos_cellcount = g_lymphocytes_functional_marker.obs.groupby([f'{x}+' for x in anti_inflammatory_markers]).count()['sample']

upsetplot.plot(functional_marker_pos_cellcount, sort_by='cardinality')
plt.savefig('figures/inflammatory-_marker_upsetplot.pdf')
plt.close()


# inflammatory cell density
pro_inflammatory_marker_expression = g_lymphocytes[:,pro_inflammatory_markers].X.mean(axis=1)
g_lymphocytes.obs['mean inflammatory marker expression'] = pro_inflammatory_marker_expression.tolist()
sc.pl.umap(g_lymphocytes, color = 'mean inflammatory marker expression', save = 'panel_h_lymphocyte_umap_pro_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_lymphocytes.obs['pro_inflammatory_marker+'] = (g_lymphocytes.obs['mean inflammatory marker expression'] > 0).ravel().tolist()


pro_inflammatory_marker_density = imc.tl.celltype_density(g_lymphocytes, celltype = 'pro_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(pro_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        pro_inflammatory_marker_density,
        save_dir=f'figures/panel_g/pro_inflammatory_cell_density/',
        palette=metadata['pathology_color'],
        pval_form=pval
    )


# regulatory cell density
anti_inflammatory_marker_expression = g_lymphocytes[:,anti_inflammatory_markers].X.mean(axis=1)
g_lymphocytes.obs['mean regulatory marker expression'] = anti_inflammatory_marker_expression.tolist()
sc.pl.umap(g_lymphocytes, color = 'mean regulatory marker expression', save = 'panel_h_lymphocyte_umap_anti_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_lymphocytes.obs['anti_inflammatory_marker+'] = (g_lymphocytes.obs['mean regulatory marker expression'] > 0).ravel().tolist()

anti_inflammatory_marker_density = imc.tl.celltype_density(g_lymphocytes, celltype = 'anti_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(anti_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        anti_inflammatory_marker_density,
        save_dir=f'figures/panel_g/anti_inflammatory_cell_density/',
        palette=metadata['pathology_color'],
        pval_form=pval
    )




import scipy as sp

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