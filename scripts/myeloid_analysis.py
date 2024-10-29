
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

'''
PANEL_H
Myeloid Analysis
'''

['Mac. (CD68+, CD163-)', 'Mac. (CD163+)', 'Mast', 'NK','PMN-MDSC', 'Neut.', 'B', 'CD4 T', 'CD8 T']
['Mac. (CD163+)', 'Mac. (CD163+, CD206+)', 'NK', 'T reg', 'CD4 T', 'CD8 T', 'B']

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


if os.path.exists('results/h_myeloid.h5ad'):
    h_myeloid = sc.read('results/h_myeloid.h5ad')
else:
    sc.pp.scale(h_myeloid)
    sc.tl.pca(h_myeloid)
    sc.pp.neighbors(h_myeloid)
    sc.tl.umap(h_myeloid)

    h_myeloid.write('results/h_myeloid.h5ad')

sc.pl.umap(h_myeloid, color = 'celltype_broad', save = 'panel_h_myeloid_umap', frameon = False, title='')
sc.pl.umap(h_myeloid, color = h_myeloid.var.index, save = 'panel_h_myeloid_umap_var', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)



cytokine = ['IFNg','IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R']
cytokine_expression = h_myeloid[:,cytokine].X.mean(axis=1)

# barplot cytokine positive cell density
h_myeloid.obs['mean cytokine expression'] = cytokine_expression.tolist()
sc.pl.umap(h_myeloid, color = 'mean cytokine expression', save = 'panel_h_myeloid_umap_cytokine', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
h_myeloid.obs['cytokine+'] = (h_myeloid.obs['mean cytokine expression'] > 0).ravel().tolist()

cytokine_density = imc.tl.celltype_density(h_myeloid, celltype = 'cytokine+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(cytokine_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        cytokine_density,
        save_dir=f'figures/panel_h/cytokine+_myeloid_density/',
        palette=metadata['pathology_color'],
        pval_form=pval
    )

# upsetplot cytokine positive cell
h_myeloid_cytokine = h_myeloid[:,cytokine].copy()
for cyto in cytokine:
    h_myeloid_cytokine.obs[f'{cyto}+'] = (h_myeloid[:,cyto].X > 0).ravel().tolist()
cytokine_pos_cellcount = h_myeloid_cytokine.obs.groupby(['IFNg+', 'IL1alpha+', 'IL1beta+', 'IL1R1+', 'IL12p40+', 'IL17A+', 'IL23p19+', 'IL23R+']).count()['sample']

import upsetplot
upsetplot.plot(cytokine_pos_cellcount, sort_by='cardinality')
plt.savefig('figures/myeloid_cytokine_upsetplot.pdf')
plt.close()




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

if os.path.exists('results/g_myeloid.h5ad'):
    g_myeloid = sc.read('results/g_myeloid.h5ad')
else:
    sc.pp.scale(g_myeloid)
    sc.tl.pca(g_myeloid)
    sc.pp.neighbors(g_myeloid)
    sc.tl.umap(g_myeloid)

    g_myeloid.write('results/g_myeloid.h5ad')

sc.pl.umap(g_myeloid, color = 'pathology', save = 'panel_g_myeloid_umap_pathology', frameon = False, title = '')
sc.pl.umap(g_myeloid, color = 'celltype_broad', save = 'panel_g_myeloid_umap', frameon = False, title = '')
sc.pl.umap(g_myeloid, color = g_myeloid.var.index, save = 'panel_g_myeloid_umap_var', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)



# inflammatory cell upset plot
g_myeloid_functional_marker = g_myeloid[:,pro_inflammatory_markers].copy()
for cyto in pro_inflammatory_markers:
    g_myeloid_functional_marker.obs[f'{cyto}+'] = (g_myeloid[:,cyto].X > 0).ravel().tolist()
functional_marker_pos_cellcount = g_myeloid_functional_marker.obs.groupby([f'{x}+' for x in pro_inflammatory_markers]).count()['sample']

upsetplot.plot(functional_marker_pos_cellcount, sort_by='cardinality')
plt.savefig('figures/myeloid_inflammatory+_marker_upsetplot.pdf')
plt.close()

# regulatory cell upset plot
g_myeloid_functional_marker = g_myeloid[:,anti_inflammatory_markers].copy()
for cyto in anti_inflammatory_markers:
    g_myeloid_functional_marker.obs[f'{cyto}+'] = (g_myeloid[:,cyto].X > 0).ravel().tolist()
functional_marker_pos_cellcount = g_myeloid_functional_marker.obs.groupby([f'{x}+' for x in anti_inflammatory_markers]).count()['sample']

upsetplot.plot(functional_marker_pos_cellcount, sort_by='cardinality')
plt.savefig('figures/myeloid_inflammatory-_marker_upsetplot.pdf')
plt.close()


# inflammatory cell density
pro_inflammatory_marker_expression = g_myeloid[:,pro_inflammatory_markers].X.mean(axis=1)
g_myeloid.obs['mean inflammatory marker expression'] = pro_inflammatory_marker_expression.tolist()
sc.pl.umap(g_myeloid, color = 'mean inflammatory marker expression', save = 'panel_h_lymphocyte_umap_pro_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_myeloid.obs['pro_inflammatory_marker+'] = (g_myeloid.obs['mean inflammatory marker expression'] > 0).ravel().tolist()

pro_inflammatory_marker_density = imc.tl.celltype_density(g_myeloid, celltype = 'pro_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(pro_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        pro_inflammatory_marker_density,
        save_dir=f'figures/panel_g/pro_inflammatory_myeloid_cell_density/',
        palette=metadata['pathology_color'],
        pval_form=pval
    )


# regulatory cell density
anti_inflammatory_marker_expression = g_myeloid[:,anti_inflammatory_markers].X.mean(axis=1)
g_myeloid.obs['mean regulatory marker expression'] = anti_inflammatory_marker_expression.tolist()
sc.pl.umap(g_myeloid, color = 'mean regulatory marker expression', save = 'panel_h_lymphocyte_umap_anti_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_myeloid.obs['anti_inflammatory_marker+'] = (g_myeloid.obs['mean regulatory marker expression'] > 0).ravel().tolist()

anti_inflammatory_marker_density = imc.tl.celltype_density(g_myeloid, celltype = 'anti_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(anti_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        anti_inflammatory_marker_density,
        save_dir=f'figures/panel_g/anti_inflammatory_myeloid_cell_density/',
        palette=metadata['pathology_color'],
        pval_form=pval
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
    plot_mwu_bar(
        macrophage_polarization,
        save_dir=f'figures/panel_h/macrophage_polarization_density/',
        palette=metadata['pathology_color'],
        pval_form=pval
    )