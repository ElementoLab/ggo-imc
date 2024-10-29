
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
import upsetplot

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
    figsize = (12,6),
    y_lim = (0,500)
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
                ax.set_ylim(y_lim[0],y_lim[1])

                n98 = np.percentile(density[ct], 98)
                n98 = y_lim[1]
                ax.set_ylim(y_lim[0] - n98*0.05, n98 * 1.2)

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

if os.path.exists('results/h_lymphocytes.h5ad'):
    h_lymphocytes = sc.read('results/h_lymphocytes.h5ad')
else:
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
h_lymphocytes = h_lymphocytes[h_lymphocytes.obs['celltype_broad'].isin(['T reg', 'CD4 T', 'CD8 T'])].copy()

cytokine_density = imc.tl.celltype_density(h_lymphocytes, celltype = 'cytokine+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(cytokine_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        cytokine_density,
        save_dir=f'figures/panel_h/cytokine+_density/',
        palette=metadata['pathology_color'],
        pval_form=pval,
        y_lim = (0,300)
    )

# upsetplot cytokine positive cell
h_lymphocytes_cytokine = h_lymphocytes[:,cytokine].copy()
for cyto in cytokine:
    h_lymphocytes_cytokine.obs[f'{cyto}+'] = (h_lymphocytes[:,cyto].X > 0).ravel().tolist()
cytokine_pos_cellcount = h_lymphocytes_cytokine.obs.groupby(['IFNg+', 'IL1alpha+', 'IL1beta+', 'IL1R1+', 'IL12p40+', 'IL17A+', 'IL23p19+', 'IL23R+']).count()['sample']

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


g_lymphocyte_markers = ['CD20', 'CD3', 'CD4', 'FoxP3', 'CD8a', 'CD45', 'CD27', 'CD45RA', 'CD45RO', 'LAG3',
'GranzymeB', 'Tbet', 'ICOS', 'CD25', 'CD28','41BB', 'CCR4', 'CCR7', 'PD1', 'TIM3', 'CTLA4',]

functional_markers = ['GranzymeB', 'Tbet', 'ICOS', 'CD25', 'CD27', 'CD28','41BB', 'CCR4', 'CCR7', 'PD1', 'TIM3', 'CTLA4', 'LAG3']
pro_inflammatory_markers = ['GranzymeB', 'Tbet', 'ICOS', 'CD27', 'CD28','41BB', 'CCR4', 'CCR7']
anti_inflammatory_markers = ['PD1', 'TIM3', 'CTLA4', 'CD25', 'LAG3']

g_lymphocyte_idx = adata_dict['PANEL_G'].obs['celltype'].isin(['T reg', 'CD4 T', 'CD8 T', 'B'])
g_lymphocytes = adata_dict['PANEL_G'][g_lymphocyte_idx,g_lymphocyte_markers].copy()

if os.path.exists('results/g_lymphocytes.h5ad'):
    g_lymphocytes = sc.read('results/g_lymphocytes.h5ad')
else:
    sc.pp.scale(g_lymphocytes)
    sc.tl.pca(g_lymphocytes)
    sc.pp.neighbors(g_lymphocytes)
    sc.tl.umap(g_lymphocytes)

    g_lymphocytes.write('results/g_lymphocytes.h5ad')

sc.pl.umap(g_lymphocytes, color = 'pathology', save = 'panel_g_lymphocyte_umap_pathology', frameon = False, title = '')
sc.pl.umap(g_lymphocytes, color = 'celltype_broad', save = 'panel_g_lymphocyte_umap', frameon = False, title = '')
sc.pl.umap(g_lymphocytes, color = g_lymphocytes.var.index, save = 'panel_g_lymphocyte_umap_var', frameon = False, use_raw = False, vmin = -0.5, vmax = 1.5, colorbar_loc = None)

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
g_tcells = g_lymphocytes[g_lymphocytes.obs['celltype_broad'].isin(['CD8 T'])].copy()

pro_inflammatory_marker_density = imc.tl.celltype_density(g_tcells, celltype = 'pro_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(pro_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        pro_inflammatory_marker_density,
        save_dir=f'figures/panel_g/pro_inflammatory_cell_density/',
        palette=metadata['pathology_color'],
        pval_form=pval,
        y_lim = (0,200)
    )


# regulatory cell density
anti_inflammatory_marker_expression = g_lymphocytes[:,anti_inflammatory_markers].X.mean(axis=1)
g_lymphocytes.obs['mean regulatory marker expression'] = anti_inflammatory_marker_expression.tolist()
sc.pl.umap(g_lymphocytes, color = 'mean regulatory marker expression', save = 'panel_h_lymphocyte_umap_anti_inflammatory_marker', title = '', frameon = False, use_raw = False, vmin = -2, vmax = 2)
g_lymphocytes.obs['anti_inflammatory_marker+'] = (g_lymphocytes.obs['mean regulatory marker expression'] > 0).ravel().tolist()
g_tregs = g_lymphocytes[g_lymphocytes.obs['celltype_broad'].isin(['T reg'])].copy()

anti_inflammatory_marker_density = imc.tl.celltype_density(g_tregs, celltype = 'anti_inflammatory_marker+', condition_keys = ['pathology', 'radio'])
imc.tl.grouped_mwu_test(anti_inflammatory_marker_density, condition_keys = ['pathology'])
for pval in ['star', 'sci_not']:
    plot_mwu_bar(
        anti_inflammatory_marker_density,
        save_dir=f'figures/panel_g/anti_inflammatory_cell_density/',
        palette=metadata['pathology_color'],
        pval_form=pval,
        y_lim = (0,50)
    )


# combined heatmap
g_lymphocytes.obs['celltype_activation'] = g_lymphocytes.obs['celltype_broad'].astype(str)
g_lymphocytes.obs.loc[g_lymphocytes.obs['pro_inflammatory_marker+'] & (g_lymphocytes.obs['celltype_broad'] == 'CD8 T'), 'celltype_activation'] = 'Activated CD8 T'
g_lymphocytes.obs.loc[g_lymphocytes.obs['pro_inflammatory_marker+'] & (g_lymphocytes.obs['celltype_broad'] == 'CD4 T'), 'celltype_activation'] = 'Activated CD4 T'
g_lymphocytes.obs.loc[g_lymphocytes.obs['anti_inflammatory_marker+'] & (g_lymphocytes.obs['celltype_broad'] == 'T reg'), 'celltype_activation'] = 'Activated T reg'
g_lymphocytes.obs['celltype_activation'] = pd.Categorical(g_lymphocytes.obs['celltype_activation'], categories = ['B', 'CD8 T', 'Activated CD8 T', 'CD4 T', 'Activated CD4 T', 'Activated T reg'])

imc.pl.celltype_heatmap(
    g_lymphocytes,
    cluster_ids = ['celltype_activation'],
    out_dir='figures/g_lymphocytes/')




adata = adata_dict['PANEL_G']
adata.obs['celltype_activation'] = adata.obs['celltype_broad'].astype(str)
adata.obs['celltype_activation'] = adata.obs['celltype_activation'].replace({'Epithelial-like (Ki67+)': 'Epithelial-like'})
adata.obs.loc[g_lymphocytes.obs.index, 'celltype_activation'] = g_lymphocytes.obs['celltype_activation'].astype(str)
adata.obs['celltype_activation'] = pd.Categorical(
    adata.obs['celltype_activation'],
    categories = [
        'Epithelial-like',
        'Fibroblast', 'Mesenchymal-like', 'Endothelial',
        'B', 'CD8 T', 'Activated CD8 T', 'CD4 T', 'Activated CD4 T', 'Activated T reg',
        'NK', 'Macrophage',
    ]
)


# distance analysis

import pandas as pd
from sklearn.neighbors import BallTree
from joblib import Parallel, delayed
import itertools

def spatial_distance (adata,
                      x_coordinate='X_centroid',
                      y_coordinate='Y_centroid',
                      z_coordinate=None,
                      phenotype='phenotype',
                      subset=None,
                      imageid='imageid',
                      verbose=True,
                      label='spatial_distance'):


    def spatial_distance_internal (adata_subset,x_coordinate,y_coordinate,z_coordinate,
                                   phenotype,subset,imageid,label):

        if verbose:
            print("Processing Image: " + str(adata_subset.obs[imageid].unique()[0]))
        # Create a dataFrame with the necessary inforamtion
        data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})

        # Function to identify shortest distance for each phenotype of interest
        def distance (pheno):
            pheno_interest = data[data['phenotype'] == pheno]

            if len(pheno_interest) == 0:
                return []
            # print(pheno)
            # Build the ball-tree for search space
            tree = BallTree(pheno_interest[['x','y']], metric='euclidean') 
            # Calculate shortest distance (if statement to account for K)
            if len(pheno_interest) > 1:
                dist, ind = tree.query(data[['x','y']], k=2, return_distance= True)
                dist = pd.DataFrame(dist)
                dist.loc[dist[0] == 0, 0]  = dist[1]
                dist = dist[0].values
            else:
                dist, ind = tree.query(data[['x','y']], k=1, return_distance= True)
                dist = list(itertools.chain.from_iterable(dist))
            
            # print(dist)
            return dist

        # Run in parallel for all phenotypes
        phenotype_list = list(data['phenotype'].unique())
        # print(phenotype_list)
        # Apply function
        final_dist = Parallel(n_jobs=-1)(delayed(distance)(pheno=i) for i in phenotype_list)     
        final_dist = pd.DataFrame(final_dist, index = phenotype_list, columns = data.index).T

        return final_dist

    # subset a particular subset of cells if the user wants else break the adata into list of anndata objects
    adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]

    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_distance_internal = lambda x: spatial_distance_internal (adata_subset=x,
                                                                       x_coordinate=x_coordinate,y_coordinate=y_coordinate, z_coordinate=z_coordinate,
                                                                       phenotype=phenotype,subset=subset,imageid=imageid,label=label) 
    all_data = list(map(r_spatial_distance_internal, adata_list)) # Apply function 

    # Merge all the results into a single dataframe    
    result = []
    for i in range(len(all_data)):
        result.append(all_data[i])
    result = pd.concat(result, join='outer')  


    # Add to anndata
    adata.uns[label] = result

    # return
    return adata

dist = spatial_distance(adata, phenotype = 'celltype_activation', imageid='roi')
dist.write('results/panel_g.scimap.h5ad')

dist_h = spatial_distance(adata_dict['PANEL_H'], phenotype = 'celltype_broad', imageid='roi')
dist_h.write('results/panel_h.scimap.h5ad')



## Seperate script

% conda activate scimap
% ipython


# library
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def spatial_distance_plot(
    adata,
    spatial_distance='spatial_distance',
    phenotype='phenotype',
    imageid='imageid',
    log=False,
    method='heatmap',
    heatmap_summarize=True,
    heatmap_na_color='grey',
    heatmap_cmap='vlag_r',
    heatmap_row_cluster=False,
    heatmap_col_cluster=False,
    heatmap_standard_scale=0,
    distance_from=None,
    distance_to=None,
    x_axis=None,
    y_axis=None,
    facet_by=None,
    plot_type=None,
    return_data=False,
    subset_col=None,
    subset_value=None,
    fileName='spatial_distance.pdf',
    saveDir=None,
    **kwargs,
):
    # set color for heatmap
    # cmap_updated = matplotlib.cm.get_cmap(heatmap_cmap)
    cmap_updated = matplotlib.colormaps[heatmap_cmap]
    cmap_updated.set_bad(color=heatmap_na_color)

    # Copy the spatial_distance results from anndata object
    try:
        diatance_map = adata.uns[spatial_distance].copy()
    except KeyError:
        raise ValueError(
            'spatial_distance not found- Please run sm.tl.spatial_distance first'
        )

    # subset the data if user requests
    if subset_col is not None:
        if isinstance(subset_value, str):
            subset_value = [subset_value]
        # find the cell names to be subsetted out
        obs = adata.obs[[subset_col]]
        cells_to_subset = obs[obs[subset_col].isin(subset_value)].index

        # subset the diatance_map
        diatance_map = diatance_map.loc[
            diatance_map.index.intersection(cells_to_subset)
        ]
        # diatance_map = diatance_map.loc[cells_to_subset]

    # Convert distance to log scale if user requests
    if log is True:
        diatance_map = np.log1p(diatance_map)

    # Method
    if method == 'heatmap':
        if heatmap_summarize is True:
            # create the necessary data
            data = pd.DataFrame({'phenotype': adata.obs[phenotype]})
            data = pd.merge(
                data, diatance_map, how='outer', left_index=True, right_index=True
            )  # merge with the distance map
            k = data.groupby(
                ['phenotype'], observed=False
            ).mean()  # collapse the whole dataset into mean expression
            d = k[k.index]
        else:
            # create new naming scheme for the phenotypes
            non_summary = pd.DataFrame(
                {'imageid': adata.obs[imageid], 'phenotype': adata.obs[phenotype]}
            )
            non_summary['imageid'] = non_summary['imageid'].astype(
                str
            )  # convert the column to string
            non_summary['phenotype'] = non_summary['phenotype'].astype(
                str
            )  # convert the column to string
            non_summary['image_phenotype'] = non_summary['imageid'].str.cat(
                non_summary['phenotype'], sep="_"
            )
            # Merge distance map with phenotype
            data = pd.DataFrame(non_summary[['image_phenotype']])
            data = pd.merge(
                data, diatance_map, how='outer', left_index=True, right_index=True
            )
            k = data.groupby(['image_phenotype'], observed=False).mean()
            d = k.sort_index(axis=1)
        # Generate the heatmap
        mask = d.isnull()  # identify the NAN's for masking
        d = d.fillna(0)  # replace nan's with 0 so that clustering will work
        # Heatmap
        plot = sns.clustermap(
            d,
            cmap=heatmap_cmap,
            row_cluster=heatmap_row_cluster,
            col_cluster=heatmap_col_cluster,
            mask=mask,
            standard_scale=heatmap_standard_scale,
            **kwargs,
        )
    else:

        # condition-1
        if distance_from is None and distance_to is None:
            raise ValueError(
                'Please include distance_from and/or distance_to parameters to use this method'
            )

        # condition-2
        if distance_from is None and distance_to is not None:
            raise ValueError('Please `distance_from` parameters to use this method')

        # condition-3
        if distance_to is not None:
            # convert input to list if needed
            if isinstance(distance_to, str):
                distance_to = [distance_to]

        # Start
        pheno_df = pd.DataFrame(
            {'imageid': adata.obs[imageid], 'phenotype': adata.obs[phenotype]}
        )  # image id and phenotype
        data = pd.merge(
            pheno_df, diatance_map, how='outer', left_index=True, right_index=True
        )  # merge with the distance map
        data = data[data['phenotype'] == distance_from]  # subset the pheno of interest

        if distance_to is not None:
            data = data[
                distance_to
            ]  # drop columns that are not requested in distance_to
        else:
            data = data.drop(
                ['phenotype', 'imageid'], axis=1
            )  # drop the phenotype column before stacking

        d = data.stack().reset_index()  # collapse everything to one column
        d.columns = ['cellid', 'group', 'distance']
        d = pd.merge(
            d, pheno_df, left_on='cellid', right_index=True
        )  # bring back the imageid and phenotype

        # Convert columns to str
        for col in ['imageid', 'group', 'phenotype']:
            d[col] = d[col].astype(str)

        # Convert columns to categorical so that it drops unused categories
        for col in ['imageid', 'group', 'phenotype']:
            d[col] = d[col].astype('category')

        # re arrange the order based on from and to list provided
        if distance_to is not None:
            d['group'] = d['group'].cat.reorder_categories(distance_to)
            d = d.sort_values('group')

        # Plotting
        if method == 'numeric':
            if (
                x_axis is None
                and y_axis is None
                and facet_by is None
                and plot_type is None
            ):
                plot = sns.catplot(
                    data=d,
                    x="distance",
                    y="group",
                    col="imageid",
                    kind="boxen",
                    **kwargs,
                )
            else:
                plot = sns.catplot(
                    data=d, x=x_axis, y=y_axis, col=facet_by, kind=plot_type, **kwargs
                )

        if method == 'distribution':
            if (
                x_axis is None
                and y_axis is None
                and facet_by is None
                and plot_type is None
            ):
                print(d)
                plot = sns.displot(
                    data=d,
                    x="distance",
                    hue="imageid",
                    col="imageid",
                    kind="kde",
                    **kwargs,
                )
            else:
                plot = sns.displot(
                    data=d, x=x_axis, hue=y_axis, col=facet_by, kind=plot_type, **kwargs
                )
                plt.xlim(0,400)
    # Saving the figure if saveDir and fileName are provided
    if saveDir:
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        full_path = os.path.join(saveDir, fileName)
        plot.savefig(full_path, dpi=300)
        plt.close()
        print(f"Saved plot to {full_path}")
    else:
        plt.show()

    # return
    if return_data is True:
        return d

import scimap as sm
dist = sc.read('results/panel_g.scimap.h5ad')
tumor_dist = dist[dist.obs['pathology'] != 'N']
for celltype in ['B','CD8 T', 'Activated CD8 T', 'CD4 T', 'Activated CD4 T', 'Activated T reg']:
    spatial_distance_plot(
        tumor_dist, phenotype = 'celltype_activation', imageid='pathology',# 'GGO ID',
        distance_from = 'Epithelial-like',
        distance_to = celltype,
        method = 'distribution',  plot_type = 'kde', 
        saveDir= 'figures/spatial_distance/', fileName = f'spatial_distance_{celltype}.pdf')



for path in dist_h.obs['pathology'].cat.categories:
    for celltype in ['Fibroblast','B','CD4 T', 'CD8 T', 'Mast', 'Endothelial']:
        d = spatial_distance_plot(
            dist_h[dist_h.obs['pathology'] == path], phenotype = 'celltype_broad', imageid='pathology',# 'GGO ID',
            distance_from = 'Tumor-like',
            distance_to = celltype,
            method = 'distribution',  plot_type = 'kde', 
            saveDir= 'figures/spatial_distance/', fileName = f'spatial_distance_{celltype}_{path}.pdf', return_data = True)



# Start
distance_from = 'Tumor-like'
celltype_name = 'celltype_broad'
celltype_name = 'celltype_activation'
adata = dist # dist_h
adata = dist[dist.obs['pathology']!='N'] # dist_h

for celltype in adata.obs[celltype_name].cat.categories:
    # 'celltype_activation'
    diatance_map = adata.uns['spatial_distance'].copy()
    pheno_df = pd.DataFrame(
        {'imageid': adata.obs['pathology'], 'phenotype': adata.obs[celltype_name]}
    )  # image id and phenotype
    data = pd.merge(
        pheno_df, diatance_map, how='outer', left_index=True, right_index=True
    )  # merge with the distance map
    data = data[data['phenotype'] == celltype]  # subset the pheno of interest

    data = data.drop(
        ['phenotype', 'imageid'], axis=1
    )  # drop the phenotype column before stacking

    d = data.stack().reset_index()  # collapse everything to one column
    d.columns = ['cellid', 'group', 'distance']
    d = pd.merge(
        d, pheno_df, left_on='cellid', right_index=True
    )  # bring back the imageid and phenotype

    # Convert columns to str
    for col in ['imageid', 'group', 'phenotype']:
        d[col] = d[col].astype(str)

    # Convert columns to categorical so that it drops unused categories
    for col in ['imageid', 'group', 'phenotype']:
        d[col] = d[col].astype('category')

    # re arrange the order based on from and to list provided
    d = d.sort_values('group')

    plot = sns.displot(
        data=d,
        x="distance",
        # hue="imageid",
        col="group",
        kind="kde",
        common_norm = False,
        log_scale=True,
        facet_kws=dict(sharey=False)
    )
    plt.xlim(0,1000)
    plot.savefig(f'figures/spatial_distance/spatial_distance_{celltype}.pdf', dpi=300)
    plt.close()
