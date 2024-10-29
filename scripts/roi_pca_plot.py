import scanpy as sc
from glob import glob
import tifffile
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
import anndata
import matplotlib
import os
import imc_analysis as imc
from scripts.load_yaml import load

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False

# load data
metadata = load('metadata/ggo_config.yml')

if os.path.exists(metadata['PCA_ANNDATA_NAME']):
    density_all = sc.read(metadata['PCA_ANNDATA_NAME'])
else:
    adata_dict = dict()
    density_dict = dict()

    for panel in ['PANEL_G', 'PANEL_H']:
        adata_dict[panel] = sc.read(metadata[panel]['AnnData']['phenotyped_file_name'])
        del adata_dict[panel].obs['ROI_area']
        for panel in adata_dict:
            for celltype in metadata['CELLTYPES']:
                #density = compute_density(
                density = imc.tl.celltype_density(
                    adata_dict[panel],
                    celltype_key = celltype,
                    condition_keys = metadata['CONDITIONS'] + ['sample', 'description', 'GGO ID'])
                # density.obs.set_index('GGO ID', inplace = True)
            density_dict[panel] = density

    # remove non unique roi_id due to failure in scanning
    for panel in ['PANEL_G', 'PANEL_H']:
        density = density_dict[panel]
        if 'GGO ID' in density.obs:
            idx_keep = density.obs[['GGO ID']].drop_duplicates().index
            density = density[idx_keep]

            density.obs.set_index('GGO ID', inplace = True)
            density_dict[panel] = density
        density_dict[panel].obs_names = density.obs.index

    # collect coexisting samples
    coexistent = sorted(list(set(density_dict['PANEL_G'].obs.index).intersection(set(density_dict['PANEL_H'].obs.index))))

    for panel in ['PANEL_G', 'PANEL_H']:
        density_dict[panel] = density_dict[panel][coexistent]

    # build concatenated PCA
    density_all = anndata.concat([density_dict['PANEL_G'], density_dict['PANEL_H']], axis = 1, join = 'outer')
    density_all.obs = density_dict[panel].obs
    density_all.write(metadata['PCA_ANNDATA_NAME'])



# density_all.obs['radio'] = pd.Categorical(density_all.obs['radio'], categories = ['N', 'PNS', 'PS', 'S'])
# density_all.uns['radio_colors'] = metadata['pred_color']# ['#EEEEEE', '#AAAAAA', '#666666', '#111111'] # metadata['radio_color']

# density_all = density_all[~density_all.obs['radio'].isna()]
density_all.obs['pathology'] = density_all.obs['pathology'].replace({'Normal':'N'}) 
density_all.obs['pathology'] = pd.Categorical(density_all.obs['pathology'], categories = ['N', 'AIS', 'MIA', 'IAC'])
density_all.uns['pathology_colors'] = metadata['pathology_color']# ['#EEEEEE', '#AAAAAA', '#666666', '#111111'] # metadata['pathology_color']

density_all = density_all[~density_all.obs['pathology'].isna()]
density_all.obs['GGO ID'] = density_all.obs['description'].str[:6]
density_all.obs.index.name = None

# Perform PCA
# Flip x axis such that normal is on the negative side
sc.pp.scale(density_all)
sc.pp.pca(density_all)
sc.pp.neighbors(density_all)
sc.tl.umap(density_all)
density_all.obsm['X_pca'][:,0] = -density_all.obsm['X_pca'][:,0]

fig, ax = plt.subplots(1,1,figsize = (12,8), dpi = 300)
sc.pl.pca(
    density_all,
    color = 'pathology',
    size = 200,
    add_outline = True,
    outline_width = (0.15, 0.),
    annotate_var_explained = True,
    title = '',
    ax = ax,
    show = False
) 
sns.despine()
plt.tight_layout()
plt.savefig('figures/roi_pca_all_combined1.pdf', bbox_inches = 'tight')
plt.close()


fig, ax = plt.subplots(1,1,figsize = (12,8), dpi = 300)
sc.pl.umap(
    density_all,
    color = 'pathology',
    size = 200,
    add_outline = True,
    outline_width = (0.15, 0.),
    title = '',
    ax = ax,
    show = False
) 
sns.despine()
plt.tight_layout()
plt.savefig('figures/roi_umap_all_combined.pdf', bbox_inches = 'tight')
plt.close()


# sns.scatterplot(
#     density_all.obs,
#     x = 'PC1',
#     y = 'PC2',
#     hue = 'pathology',
#     color = density_all.uns['pathology_colors']
# )
# sns.set_theme()
# sns.despine()
# plt.tight_layout()
# plt.savefig('figures/roi_pca_all_combined_with_grids.pdf', bbox_inches = 'tight')
# plt.close()

# density_all.obs['PC1'] = density_all.obsm['X_pca'][:,0]
# density_all.obs['PC2'] = density_all.obsm['X_pca'][:,1]
# density_all.obs.to_csv('metadata/roi_pca_obs.csv')


# color by group
a = sc.read(metadata['patient_celltype_broad_clustered'])
#density_all.obs = density_all.obs.merge(a.obs[['Group']], left_on = 'GGO ID', right_index = True)
tmp = density_all.obs.merge(a.obs[['Group']], left_on = 'GGO ID', right_index = True)
tmp['Group'] = tmp['Group'].astype(str)
tmp.loc[tmp['pathology']=='N','Group'] = 'N'
tmp['Group'] = pd.Categorical(tmp['Group'], categories = ['N', 'Group1', 'Group2', 'Group3', 'Group4'], ordered = True)
da = density_all[tmp.index]
da.obs = tmp
da.uns['Group_colors'] = ['#2AC674', '#8DAC91', '#F9C8AB', '#DB4962', '#252F32']

fig, ax = plt.subplots(1,1,figsize = (6,5), dpi = 300)
sc.pl.pca(
    da,
    color = 'Group',
    size = 200,
    add_outline = True,
    outline_width = (0.15, 0.),
    annotate_var_explained = True,
    title = '',
    legend_loc = None,
    ax = ax,
    show = False
) 

# display mean and covariance matrix as ellipse on pca plot
for i,group in enumerate(da.obs['Group'].cat.categories):
    da_group = da[da.obs['Group'] == group]
    mean = da_group.obsm['X_pca'][:,:2].mean(axis = 0)
    cov = np.cov(da_group.obsm['X_pca'][:,:2].T)

    from matplotlib.patches import Ellipse

    eigenvalues, eigenvectors = np.linalg.eig(cov)
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
    width, height = 2 * np.sqrt(eigenvalues)
    ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle, facecolor=da.uns['Group_colors'][i], edgecolor = None, alpha = 0.3)
    plt.gca().add_patch(ellipse)


sns.despine()
plt.tight_layout()
plt.savefig('figures/roi_pca_group.pdf', bbox_inches = 'tight')
plt.close()

