import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt

import os
import imc_analysis as imc

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

# def spatial_distance (adata,
#                       x_coordinate='X_centroid',
#                       y_coordinate='Y_centroid',
#                       z_coordinate=None,
#                       phenotype='phenotype',
#                       subset=None,
#                       imageid='imageid',
#                       verbose=True,
#                       label='spatial_distance'):

#     import pandas as pd
#     from sklearn.neighbors import BallTree
#     from joblib import Parallel, delayed
#     import itertools

#     def spatial_distance_internal (adata_subset,x_coordinate,y_coordinate,z_coordinate,
#                                    phenotype,subset,imageid,label):

#         if verbose:
#             print("Processing Image: " + str(adata_subset.obs[imageid].unique()[0]))
#         # Create a dataFrame with the necessary inforamtion
#         data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})

#         # Function to identify shortest distance for each phenotype of interest
#         def distance (pheno):
#             pheno_interest = data[data['phenotype'] == pheno]

#             if len(pheno_interest) == 0:
#                 return []
#             # print(pheno)
#             # Build the ball-tree for search space
#             tree = BallTree(pheno_interest[['x','y']], metric='euclidean') 
#             # Calculate shortest distance (if statement to account for K)
#             if len(pheno_interest) > 1:
#                 dist, ind = tree.query(data[['x','y']], k=2, return_distance= True)
#                 dist = pd.DataFrame(dist)
#                 dist.loc[dist[0] == 0, 0]  = dist[1]
#                 dist = dist[0].values
#             else:
#                 dist, ind = tree.query(data[['x','y']], k=1, return_distance= True)
#                 dist = list(itertools.chain.from_iterable(dist))
            
#             # print(dist)
#             return dist

#         # Run in parallel for all phenotypes
#         phenotype_list = list(data['phenotype'].unique())
#         # print(phenotype_list)
#         # Apply function
#         final_dist = Parallel(n_jobs=-1)(delayed(distance)(pheno=i) for i in phenotype_list)     
#         final_dist = pd.DataFrame(final_dist, index = phenotype_list, columns = data.index).T

#         return final_dist

#     # subset a particular subset of cells if the user wants else break the adata into list of anndata objects
#     adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]

#     # Apply function to all images and create a master dataframe
#     # Create lamda function 
#     r_spatial_distance_internal = lambda x: spatial_distance_internal (adata_subset=x,
#                                                                        x_coordinate=x_coordinate,y_coordinate=y_coordinate, z_coordinate=z_coordinate,
#                                                                        phenotype=phenotype,subset=subset,imageid=imageid,label=label) 
#     all_data = list(map(r_spatial_distance_internal, adata_list)) # Apply function 

#     # Merge all the results into a single dataframe
#     result = []
#     for i in range(len(all_data)):
#         result.append(all_data[i])
#     result = pd.concat(result, join='outer')

#     # Add to anndata
#     adata.uns[label] = result

#     # return
#     return adata

# def adata_to_dist_df(
#     adata: anndata.AnnData,
#     dist_key: str = 'spatial_distance',
#     image_id: str = 'roi',
#     phenotype: str = 'celltype',
#     celltype: str = '', # celltype of interest. if not in phenotype, ignore
#     ):

#     diatance_map = adata.uns['spatial_distance'].copy()
#     pheno_df = pd.DataFrame(
#         {'imageid': adata.obs[image_id], 'phenotype': adata.obs[celltype_name]}
#     )  # image id and phenotype
#     data = pd.merge(
#         pheno_df, diatance_map, how='outer', left_index=True, right_index=True
#     )  # merge with the distance map

#     if celltype in data['phenotype'].cat.categories:
#         data = data[data['phenotype'] == celltype]  # subset the pheno of interest
#         print(f'Calculating distance density for {celltype}')
#     else:
#         print('Calculating distance density for all cells')

#     data = data.drop(
#         ['phenotype', 'imageid'], axis=1
#     )  # drop the phenotype column before stacking

#     d = data.stack().reset_index()  # collapse everything to one column
#     d.columns = ['cellid', 'group', 'distance']
#     d = pd.merge(
#         d, pheno_df, left_on='cellid', right_index=True
#     )  # bring back the imageid and phenotype

#     # Convert columns to str
#     for col in ['imageid', 'group', 'phenotype']:
#         d[col] = d[col].astype(str)

#     # Convert columns to categorical so that it drops unused categories
#     for col in ['imageid', 'group', 'phenotype']:
#         d[col] = d[col].astype('category')

#     # re arrange the order based on from and to list provided
#     d = d.sort_values('group')
#     return d
    
# # library
# def spatial_distance_plot(
#     adata,
#     spatial_distance='spatial_distance',
#     phenotype='phenotype',
#     imageid='imageid',
#     log=False,
#     method='heatmap',
#     heatmap_summarize=True,
#     heatmap_na_color='grey',
#     heatmap_cmap='vlag_r',
#     heatmap_row_cluster=False,
#     heatmap_col_cluster=False,
#     heatmap_standard_scale=0,
#     distance_from=None,
#     distance_to=None,
#     x_axis=None,
#     y_axis=None,
#     facet_by=None,
#     plot_type=None,
#     return_data=False,
#     subset_col=None,
#     subset_value=None,
#     fileName='spatial_distance.pdf',
#     saveDir=None,
#     **kwargs,
# ):
    
#     import pandas as pd
#     import matplotlib
#     import matplotlib.pyplot as plt
#     import numpy as np
#     import seaborn as sns
#     # set color for heatmap
#     # cmap_updated = matplotlib.cm.get_cmap(heatmap_cmap)
#     cmap_updated = matplotlib.colormaps[heatmap_cmap]
#     cmap_updated.set_bad(color=heatmap_na_color)

#     # Copy the spatial_distance results from anndata object
#     try:
#         diatance_map = adata.uns[spatial_distance].copy()
#     except KeyError:
#         raise ValueError(
#             'spatial_distance not found- Please run sm.tl.spatial_distance first'
#         )

#     # subset the data if user requests
#     if subset_col is not None:
#         if isinstance(subset_value, str):
#             subset_value = [subset_value]
#         # find the cell names to be subsetted out
#         obs = adata.obs[[subset_col]]
#         cells_to_subset = obs[obs[subset_col].isin(subset_value)].index

#         # subset the diatance_map
#         diatance_map = diatance_map.loc[
#             diatance_map.index.intersection(cells_to_subset)
#         ]
#         # diatance_map = diatance_map.loc[cells_to_subset]

#     # Convert distance to log scale if user requests
#     if log is True:
#         diatance_map = np.log1p(diatance_map)

#     # Method
#     if method == 'heatmap':
#         if heatmap_summarize is True:
#             # create the necessary data
#             data = pd.DataFrame({'phenotype': adata.obs[phenotype]})
#             data = pd.merge(
#                 data, diatance_map, how='outer', left_index=True, right_index=True
#             )  # merge with the distance map
#             k = data.groupby(
#                 ['phenotype'], observed=False
#             ).mean()  # collapse the whole dataset into mean expression
#             d = k[k.index]
#         else:
#             # create new naming scheme for the phenotypes
#             non_summary = pd.DataFrame(
#                 {'imageid': adata.obs[imageid], 'phenotype': adata.obs[phenotype]}
#             )
#             non_summary['imageid'] = non_summary['imageid'].astype(
#                 str
#             )  # convert the column to string
#             non_summary['phenotype'] = non_summary['phenotype'].astype(
#                 str
#             )  # convert the column to string
#             non_summary['image_phenotype'] = non_summary['imageid'].str.cat(
#                 non_summary['phenotype'], sep="_"
#             )
#             # Merge distance map with phenotype
#             data = pd.DataFrame(non_summary[['image_phenotype']])
#             data = pd.merge(
#                 data, diatance_map, how='outer', left_index=True, right_index=True
#             )
#             k = data.groupby(['image_phenotype'], observed=False).mean()
#             d = k.sort_index(axis=1)
#         # Generate the heatmap
#         mask = d.isnull()  # identify the NAN's for masking
#         d = d.fillna(0)  # replace nan's with 0 so that clustering will work
#         # Heatmap
#         plot = sns.clustermap(
#             d,
#             cmap=heatmap_cmap,
#             row_cluster=heatmap_row_cluster,
#             col_cluster=heatmap_col_cluster,
#             mask=mask,
#             standard_scale=heatmap_standard_scale,
#             **kwargs,
#         )
#     else:

#         # condition-1
#         if distance_from is None and distance_to is None:
#             raise ValueError(
#                 'Please include distance_from and/or distance_to parameters to use this method'
#             )

#         # condition-2
#         if distance_from is None and distance_to is not None:
#             raise ValueError('Please `distance_from` parameters to use this method')

#         # condition-3
#         if distance_to is not None:
#             # convert input to list if needed
#             if isinstance(distance_to, str):
#                 distance_to = [distance_to]

#         # Start
#         pheno_df = pd.DataFrame(
#             {'imageid': adata.obs[imageid], 'phenotype': adata.obs[phenotype]}
#         )  # image id and phenotype
#         data = pd.merge(
#             pheno_df, diatance_map, how='outer', left_index=True, right_index=True
#         )  # merge with the distance map
#         data = data[data['phenotype'] == distance_from]  # subset the pheno of interest

#         if distance_to is not None:
#             data = data[
#                 distance_to
#             ]  # drop columns that are not requested in distance_to
#         else:
#             data = data.drop(
#                 ['phenotype', 'imageid'], axis=1
#             )  # drop the phenotype column before stacking

#         d = data.stack().reset_index()  # collapse everything to one column
#         d.columns = ['cellid', 'group', 'distance']
#         d = pd.merge(
#             d, pheno_df, left_on='cellid', right_index=True
#         )  # bring back the imageid and phenotype

#         # Convert columns to str
#         for col in ['imageid', 'group', 'phenotype']:
#             d[col] = d[col].astype(str)

#         # Convert columns to categorical so that it drops unused categories
#         for col in ['imageid', 'group', 'phenotype']:
#             d[col] = d[col].astype('category')

#         # re arrange the order based on from and to list provided
#         if distance_to is not None:
#             d['group'] = d['group'].cat.reorder_categories(distance_to)
#             d = d.sort_values('group')

#         # Plotting
#         if method == 'numeric':
#             if (
#                 x_axis is None
#                 and y_axis is None
#                 and facet_by is None
#                 and plot_type is None
#             ):
#                 plot = sns.catplot(
#                     data=d,
#                     x="distance",
#                     y="group",
#                     col="imageid",
#                     kind="boxen",
#                     **kwargs,
#                 )
#             else:
#                 plot = sns.catplot(
#                     data=d, x=x_axis, y=y_axis, col=facet_by, kind=plot_type, **kwargs
#                 )

#         if method == 'distribution':
#             if (
#                 x_axis is None
#                 and y_axis is None
#                 and facet_by is None
#                 and plot_type is None
#             ):
#                 print(d)
#                 plot = sns.displot(
#                     data=d,
#                     x="distance",
#                     hue="imageid",
#                     col="imageid",
#                     kind="kde",
#                     **kwargs,
#                 )
#             else:
#                 plot = sns.displot(
#                     data=d, x=x_axis, hue=y_axis, col=facet_by, kind=plot_type, **kwargs
#                 )
#                 plt.xlim(0,400)
#     # Saving the figure if saveDir and fileName are provided
#     if saveDir:
#         if not os.path.exists(saveDir):
#             os.makedirs(saveDir)
#         full_path = os.path.join(saveDir, fileName)
#         plot.savefig(full_path, dpi=300)
#         plt.close()
#         print(f"Saved plot to {full_path}")
#     else:
#         plt.show()

#     # return
#     if return_data is True:
#         return d

adata_dict = dict()

for p in ['PANEL_G', 'PANEL_H']:
    adata_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])

dist_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    dist_dict[p] = imc.tl.spatial_distance(adata_dict[p], subset = ['03132023_CTMA123G_panelH_1-02'], celltype = 'celltype_broad', roi_key='GGO ID')

dist_dict[p]
adata_dict[p].obs[['GGO ID', 'pathology']]
dist_dict[p].obs['imageid']
# distance_from = 'Tumor-like'
celltype_name = 'celltype_broad'

dist_distribution = dict()

# establishing baseline for proximity
for p in dist_dict:
    # adata = dist_dict[p] # dist_h
    adata = dist_dict[p][dist_dict[p].obs['pathology']!='N'] # dist_h

    d = imc.utils.adata_to_dist_df(adata, roi_key = 'GGO ID', celltype = 'celltype_broad')
    d = d[d['group'] != d['phenotype']]
    d['panel'] = p
    dist_distribution[p] = d

def categorize(value):
    if value < 16.5:
        return 'neighboring'
    elif value < 26:
        return 'proximal'
    else:
        return 'distal'

dist = pd.concat(dist_distribution.values())
dist['proximity'] = dist['distance'].apply(categorize)
distance = dist.groupby('proximity')[['cellid']].count()
distance_prop = distance/distance.sum()
distance_prop_str = f''
for d_ in ['neighboring', 'proximal', 'distal']:
    distance_prop_str += f"{d_} : {distance_prop.loc[d_].values[0]:>5.2%}\n"

plot = sns.displot(
    data=dist,
    x="distance",
    hue="panel",
    kind="kde",
    common_norm = False,
    log_scale=True,
    facet_kws=dict(sharey=False)
)
plt.xlim(5,1000)
plt.vlines(x=16.5,ymin=0, ymax=1, color = 'darkred',linestyles = 'dashed')
plt.vlines(x=26,ymin=0, ymax=1, color = 'darkblue',linestyles = 'dashed')
plt.text(s = distance_prop_str, y = 1, x = 1000, va = 'top', ha = 'right')
sns.despine()
plot.savefig(f'figures/spatial_distance_all.pdf', dpi=300)
plt.close()

# for celltypes
for p in dist_dict:
    # adata = dist_dict[p] # dist_h
    adata = dist_dict[p][dist_dict[p].obs['pathology']!='N'] # dist_h

    for celltype in adata.obs[celltype_name].cat.categories:
        
        d = adata_to_dist_df(adata, image_id = 'GGO ID', phenotype = 'celltype_broad', celltype = celltype)
        d['proximity'] = d['distance'].apply(categorize)
        
        plot = sns.displot(
            data=d,
            x="distance",
            # hue="imageid",
            col="group",
            kind="kde",
            common_norm = False,
            log_scale=True,
            facet_kws=dict(sharey=False),
        )

        plt.xlim(5,1000)

        for i, ax in enumerate(plot.axes.flat):
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            ax.vlines(x=16.5,ymin=0, ymax=ymax, color = 'darkred',linestyles = 'dashed')
            ax.vlines(x=26,ymin=0, ymax=ymax, color = 'darkblue',linestyles = 'dashed')

            cell_type_ = d['group'].cat.categories[i]
            d_ = d[d['group'] == cell_type_].copy()
            distance = d_.groupby('proximity')[['cellid']].count()
            distance_prop = distance/distance.sum()

            distance_prop_str = f''

            for d_ in ['neighboring', 'proximal', 'distal']:
                distance_prop_str += f"{d_} : {distance_prop.loc[d_].values[0]:>5.2%}\n"
                dist_df_dict[d_].loc[celltype, cell_type_] = distance_prop.loc[d_].values[0]

            ax.text(s = distance_prop_str, y = ymax, x = xmax, va = 'top', ha = 'right')
            ax.set_xlabel('distance (μm)')

        sns.despine()
        plot.axes.flat[0].set_ylabel(f'{celltype} cell density distribution')
        os.makedirs(f'figures/{p}/spatial_distance/', exist_ok = True)
        plot.savefig(f'figures/{p}/spatial_distance/spatial_distance_{celltype}.pdf', dpi=300)
        plt.close()


# distance distributions
celltypes = adata.obs[celltype_name].cat.categories
dist_df_dict = dict()
for d_ in ['neighboring', 'proximal', 'distal']:
    dist_df_dict[d_] = pd.DataFrame(index = celltypes, columns = celltypes)

for p in dist_distribution:
    dist = dist_distribution[p]
    dist['proximity'] = dist['distance'].apply(categorize)

    dist.groupby(['group', 'phenotype', 'proximity'])
    for d_ in ['neighboring', 'proximal', 'distal']:
        dist_df_dict[d_].loc[celltype, cell_type_] = distance_prop.loc[d_].values[0]


dist_df_dict['neighboring+proximal'] = dist_df_dict['neighboring'] + dist_df_dict['proximal']
sns.clustermap(dist_df_dict['neighboring+proximal'].astype(float), cmap = 'Spectral_r', method = '')
plt.show()
sns.heatmap(dist_df_dict['neighboring+proximal'].astype(float), cmap = 'Spectral_r')
plt.show()





# dist = spatial_distance(adata, phenotype = 'celltype', imageid='roi')
# dist

# adata_dict = dict()

# for p in ['PANEL_G', 'PANEL_H']:
#     adata_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])

# all_dist_dict = dict()
# for p in ['PANEL_G', 'PANEL_H']:
#     all_dist_dict[p] = spatial_distance(adata_dict[p], subset = ['03132023_CTMA123G_panelH_1-02'], phenotype = 'celltype_broad', imageid='GGO ID')

# dist_dict = dict()
# for p in ['PANEL_G', 'PANEL_H']:
#     dist_dict[p] = spatial_distance(adata, subset = ['03132023_CTMA123G_panelH_1-02'], phenotype = 'celltype_broad', imageid='GGO ID')

# ['Tumor-like', 'Epithelial-like', 'Endothelial', 'Fibroblast', 'Mast']

# celltype = 'Tumor-like' # 'Epithelial-like' # 'Tumor-like'

# for p in dist_dict:
#     data = dist_dict[p][dist_dict[p].obs['celltype_broad'] == celltype]
#     fig, axes = plt.subplots(
#         1, len(data.uns['spatial_distance'].columns),
#         figsize = (3 * len(data.uns['spatial_distance'].columns), 3),
#         dpi = 300
#     )
#     for i, cell in enumerate(data.uns['spatial_distance'].columns):
#         sns.kdeplot(data.uns['spatial_distance'][[cell]], color = data.obs['pathology'], ax = axes[i])
#         axes[i].set_xscale('log')
#         axes[i].set_xlim(1,1000)
#         axes[i].set_title(cell)
#         if i!=0:
#             axes[i].set_ylabel('')
#         if i < len(data.uns['spatial_distance'].columns) - 1:
#             axes[i].legend().remove()
#     plt.show()