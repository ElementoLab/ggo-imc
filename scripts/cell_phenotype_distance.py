
import os
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt

import scimap as sm

import warnings
warnings.simplefilter("ignore", UserWarning)

#import imc_analysis as imc
from scripts.load_yaml import load

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

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

dist = spatial_distance(adata, phenotype = 'celltype', imageid='roi')
dist

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



dist_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    for cond in adata_dict[p].obs['pathology'].cat.categories:
        adata = adata_dict[p][adata_dict[p].obs['pathology'] == cond].copy()

        dist_dict[(p,cond)] = spatial_distance(adata, subset = ['03132023_CTMA123G_panelH_1-02'], phenotype = 'celltype_broad', imageid='roi')

['Tumor-like', 'Epithelial-like', 'Endothelial', 'Fibroblast']

celltype = 'Tumor-like' # 'Epithelial-like' # 'Tumor-like'

for k in dist_dict:

    data = dist_dict[k]
    data[data.obs['celltype_broad'] == celltype]
    print(k)
    print(data.uns['spatial_distance'].describe())