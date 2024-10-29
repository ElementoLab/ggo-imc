
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

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')


import yaml
c = 'celltype'
cluster_res = 'cluster_0.5'

metadata_filename = 'metadata/ggo_config.yml'

with open(metadata_filename, "r") as stream:
    try:
        metadata = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_file_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])
    adata_dict[p].obs['PID'] = adata_dict[p].obs['description'].str[:6]
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)


feature_count = dict()
for p in ['PANEL_G', 'PANEL_H']:
    feature_count[p] = adata_dict[p].obs[['PID','radio','pathology']].drop_duplicates()
feature_count = pd.concat(feature_count.values()).drop_duplicates()
print('Number of total patients:',len(feature_count['PID'].unique()))
for f in ['radio', 'pathology', ['radio', 'pathology']]:
    print(f, ":",feature_count.drop_duplicates()[f].value_counts())