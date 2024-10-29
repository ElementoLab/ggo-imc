
import os
from glob import glob
from pathlib import Path
import tifffile
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

from skimage.exposure import adjust_gamma
from skimage import filters
import scipy.ndimage as ndi
import scipy

import seaborn as sns
import matplotlib.pyplot as plt

import anndata
import warnings
warnings.simplefilter("ignore", UserWarning)

import yaml
metadata_filename = 'metadata/ggo_config.yml'

c = 'celltype'
cluster_res = 'cluster_1.0'

with open(metadata_filename, "r") as stream:
    try:
        metadata = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

print(f"Reading {metadata['panel_g_phenotyped_file_name']}...")
pg_phenotyped = sc.read(metadata['panel_g_phenotyped_file_name'])
print(f"Reading {metadata['panel_h_phenotyped_file_name']}...")
ph_phenotyped = sc.read(metadata['panel_h_phenotyped_file_name'])

pg_phenotyped.obs['PID'] = pg_phenotyped.obs['sample_y'].str[:-1]
ph_phenotyped.obs['PID'] = ph_phenotyped.obs['sample_y'].str[:-1]

pg_phenotyped.obs['celltype_broad'] = pg_phenotyped.obs['celltype'].str.split(' (', regex = False).str[0]
ph_phenotyped.obs['celltype_broad'] = ph_phenotyped.obs['celltype'].str.split(' (', regex = False).str[0]

pg_phenotyped.obs['ggo'] = pd.Categorical(pg_phenotyped.obs['ggo'], categories = ['Normal', 'Pure', '<25%', '25-50%', '>50%'])
pg_phenotyped.obs.loc[pg_phenotyped.obs['sample status'] == 'Normal','ggo'] = 'Normal'

ph_phenotyped.obs['ggo'] = pd.Categorical(ph_phenotyped.obs['ggo'], categories = ['Normal', 'Pure', '<25%', '25-50%', '>50%'])
ph_phenotyped.obs.loc[ph_phenotyped.obs['sample status'] == 'Normal','ggo'] = 'Normal'

solidity_mapper = {c: i for i, c in enumerate(ph_phenotyped.obs['ggo'].cat.categories)}
ph_phenotyped.obs['solidity_num'] = ph_phenotyped.obs['ggo'].replace(solidity_mapper)
pg_phenotyped.obs['solidity_num'] = pg_phenotyped.obs['ggo'].replace(solidity_mapper)


pg_phenotyped.obs['Gender'] = pd.Categorical(pg_phenotyped.obs['Gender'].str.upper().replace(metadata['gender_map']))
ph_phenotyped.obs['Gender'] = pd.Categorical(ph_phenotyped.obs['Gender'].str.upper().replace(metadata['gender_map']))

pg_phenotyped.write(metadata['panel_g_phenotyped_labeled_file_name'])
ph_phenotyped.write(metadata['panel_h_phenotyped_labeled_file_name'])