
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

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
#matplotlib.use('Agg')

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


adata_dict['PANEL_G'].obs['major_minor_ratio'] = adata_dict['PANEL_G'].obs['major_axis_length'] / adata_dict['PANEL_G'].obs['minor_axis_length']
adata_epi = adata_dict['PANEL_G'][adata_dict['PANEL_G'].obs['celltype_broad'].isin(['Epi. Prol.', 'Epi.-like'])]
adata_emt = adata_dict['PANEL_G'][adata_dict['PANEL_G'].obs['celltype_broad'].isin(['Epi.-like','Mesen.-like'])]
morph_features = ['area', 'major_minor_ratio']

sc.pl.violin(adata_epi, keys = morph_features, groupby = 'pathology', use_raw = False, jitter = False, inner = 'box')
sc.pl.violin(adata_emt, keys = morph_features, groupby = 'pathology', use_raw = False, jitter = False, inner = 'box')
