import os
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import imc_analysis as imc

# ─────────────────────────────────────────────────────────────
# Load metadata
# ─────────────────────────────────────────────────────────────
metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

# ─────────────────────────────────────────────────────────────
# Load AnnData per panel
# ─────────────────────────────────────────────────────────────
def load_adata_panel(panel_key: str):
    panel_meta = metadata[panel_key]['AnnData']
    print(f"Reading {panel_meta['phenotyped_umap_name']}...")
    return sc.read(
        panel_meta['phenotyped_umap_name'],
        backup_url=panel_meta['backup_url']
    )

adata_dict = {
    panel: load_adata_panel(panel) for panel in ['PANEL_G', 'PANEL_H']
}

# ─────────────────────────────────────────────────────────────
# Myeloid analysis helper
# ─────────────────────────────────────────────────────────────
def extract_and_cluster_myeloids(
    adata: sc.AnnData,
    markers: list,
    celltype_filter: list,
    output_path: str,
    backup_url: str,
    key_added: str = 'mac_cluster',
    override: bool = False
):
    idx = adata.obs['celltype'].isin(celltype_filter)
    adata_myeloid = adata[idx, markers].copy()

    if 'celltype_broad' in adata_myeloid.obs:
        adata_myeloid.obs['celltype_broad'] = (
            adata_myeloid.obs['celltype_broad']
            .astype(str)
            .replace({'Neutrophil': 'Monocyte'})
        )
        adata_myeloid.obs['celltype_broad'] = pd.Categorical(
            adata_myeloid.obs['celltype_broad'],
            categories=['Monocyte', 'Macrophage', 'Mast', 'NK', 'PMN-MDSC']
        )

    if os.path.exists(output_path) and not override:
        return sc.read(output_path, backup_url=backup_url)

    sc.pp.scale(adata_myeloid)
    sc.tl.pca(adata_myeloid)
    sc.pp.neighbors(adata_myeloid)
    sc.tl.umap(adata_myeloid)
    sc.tl.leiden(adata_myeloid, key_added=key_added, resolution = 0.5)

    # Uncomment to write result
    # adata_myeloid.write(output_path)

    return adata_myeloid

# ─────────────────────────────────────────────────────────────
# Panel H Myeloid Analysis
# ─────────────────────────────────────────────────────────────
h_markers = [
    'VISTA', 'CD14', 'CD16', 'CD68', 'CD163', 'HLADR', 'CD33',
    'IL1beta', 'IL12p40', 'IL23p19',
]
h_celltypes = ['Mac. (CD68+, CD163-)', 'Mac. (CD163+)', 'NK', 'Neut.']
panel_h_meta = metadata['PANEL_H']['AnnData']

h_myeloid = extract_and_cluster_myeloids(
    adata=adata_dict['PANEL_H'],
    markers=h_markers,
    celltype_filter=h_celltypes,
    output_path=panel_h_meta['myeloids'],
    backup_url=panel_h_meta['myeloids_url'],
    override=True,
)

# ─────────────────────────────────────────────────────────────
# Panel G Myeloid Analysis
# ─────────────────────────────────────────────────────────────
g_markers = ['HLADR', 'CD14', 'CD16', 'CD68', 'CD163', 'CD206', 'PDL1']
g_celltypes = ['Mac. (CD163+)', 'Mac. (CD163+, CD206+)', 'NK']
panel_g_meta = metadata['PANEL_G']['AnnData']

g_myeloid = extract_and_cluster_myeloids(
    adata=adata_dict['PANEL_G'],
    markers=g_markers,
    celltype_filter=g_celltypes,
    output_path=panel_g_meta['myeloids'],
    backup_url=panel_g_meta['myeloids_url'],
    override=True,
)


sc.pl.umap(h_myeloid, color = ['mac_cluster', 'pathology'] + h_myeloid.var.index.tolist(), use_raw = False)


mac_map = {
    '0': 'M2-like TAMs (HLA-DR+, CD14+, CD16+, CD163+)',
    '1': 'M-MDSCs (HLA-DR-, CD14+, CD16-)',
    '2': 'Immature Myeloid Cells (HLA-DR+, CD14+, CD16+, VISTA+, CD33+)',
    '3': 'Immature Myeloid Cells (HLA-DR+, CD14+, CD16+, VISTA+, CD33+)',
    '4': 'Alveolar Macrophages (HLA-DR+, CD14-, CD16+, CD163+)', #'M-MDSCs (HLA-DR+, CD14+, CD16-, CD68+)',
    '5': 'M-MDSCs (HLA-DR-, CD14+, CD16-, VISTA+, IL23p19+)',
    '6': 'Alveolar Macrophages (HLA-DR+, CD14-, CD16+, CD163+)',
    '8': 'Immature Myeloid Cells (HLA-DR+, CD14+, CD16+, VISTA+, CD33+)',
    '9': 'Immature Myeloid Cells (HLA-DR+, CD14+, CD16+, VISTA+, CD33+)'
}

h_myeloid.obs['mac_cluster_celltype'] = h_myeloid.obs['mac_cluster'].replace(mac_map)
h_myeloid = h_myeloid[h_myeloid.obs['mac_cluster_celltype'].isin(mac_map.values())].copy()
h_myeloid.obs['mac_cluster_celltype_broad'] = h_myeloid.obs['mac_cluster_celltype'].str.split(r' \(').str[0]
h_myeloid.obs['mac_cluster_celltype_broad'] = pd.Categorical(
    h_myeloid.obs['mac_cluster_celltype_broad'],
    categories = ['Alveolar Macrophages', 'M-MDSCs', 'Immature Myeloid Cells', 'M2-like TAMs'],
    ordered = True
)
sc.pl.umap(h_myeloid, color = ['mac_cluster_celltype_broad', 'mac_cluster_celltype'] + h_myeloid.var.index.tolist(), use_raw = False)


for celltype in ['mac_cluster_celltype', 'mac_cluster_celltype_broad']:
    # computing celltype density
    density = imc.tl.celltype_density(
        h_myeloid,
        celltype = celltype,
        condition_keys = ['pathology'])

    # statistical testing
    imc.tl.grouped_mwu_test(
        density,
        condition_keys = ['pathology']
    )

    for pval_form in ['star', 'sci_not']:
        imc.pl.plot_mwu(
            density,
            kind = 'box-line',
            save_dir=f'figures/mac_densities_pathology/{celltype}',
            pval_form=pval_form
        )


adata = adata_dict['PANEL_H'].copy()
adata.obs['celltype'] = adata.obs['celltype'].astype(str)
adata.obs.loc[h_myeloid.obs.index, 'celltype'] = h_myeloid.obs['mac_cluster_celltype_broad']

# computing celltype density
density = imc.tl.celltype_density(
    adata,
    celltype = 'celltype',
    condition_keys = ['pathology'])

imc.pl.celltype_heatmap(adata, cluster_ids = ['celltype'], out_dir='figures/celltype/heatmap/', var_names = metadata['PANEL_H']['var_celltype_groups'])

# statistical testing
imc.tl.grouped_mwu_test(
    density,
    condition_keys = ['pathology']
)

for pval_form in ['star', 'sci_not']:
    imc.pl.plot_mwu(
        density,
        kind = 'box-line',
        save_dir=f'figures/mac_densities_pathology/{celltype}_all/',
        pval_form=pval_form
    )
