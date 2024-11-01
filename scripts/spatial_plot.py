import squidpy as sq
import scanpy as sc

import matplotlib.pyplot as plt
import imc_analysis as imc

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')


adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['utag_labeled_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])


# for panel in ['PANEL_G', 'PANEL_H']:
#     print(f'Producing files for {panel}')
#     for roi in tqdm(adata_dict[panel].obs['roi'].unique()):
        
#         a = adata_dict[panel][adata_dict[panel].obs['roi'] == roi].copy()
        
#         sq.gr.spatial_neighbors(a, coord_type = 'generic', radius = 40)
        
#         id_, r, p = a.obs[['description', 'radio', 'pathology']].iloc[0]
        
#         sc.pl.spatial(
#             a,
#             color = ['uE_broad', 'celltype', 'celltype_broad'],
#             spot_size = 10,
#             edges = True,
#             neighbors_key = 'spatial_neighbors',
#             frameon = False,
#             show = False
#         )
        
#         os.makedirs(f'figures/spatial/{r}/', exist_ok = True)
#         plt.suptitle(f'ID: {id_}, Radiology: {r}, Histology: {p}')
#         plt.tight_layout()
#         plt.savefig(f'figures/spatial/{r}/{id_}_{panel}.pdf')
#         plt.close()


for panel in ['PANEL_G', 'PANEL_H']:
    print(f'Producing files for {panel}')
    #for roi in tqdm(adata_dict[panel].obs['roi'].unique()):
    for roi in ['10202022_CTMA123D_PanelG_1-08', '10252022_CTMA123D_PanelG_1-08']:
        
        a = adata_dict[panel][adata_dict[panel].obs['roi'] == roi].copy()
        
        if len(a):
            sq.gr.spatial_neighbors(a, coord_type = 'generic', radius = 40)
            
            id_, r, p = a.obs[['description', 'radio', 'pathology']].iloc[0]
            
            sc.pl.spatial(
                a,
                color = ['uE_broad', 'celltype', 'celltype_broad'],
                spot_size = 10,
                edges = True,
                neighbors_key = 'spatial_neighbors',
                frameon = False,
                show = False
            )
            
            os.makedirs(f'figures/spatial/{r}/', exist_ok = True)
            plt.suptitle(f'ID: {id_}, Radiology: {r}, Histology: {p}')
            plt.tight_layout()
            plt.savefig(f'figures/spatial/{r}/{id_}_{panel}.pdf')
            plt.close()