from glob import glob
from tqdm import tqdm
import yaml
import pandas as pd
import anndata
import scanpy as sc
import pandas as pd
from tqdm import tqdm

def phenotyping(
    a: anndata.AnnData,
    channels_include = None,
    channels_exclude = None,
    filter_cells: bool = True,
    z_score: bool = True,
    z_score_per: str = "roi",
    z_score_cap: float = 3.0,
    remove_batch: bool = True,
    batch_variable: str = "sample",
    dim_res_algos = ["umap",],
    clustering_method: str = "leiden",
    clustering_resolutions =[0.5],
    save_name = 'quant'
):
    if isinstance(a, str):
        print(f"Reading h5ad file: '{a}'.")
        a = sc.read(a)
    if remove_batch:
        if a.obs[batch_variable].nunique() <= 1:
            print(
                "Batch correction not possible as only one batch detected. "
                "Check `batch_variable` keyord argument."
            )
            remove_batch = False
    if "sample" not in a.obs.columns:
        a.obs["sample"] = a.obs["roi"].str.extract(r"(.*)-\d+")[0].fillna("")
    if a.raw is None:
        a.raw = a
    # Add morphological variables to obs
    sel = a.var.index.str.contains(r"\(")
    v = a.var.index[~sel]
    for col in v:
        a.obs[col] = a[:, col].X.flatten().tolist()
    a = a[:, sel]
    # Filter out channels
    if channels_exclude is not None:
        a = a[:, ~a.var.index.isin(channels_exclude)]
    if channels_include is not None:
        a = a[:, channels_include]
    a = a.copy()
    # # reduce DNA chanels to one, and move to obs
    dnas = a.var.index[a.var.index.str.contains(r"DNA\d")]
    a.obs["DNA"] = a[:, dnas].X.mean(1)
    a = a[:, ~a.var.index.isin(dnas)]
    # Filter out cells
    if filter_cells:
        if "solidity" not in a.obs.columns:
            print(
                "Could not filter cells based on solidity likely because morphological quantification was not performed!"
            )
        else:
            exclude = a.obs["solidity"] == 1
            p = (exclude).sum() / a.shape[0] * 100
            print(f"Filtered out {exclude.sum()} cells ({p:.2f} %)")
    # Scaling/Normalization
    print("Performing data scaling/normalization.")
    sc.pp.log1p(a)
    if z_score:
        _ads = list()
        for roi_name in a.obs["roi"].unique():
            try:
                a2 = a[a.obs["roi"] == roi_name, :].copy()
                sc.pp.scale(a2, max_value=z_score_cap)
                a2.X[a2.X < -z_score_cap] = -z_score_cap
                # print(a2.X.min(), a2.X.max())
                _ads.append(a2)
            except ZeroDivisionError:
                print(f'ZeroDivisionError: skipping {roi_name}.')
        a = anndata.concat(_ads)
        sc.pp.scale(a)
    if remove_batch:
        sc.pp.pca(a)
        sc.external.pp.harmony_integrate(a, batch_variable, max_iter_harmony = 20)
    a.write(f'results/{save_name}.harmony.h5ad')
    # Dimensionality reduction
    print("Performing dimensionality reduction.")
    sc.pp.pca(a)
    if remove_batch:
        sc.external.pp.bbknn(a, batch_key=batch_variable)
    else:
        sc.pp.neighbors(a)

    
    # Clustering
    print("Performing clustering.")
    if clustering_method == "leiden":
        for res in clustering_resolutions:
            sc.tl.leiden(a, resolution=res, key_added=f"cluster_{res}")
            a.obs[f"cluster_{res}"] = pd.Categorical(
                a.obs[f"cluster_{res}"].astype(int) + 1
            )
        a.write(f'results/{save_name}.harmony.leiden.h5ad')
    elif clustering_method == "parc":
        for res in clustering_resolutions:
            p = PARC(
                a.X,
                neighbor_graph=a.obsp["connectivities"],
                random_seed=42,
                resolution_parameter=res,
            )
            p.run_PARC()
            a.obs[f"cluster_{res}"] = pd.Categorical(pd.Series(p.labels) + 1)
        a.write(f'results/{save_name}.harmony.parc.h5ad')
    print("Finished phenotyping.")
    if "umap" in dim_res_algos:
        sc.pp.subsample(a, n_obs=30000)
        sc.pp.neighbors(a)
        sc.tl.umap(a, gamma=3)
        a.write(f'results/{save_name}.harmony.umap.h5ad')
    return a

features_extract = [
	'roi',
	'description',
	'sample',
	'Radiology',
	'SOLID COMP',
	'clinical tumor SIZE (cm)',
	'Ratio (100%)',
	'Race',
	'Smoking',
	'Smoking Status',
	'Packyrs',
	'Race.1',
	'Age',
	'Gender',
	'PATH STAGE',
	'pred',
	'MM Class',
	'minor',
	'INV SIZE (MM)',
	'Dense scar (MM)',
	'Comment',
	'Molecular',
	'aah',
	'TMA: GGO# (T)',
	'TMA: GGO#(N)',
	'GGO ID',
	'core',
	'panel',
]

exclude_channel = [
	'80ArAr(ArAr80)',
	'131Xe(Xe131)',
	'138Ba(Ba138)',
	'190BCKG(BCKG190)',
	'<EMPTY>(Pt194)',
	'<EMPTY>(Pt195)',
	'<EMPTY>(Pt196)',
	'<EMPTY>(Pt198)',
	'208Pb(Pb208)',
]

panel_map = {
	'Panel_G': 'PanelG',
	'Panel_H': 'PanelH'
}

# read yaml files
yaml_files = glob('data/*.yaml')
yaml_metadata = []

for f in tqdm(yaml_files):
	with open(f,'r') as file:
		yaml_content = yaml.safe_load(file)
		acq = pd.DataFrame(yaml_content['acquisitions']).T
		acq['source_path'] = acq['source_path'].str.replace('added_data/','data/', regex = False)
		acq.index = acq['source_path'].str.replace('data/','', regex = False).str.replace('.mcd','', regex = False) + '-' + acq.index.astype(str).str.zfill(2)
		yaml_metadata.append(
			acq
		)


metadata_filename = 'metadata/ggo_config.yml'

with open(metadata_filename, "r") as stream:
	try:
		metadata = yaml.safe_load(stream)
	except yaml.YAMLError as exc:
		print(exc)


# create csv out of yaml metadata
df = pd.concat(yaml_metadata)
df['description'] = df['description'].replace('GG0', 'GGO', regex = True)
df['sample'] = df['description'].str.split('_').str[0]
df = df.reset_index(names = 'roi')
df.to_csv('metadata/roi_map.csv')

# filter samples with GGO ID
df[~df['description'].str.contains('GGO')].to_csv('metadata/no_ggo_id_data.csv')


clin = pd.read_excel(metadata['clinical_data'])
T = df.merge(clin, left_on = 'sample', right_on = 'TMA: GGO# (T)')
T['GGO ID'] = T['description']
N = df.merge(clin, left_on = 'sample', right_on = 'TMA: GGO#(N)')
N['GGO ID'] = N['description']
mapped_df = pd.concat([T,N])

mapped_df['core'] = mapped_df['source_path'].str.split('_').str[1]
mapped_df['panel'] = mapped_df['source_path'].replace(panel_map, regex = True).str.contains('PanelG')
mapped_df['panel'] = mapped_df['panel'].replace({True:'G',False:'H'})

print(f'Mapped {len(mapped_df)}/{len(df)} ({round(len(mapped_df)/len(df),2) * 100}%) samples.')


# save results
mapped_df.to_csv('metadata/roi_mapped_clinical.csv')
df[~df['roi'].isin(mapped_df['roi'])].to_csv('metadata/roi_no_matched_clinical_annotation.csv')


quant_file_list = glob('processed/quantification/quantification_*_full.h5ad')
quant_file_exists = ('processed/quantification/quantification_' + mapped_df['roi'] + '_full.h5ad').isin(quant_file_list)
mapped_df[~quant_file_exists].to_csv('metadata/no_quant_found.csv')

mapped_df = mapped_df[quant_file_exists]
dup = mapped_df.groupby('roi').count()
dup_idx = dup[dup['slide_id'] == 2].index
dupped_df = mapped_df[mapped_df['roi'].isin(dup_idx)]

sample_count = mapped_df.sort_values(['core','panel']).groupby(['core', 'panel', 'GGO ID']).count().reset_index().pivot(index = ['core', 'GGO ID'], columns = 'panel', values = 'roi').fillna(0)
sample_count.to_csv('metadata/roi_count(missing).csv')
mapped_df.to_csv('metadata/roi_mapped_clinical_quant_exists.csv')

mapped_df['Smoking Status'] = pd.Categorical(mapped_df['Smoking Status'])
mapped_df['Smoking Status'] = mapped_df['Smoking Status'].replace({0:'Unknown'})

mapped_df['quant_file'] = 'processed/quantification/quantification_' + mapped_df['roi'] + '_full.h5ad'

panel_g = mapped_df[mapped_df['quant_file'].str.lower().str.contains('panel[ ]?_?g', regex = True)]['quant_file'].tolist()
panel_h = mapped_df[mapped_df['quant_file'].str.lower().str.contains('panel[ ]?_?h', regex = True)]['quant_file'].tolist()

import anndata
panel_g = anndata.concat([anndata.read(f) for f in panel_g])
panel_h = anndata.concat([anndata.read(f) for f in panel_h])


pgobs = panel_g.obs.merge(mapped_df[features_extract], how = 'left', left_on = 'roi', right_on = 'roi')
phobs = panel_h.obs.merge(mapped_df[features_extract], how = 'left',  left_on = 'roi', right_on = 'roi')

to_add = set(pgobs.columns) - set(panel_g.obs.columns)
for s in to_add:
	panel_g.obs[s] = pgobs[s].tolist()


to_add = set(phobs.columns) - set(panel_h.obs.columns)
for s in to_add:
	panel_h.obs[s] = phobs[s].tolist()

panel_g.write('results/panel_g.raw.labeled.h5ad')
panel_h.write('results/panel_h.raw.labeled.h5ad')

pg = phenotyping(panel_g, channels_exclude = exclude_channel, batch_variable = 'GGO ID', save_name = 'panel_g', clustering_method = 'parc')
#ph = phenotyping(panel_h, channels_exclude = exclude_channel, batch_variable = 'GGO ID', save_name = 'panel_h', clustering_method = 'parc')



# area1 = pd.read_csv(metadata['PANEL_G']['image_metadata']['ROI_AREA_FILE'], index_col = 0)
# area2 = pd.read_csv(metadata['PANEL_H']['image_metadata']['ROI_AREA_FILE'], index_col = 0)
# area = pd.concat([area1, area2], ignore_index = True).drop_duplicates()
