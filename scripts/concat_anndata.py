from glob import glob
import anndata

outdir = 'results/'

patterns = {
	"panel_G": ["PanelG", "Panel G"],
	"panel_H": ["PanelH", "Panel H"]
}

# patterns = {
# 	"panel1": ["[!Panel2]", "[!panel2]"],
# 	"panel2": ["Panel2", "panel2"]
# }

adata_dict = dict()
for ps in patterns:
	files = []
	for p in patterns[ps]:
		files = files + glob(f'processed/quantification/*{p}*.h5ad')

	adatas = [anndata.read(f) for f in files]
	adata_dict[ps] = anndata.concat(adatas)
	adata_dict[ps].write(f"{outdir}/{ps}.raw.h5ad")

