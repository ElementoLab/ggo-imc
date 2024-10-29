# ggo-imc
GGO IMC analysis repository to reproduce figures from manuscript

Spatial landscape of ground-glass opacity in early-stage lung adenocarcinoma

## Setup
```
# setting up conda environment
conda create -n ggo-imc python==3.9.10 -y
conda activate ggo-imc
pip install -r requirements.txt

# downloading data
make download
```

## Analysis Pipeline
For streamlined analysis to produce figures use the command
```
make run
```

Alternatively for individual plots:
```
# plot celltype heatmaps
make celltype

# to plot differential abundance
make diff_abundance

# plot interactions
make interaction

```