# External Integrations

**Analysis Date:** 2026-03-19

## APIs & External Services

**Box (Weill Cornell Shared Storage):**
- Configuration file download: `https://wcm.box.com/shared/static/mdntp2xf9tjobxeidkw93mg8jysb7nh9.yml`
- Panel G backup datasets (multiple h5ad files via URLs):
  - Phenotyped UMAP: `https://wcm.box.com/shared/static/u7po2u5f3stjxl27rw7qxble4j7kur1g.h5ad`
  - Lymphocytes: `https://wcm.box.com/shared/static/qxxxd62zzmme42uhks5dy5jsxnctb0al.h5ad`
  - Myeloids: `https://wcm.box.com/shared/static/0cplylp3n6t4udzvuqtnnd31q3rhyq3y.h5ad`
  - UTAG: `https://wcm.box.com/shared/static/0tj5jcl4m37muhm9afnvks4b3sd9lpsz.h5ad`
- Panel H backup datasets (multiple h5ad files via URLs):
  - Phenotyped UMAP: `https://wcm.box.com/shared/static/ueuwm6wufc479sbhd5utc6p6vkgiated.h5ad`
  - Lymphocytes: `https://wcm.box.com/shared/static/puxd3c4eevitsk55x7nhbhdie81muygf.h5ad`
  - Myeloids: `https://wcm.box.com/shared/static/qdj6q7rcgumfucil6dzhxcw06lsk7x20.h5ad`
- PCA anndata: `https://wcm.box.com/shared/static/faeuq6ojwqzm2wb7mlsv7hgi1jmg5yst.h5ad`
- Patient celltype grouped: `https://wcm.box.com/shared/static/ogypsajbwwlrpjuifpri4gng1rf1rc2z.h5ad`
  - SDK/Client: `requests` library
  - Fallback mechanism: `sc.read(..., backup_url=...)` in scanpy handles automatic download if local file missing

## Data Storage

**Databases:**
- None - Analysis uses file-based HDF5 (AnnData) format stored locally

**File Storage:**
- Local filesystem only
  - Local data directory: `data/` (processed IMC images)
  - Results: `results/` (HDF5 output)
  - Processed images: `processed/PANEL_G/` and `processed/PANEL_H/` (TIFF images and masks)
  - Metadata: `metadata/` (CSV ROI data, Excel clinical annotations, downloadable YAML config)

**Remote Backup:**
- Box.com for h5ad backup files
- Automatic fallback: if local file missing, scanpy's `backup_url` parameter triggers HTTP download

**Caching:**
- None - No explicit caching layer; Snakemake manages intermediate outputs via touch files

## Authentication & Identity

**Auth Provider:**
- None - Box URLs are public/shared links (no credentials required)
- Configuration stored in code: URLs hardcoded in `metadata/ggo_config.yml`

## Monitoring & Observability

**Error Tracking:**
- None detected

**Logs:**
- Console output via Snakemake and scanpy logging
- Snakemake log: `logs/` directory (managed by workflow)
- scanpy logging: `sc.logging.info(...)` calls in pipeline scripts (e.g., `scripts/utils.py`)

## CI/CD & Deployment

**Hosting:**
- Local execution (Weill Cornell institutional environment)
- No cloud deployment detected

**CI Pipeline:**
- None - Manual execution via Snakemake
- Containerization: Docker (macOS) / Apptainer (Linux) for reproducibility
- Container built externally as `containers/sc_tools.sif`

## Environment Configuration

**Required env vars:**
- `SC_TOOLS_RUNTIME` - Controls container runtime:
  - Default: auto-detect (Docker on macOS, Apptainer on Linux)
  - `SC_TOOLS_RUNTIME=none` - Native conda environment without container

**Configuration files:**
- `config.yaml` - Snakemake config (repo paths, container image)
- `metadata/ggo_config.yml` - Panel-specific paths, color schemes, clinical metadata (downloaded at runtime)

**Secrets location:**
- Not applicable - No API keys or credentials required
- Box URLs are public shared links

## Webhooks & Callbacks

**Incoming:**
- None

**Outgoing:**
- None

## Data Download Flow

**Configuration Download (`scripts/download_yaml.py`):**
1. HTTP GET to Box URL via `requests.get()`
2. Stream response with `tqdm` progress bar
3. Save to `metadata/ggo_config.yml`
4. Called by Snakemake `download` rule on startup

**AnnData Download (implicit in `sc.read(..., backup_url=...)`):**
1. Scanpy checks for local file first
2. If missing, HTTP GET to Box backup URL
3. Save and load as HDF5

**Image Data:**
- Pre-processed TIFF files stored locally in `processed/PANEL_G/*/tiffs/` and `processed/PANEL_H/*/tiffs/`
- No external download mechanism (loaded from local filesystem)

## Clinical & Metadata Integration

**Clinical Data:**
- Source: `metadata/De-identified J&J Clinical Annotations_UpdatedAAH031523 (1).xlsx`
- Format: Excel spreadsheet with de-identified patient metadata
- Integration: Loaded in patient analysis scripts (e.g., `scripts/patient_group.py`)
- Fields: pathology, radiology, smoking status, gender, race, risk predictions

**ROI Metadata:**
- CSV files for each panel (Panel_G_ROI_area.csv, Panel_H_ROI_mean.csv, etc.)
- Contains: ROI areas, mean intensity, variance, nucleus/cytoplasm masks, spillover matrices
- Used for spatial analysis and quality control

---

*Integration audit: 2026-03-19*
