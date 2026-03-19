# External Integrations

**Analysis Date:** 2026-03-19

## APIs & External Services

**Data Download:**
- Box (Weill Cornell Medicine shared storage)
  - Service: Cloud file storage and sharing
  - What it's used for: Distribution of metadata YAML config file
  - SDK/Client: `requests` (HTTP GET)
  - URL: `https://wcm.box.com/shared/static/mdntp2xf9tjobxeidkw93mg8jysb7nh9.yml`
  - Invoked by: `scripts/download_yaml.py` via `utils.load_config()`

**Data Repository:**
- Zenodo
  - Service: Open-access scientific data repository
  - What it's used for: Long-term archival of processed and raw IMC data
  - DOI references documented in README.md but accessed via external tools (wget/curl, not within Python)
  - Processed Data DOI: 10.5281/zenodo.14822106
  - Raw Data DOI: 10.5281/zenodo.14822106

## Data Storage

**Databases:**
- Not applicable — no SQL/NoSQL databases used
- Data format: HDF5 (via AnnData .h5ad files)
  - Location: Referenced in metadata YAML under `AnnData` paths
  - Client: `scanpy.read()` (wraps h5py/tables)
  - Backup URLs: Metadata includes `backup_url` for each panel in case primary path is unavailable

**File Storage:**
- Local filesystem only
  - Data: AnnData matrices (.h5ad files) stored locally after download
  - Processed outputs: Figures stored in `figures/figure{1-5}/` subdirectories
  - Artifacts: Sentinel files (*.done) track Snakemake rule completion
  - Sample outputs: `.pdf`, `.png`, `.svg` formats for publication

**Caching:**
- Not detected — no explicit caching layer
- Snakemake uses sentinel files (*.done) to avoid re-running completed rules
- Backup URLs in metadata YAML serve as failover mechanism for data loading

## Authentication & Identity

**Auth Provider:**
- Not applicable — no user authentication system
- Data access via Box: Public shared link (no credentials required in code)
- Zenodo access: Direct download (no authentication); restricted access may be enforced by repository until publication

## Monitoring & Observability

**Error Tracking:**
- Not detected — no error reporting service (Sentry, etc.)
- Errors surfaced via stdout/stderr from Python scripts

**Logs:**
- Approach: Standard Python logging (via scanpy.logging)
- scanpy.logging outputs to stdout/stderr with configurable verbosity
- Snakemake logs rule execution to .log files (not detected in this project)

## CI/CD & Deployment

**Hosting:**
- GitHub repository (source of truth for code)
- Zenodo (permanent archival of data)
- No continuous deployment — figures are generated locally by researchers or on HPC

**CI Pipeline:**
- Pytest-based smoke testing
  - Config: `pyproject.toml` [tool.pytest.ini_options]
  - Test files: `tests/test_imports.py`, `tests/test_units.py`
  - Fixtures: `tests/conftest.py` with mocked imc_analysis, scanpy, and seaborn
  - Purpose: Verify all 10 pipeline scripts import without errors and execute basic workflows
- No automated GitHub Actions detected

**Snakemake Orchestration:**
- Workflow definition: `Snakefile`
- Invocation: `snakemake -d . -s Snakefile [target]`
- Targets: `figure1` through `figure5`, `all`, `download`, `clean`

## Environment Configuration

**Required env vars:**
- Not detected — no environment variables required for normal operation
- Optional: `SC_TOOLS_RUNTIME` (controls container mode; default: auto-detect via run_container.sh)

**Secrets location:**
- Not applicable — no API keys, tokens, or credentials in code
- Box URL is public shared link
- Zenodo access is direct download (no authentication)

## Webhooks & Callbacks

**Incoming:**
- Not detected — no webhook endpoints

**Outgoing:**
- Not detected — no callbacks to external services

## Data Sources & Workflows

**Download Workflow:**
```
1. User runs: snakemake figure1 (or any rule)
2. Rule invokes: python scripts/download_yaml.py
3. download_yaml.py calls: requests.get(CONFIG_URL)
4. Config downloaded to: metadata/ggo_config.yml
5. All scripts call: utils.load_config() → reads local YAML
6. Metadata contains AnnData paths (local .h5ad files)
7. Scripts load panels: utils.load_panels(metadata) → sc.read(path, backup_url=url)
```

**Backup URL Pattern:**
- Each panel in metadata YAML includes a `backup_url` field
- If primary path does not exist on filesystem, scanpy.read() falls back to backup_url
- Backup URLs are HTTP(S) endpoints serving .h5ad files (via Box or cloud storage)

**Container Execution:**
```
1. Snakemake rule: shell: run_container("scripts/celltype_heatmap_info.py")
2. run_container() function wraps: ./scripts/run_container.sh {project} python {script}
3. run_container.sh auto-detects runtime (Docker on macOS, Apptainer on Linux)
4. Spins up sc_tools container with mounted project directory
5. Python script executes inside container with isolated dependencies
6. Output files written to mounted project directory (visible on host)
```

---

*Integration audit: 2026-03-19*
