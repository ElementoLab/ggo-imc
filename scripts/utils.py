"""
utils.py — Shared utilities for the ggo-imc pipeline scripts.

All pipeline scripts (figure1–figure5) import from here to avoid
repeating config loading, AnnData loading, and figure-saving boilerplate.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scanpy as sc
import imc_analysis as imc

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CONFIG_PATH = 'metadata/ggo_config.yml'
CONFIG_URL = 'https://wcm.box.com/shared/static/mdntp2xf9tjobxeidkw93mg8jysb7nh9.yml'

PANELS = ['PANEL_G', 'PANEL_H']

CYTOKINE: list[str] = ['IFNg', 'IL1alpha', 'IL1beta', 'IL1R1', 'IL12p40', 'IL17A', 'IL23p19', 'IL23R']


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(
    config_path: str = CONFIG_PATH,
    *,
    download_url: str | None = CONFIG_URL,
) -> dict:
    """Load the project YAML config, optionally downloading it first.

    Parameters
    ----------
    config_path:
        Path to ``ggo_config.yml`` (relative to the working directory).
    download_url:
        If provided and the file does not already exist, download from this
        URL before parsing.  Pass ``None`` to skip the download step.

    Returns
    -------
    dict
        Parsed YAML content.
    """
    if download_url is not None:
        imc.utils.download(download_url, config_path)
    return imc.utils.parse_yaml(config_path)


# ---------------------------------------------------------------------------
# AnnData loading
# ---------------------------------------------------------------------------

def load_panels(
    metadata: dict,
    panels: list[str] = PANELS,
    *,
    anndata_key: str = 'phenotyped_umap_name',
    backup_key: str = 'backup_url',
) -> dict:
    """Load per-panel AnnData objects from the config.

    Reads ``metadata[panel]['AnnData'][anndata_key]`` for each panel and
    falls back to ``metadata[panel]['AnnData'][backup_key]`` if the local
    file is absent.

    Parameters
    ----------
    metadata:
        Parsed config dict (from :func:`load_config`).
    panels:
        Which panels to load, e.g. ``['PANEL_G', 'PANEL_H']``.
    anndata_key:
        Key inside ``metadata[panel]['AnnData']`` for the local file path.
    backup_key:
        Key inside ``metadata[panel]['AnnData']`` for the remote backup URL.

    Returns
    -------
    dict
        Mapping ``panel -> AnnData``.

    Raises
    ------
    ValueError
        If a panel name or required key is missing from *metadata*.
    """
    for p in panels:
        if p not in metadata:
            raise ValueError(
                f"Panel '{p}' not found in config. "
                f"Available keys: {sorted(metadata.keys())}"
            )
        anndata_section = metadata[p].get('AnnData')
        if anndata_section is None:
            raise ValueError(
                f"Config section metadata['{p}']['AnnData'] is missing."
            )
        for key in (anndata_key, backup_key):
            if key not in anndata_section:
                raise ValueError(
                    f"Key '{key}' not found in metadata['{p}']['AnnData']. "
                    f"Available keys: {sorted(anndata_section.keys())}"
                )

    adata_dict: dict = {}
    for p in panels:
        path = metadata[p]['AnnData'][anndata_key]
        url = metadata[p]['AnnData'][backup_key]
        sc.logging.info(f"Reading {path}...")
        adata_dict[p] = sc.read(path, backup_url=url)
    return adata_dict


def load_single_panel(
    metadata: dict,
    panel: str,
    *,
    anndata_key: str = 'phenotyped_umap_name',
    backup_key: str = 'backup_url',
):
    """Load one panel's AnnData object.

    Convenience wrapper around :func:`load_panels` for scripts that only
    need a single panel.

    Parameters
    ----------
    metadata:
        Parsed config dict.
    panel:
        Panel name, e.g. ``'PANEL_H'``.
    anndata_key, backup_key:
        Same as in :func:`load_panels`.

    Returns
    -------
    AnnData
    """
    return load_panels(metadata, [panel], anndata_key=anndata_key, backup_key=backup_key)[panel]


# ---------------------------------------------------------------------------
# Figure saving
# ---------------------------------------------------------------------------

def save_figure(
    path: str,
    *,
    fig: plt.Figure | None = None,
    bbox_inches: str = 'tight',
    close: bool = True,
    **savefig_kwargs,
) -> None:
    """Save the current (or a given) figure, creating parent directories first.

    Parameters
    ----------
    path:
        Destination file path, e.g. ``'figures/figure1/foo.pdf'``.
    fig:
        Matplotlib figure to save.  If ``None``, uses ``plt.gcf()``.
    bbox_inches:
        Passed to ``savefig``.  Defaults to ``'tight'``.
    close:
        If ``True`` (default), call ``plt.close()`` after saving.
    **savefig_kwargs:
        Any additional keyword arguments forwarded to ``savefig``.
    """
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    target = fig if fig is not None else plt.gcf()
    target.savefig(path, bbox_inches=bbox_inches, **savefig_kwargs)
    if close:
        plt.close(target)


def ensure_dir(path: str) -> str:
    """Create *path* (and any missing parents) and return it unchanged.

    Thin wrapper around ``os.makedirs`` kept here so scripts don't need to
    import ``os`` just for directory creation.
    """
    os.makedirs(path, exist_ok=True)
    return path
