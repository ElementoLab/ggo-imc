#!/usr/bin/env python
"""Integration benchmark: log1p+zscore+Harmony vs arcsinh+Harmony vs arcsinh+CytoVI.

For each panel (G and H), subsamples 50k cells (stratified by celltype),
runs three integration pipelines, and generates an HTML report with
scib-style metrics.

Usage:
    cd projects/imc/ggo-imc
    python scripts/benchmark_integration.py [--n-cells 50000] [--seed 42]
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata
import numpy as np
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent

# Files
PANEL_FILES = {
    "G": PROJECT_DIR / "results" / "panel_g.phenotyped.umap.h5ad",
    "H": PROJECT_DIR / "results" / "panel_h.phenotyped.umap.h5ad",
}

BATCH_KEY = "roi"
CELLTYPE_KEY = "celltype"
BATCH_VARIABLE_HARMONY = "GGO ID"  # original pipeline uses GGO ID for Harmony

# Channels to exclude (same as yaml_to_csv.py)
EXCLUDE_CHANNELS = [
    "80ArAr(ArAr80)",
    "131Xe(Xe131)",
    "138Ba(Ba138)",
    "190BCKG(BCKG190)",
    "<EMPTY>(Pt194)",
    "<EMPTY>(Pt195)",
    "<EMPTY>(Pt196)",
    "<EMPTY>(Pt198)",
    "208Pb(Pb208)",
]


def subsample_stratified(
    adata: anndata.AnnData,
    n_cells: int,
    celltype_key: str,
    seed: int = 42,
) -> anndata.AnnData:
    """Subsample preserving cell type proportions."""
    rng = np.random.default_rng(seed)
    obs = adata.obs

    if n_cells >= len(obs):
        logger.info("Requested %d cells but only %d available; using all", n_cells, len(obs))
        return adata.copy()

    proportions = obs[celltype_key].value_counts(normalize=True)
    indices = []
    for ct, frac in proportions.items():
        ct_idx = obs.index[obs[celltype_key] == ct].tolist()
        n_take = max(1, int(round(frac * n_cells)))
        n_take = min(n_take, len(ct_idx))
        chosen = rng.choice(ct_idx, size=n_take, replace=False)
        indices.extend(chosen.tolist())

    # Trim or pad to exactly n_cells
    if len(indices) > n_cells:
        indices = rng.choice(indices, size=n_cells, replace=False).tolist()

    return adata[indices].copy()


def get_raw_counts(adata: anndata.AnnData) -> anndata.AnnData:
    """Extract raw counts into a fresh AnnData, filtering excluded channels."""
    if adata.raw is not None:
        raw_adata = anndata.AnnData(
            X=adata.raw.X.copy() if hasattr(adata.raw.X, "copy") else adata.raw.X.toarray().copy(),
            obs=adata.obs.copy(),
            var=adata.raw.var.copy(),
        )
    else:
        raw_adata = adata.copy()

    # Filter excluded channels
    keep = ~raw_adata.var.index.isin(EXCLUDE_CHANNELS)
    raw_adata = raw_adata[:, keep].copy()

    # Remove DNA channels (store in obs like the original pipeline)
    dna_mask = raw_adata.var.index.str.contains(r"DNA\d")
    if dna_mask.any():
        raw_adata.obs["DNA"] = np.asarray(raw_adata[:, dna_mask].X.mean(axis=1)).flatten()
        raw_adata = raw_adata[:, ~dna_mask].copy()

    # Keep only channels with parentheses (marker channels)
    marker_mask = raw_adata.var.index.str.contains(r"\(")
    raw_adata = raw_adata[:, marker_mask].copy()

    return raw_adata


def run_method_a(adata: anndata.AnnData, batch_key: str) -> anndata.AnnData:
    """Method A: log1p + per-ROI z-score (cap 3) + global z-score + PCA + Harmony.

    Replicates the original yaml_to_csv.py phenotyping pipeline exactly.
    """
    a = adata.copy()
    a.raw = a.copy()

    # log1p
    sc.pp.log1p(a)

    # Per-ROI z-score with cap
    z_cap = 3.0
    ads = []
    for roi_name in a.obs[BATCH_KEY].unique():
        try:
            a2 = a[a.obs[BATCH_KEY] == roi_name, :].copy()
            sc.pp.scale(a2, max_value=z_cap)
            a2.X[a2.X < -z_cap] = -z_cap
            ads.append(a2)
        except ZeroDivisionError:
            logger.warning("ZeroDivisionError: skipping ROI %s", roi_name)
    a = anndata.concat(ads)

    # Global z-score
    sc.pp.scale(a)

    # PCA + Harmony
    sc.pp.pca(a)
    sc.external.pp.harmony_integrate(a, BATCH_VARIABLE_HARMONY, max_iter_harmony=20)

    a.obsm["X_method_a"] = a.obsm["X_pca_harmony"].copy()
    return a


def run_method_b(adata: anndata.AnnData, batch_key: str) -> anndata.AnnData:
    """Method B: arcsinh(cofactor=5) + scale + PCA + Harmony."""
    from sc_tools.pp import arcsinh_transform, pca, run_harmony, scale

    a = adata.copy()
    a.raw = a.copy()

    arcsinh_transform(a, cofactor=5)
    scale(a, max_value=10)
    pca(a)
    run_harmony(a, batch_key=BATCH_VARIABLE_HARMONY, max_iter_harmony=20)

    a.obsm["X_method_b"] = a.obsm["X_pca_harmony"].copy()
    return a


def run_method_c(adata: anndata.AnnData, batch_key: str) -> anndata.AnnData:
    """Method C: arcsinh(cofactor=5) + CytoVI."""
    from sc_tools.pp import arcsinh_transform, run_cytovi

    a = adata.copy()
    a.raw = a.copy()

    arcsinh_transform(a, cofactor=5)
    run_cytovi(a, batch_key=batch_key)

    a.obsm["X_method_c"] = a.obsm["X_cytovi"].copy()
    return a


def run_benchmark(n_cells: int = 50_000, seed: int = 42):
    """Run the integration benchmark for both panels."""
    from sc_tools.bm import compare_integrations
    from sc_tools.bm.report import generate_integration_report

    output_dir = PROJECT_DIR / "figures" / "benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    for panel, h5ad_path in PANEL_FILES.items():
        logger.info("=== Panel %s ===", panel)

        if not h5ad_path.exists():
            logger.warning("File not found: %s — skipping panel %s", h5ad_path, panel)
            continue

        # Load and subsample
        logger.info("Loading %s...", h5ad_path)
        adata_full = sc.read_h5ad(str(h5ad_path))
        logger.info("Full dataset: %d cells x %d features", adata_full.n_obs, adata_full.n_vars)

        if CELLTYPE_KEY not in adata_full.obs.columns:
            logger.error(
                "Column '%s' not found in obs. Available: %s",
                CELLTYPE_KEY,
                list(adata_full.obs.columns),
            )
            continue

        adata_sub = subsample_stratified(adata_full, n_cells, CELLTYPE_KEY, seed=seed)
        logger.info("Subsampled to %d cells", adata_sub.n_obs)

        # Get raw counts
        raw = get_raw_counts(adata_sub)
        logger.info("Raw marker matrix: %d cells x %d markers", raw.n_obs, raw.n_vars)

        # Ensure batch keys exist
        for key in [BATCH_KEY, BATCH_VARIABLE_HARMONY, CELLTYPE_KEY]:
            if key in adata_sub.obs.columns and key not in raw.obs.columns:
                raw.obs[key] = adata_sub.obs.loc[raw.obs.index, key].values

        # --- Method A: log1p + zscore + Harmony ---
        logger.info("Running Method A (log1p + zscore + Harmony)...")
        try:
            result_a = run_method_a(raw, BATCH_KEY)
            logger.info("Method A: embedding shape %s", result_a.obsm["X_method_a"].shape)
        except Exception:
            logger.exception("Method A failed")
            result_a = None

        # --- Method B: arcsinh + Harmony ---
        logger.info("Running Method B (arcsinh + Harmony)...")
        try:
            result_b = run_method_b(raw, BATCH_KEY)
            logger.info("Method B: embedding shape %s", result_b.obsm["X_method_b"].shape)
        except Exception:
            logger.exception("Method B failed")
            result_b = None

        # --- Method C: arcsinh + CytoVI ---
        logger.info("Running Method C (arcsinh + CytoVI)...")
        try:
            result_c = run_method_c(raw, BATCH_KEY)
            logger.info("Method C: embedding shape %s", result_c.obsm["X_method_c"].shape)
        except Exception:
            logger.exception("Method C failed")
            result_c = None

        # Merge embeddings into one AnnData for comparison
        # Use result_a as base (or whichever succeeded)
        base = result_a or result_b or result_c
        if base is None:
            logger.error("All methods failed for panel %s", panel)
            continue

        combined = base.copy()
        embeddings = {}
        if result_a is not None:
            combined.obsm["X_method_a"] = result_a.obsm["X_method_a"]
            embeddings["log1p+zscore+Harmony"] = "X_method_a"
        if result_b is not None:
            # Align indices in case concat changed ordering
            common_idx = combined.obs.index.intersection(result_b.obs.index)
            combined_loc = combined.obs.index.get_indexer(common_idx)
            b_loc = result_b.obs.index.get_indexer(common_idx)
            emb = np.zeros((combined.n_obs, result_b.obsm["X_method_b"].shape[1]))
            emb[combined_loc] = result_b.obsm["X_method_b"][b_loc]
            combined.obsm["X_method_b"] = emb
            embeddings["arcsinh+Harmony"] = "X_method_b"
        if result_c is not None:
            common_idx = combined.obs.index.intersection(result_c.obs.index)
            combined_loc = combined.obs.index.get_indexer(common_idx)
            c_loc = result_c.obs.index.get_indexer(common_idx)
            emb = np.zeros((combined.n_obs, result_c.obsm["X_method_c"].shape[1]))
            emb[combined_loc] = result_c.obsm["X_method_c"][c_loc]
            combined.obsm["X_method_c"] = emb
            embeddings["arcsinh+CytoVI"] = "X_method_c"

        if len(embeddings) < 2:
            logger.warning("Only %d method(s) succeeded for panel %s", len(embeddings), panel)

        # Run comparison
        logger.info("Computing integration metrics...")
        try:
            comparison_df = compare_integrations(
                combined,
                embeddings=embeddings,
                batch_key=BATCH_KEY,
                celltype_key=CELLTYPE_KEY,
                include_unintegrated=False,
            )
        except Exception:
            logger.exception("compare_integrations failed for panel %s", panel)
            continue

        # Save results
        csv_path = output_dir / f"integration_metrics_panel_{panel.lower()}.csv"
        comparison_df.to_csv(csv_path, index=False)
        logger.info("Metrics saved to %s", csv_path)

        # Generate HTML report
        try:
            report_path = generate_integration_report(
                comparison_df,
                adata=combined,
                embeddings=embeddings,
                batch_key=BATCH_KEY,
                celltype_key=CELLTYPE_KEY,
                output_path=output_dir / f"integration_benchmark_report_panel_{panel.lower()}.html",
                title=f"IMC Integration Benchmark - Panel {panel}",
            )
            logger.info("HTML report: %s", report_path)
        except Exception:
            logger.exception("Report generation failed (CSV is still saved)")

        # Print summary
        print(f"\n=== Panel {panel} Integration Benchmark ===")
        print(
            f"Cells: {combined.n_obs}, Batches ({BATCH_KEY}): {combined.obs[BATCH_KEY].nunique()}, "
            f"Cell types: {combined.obs[CELLTYPE_KEY].nunique()}"
        )
        print()
        for _, row in comparison_df.iterrows():
            print(f"  {row['method']}:")
            print(
                f"    Overall: {row.get('overall_score', 0):.3f}  "
                f"(batch: {row.get('batch_score', 0):.3f}, bio: {row.get('bio_score', 0):.3f})"
            )
            for metric in ["asw_batch", "pcr", "graph_connectivity", "asw_celltype", "ari", "nmi"]:
                if metric in row.index:
                    print(f"    {metric}: {row[metric]:.3f}")
            print()


def main():
    parser = argparse.ArgumentParser(description="Run IMC integration benchmark")
    parser.add_argument(
        "--n-cells", type=int, default=50000, help="Cells to subsample (default: 50000)"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    run_benchmark(n_cells=args.n_cells, seed=args.seed)


if __name__ == "__main__":
    main()
