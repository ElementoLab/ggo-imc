#!/usr/bin/env python
"""Segmentation benchmark: DeepCell vs Cellpose vs StarDist on 20 random ROIs.

Selects 10 ROIs from Panel G and 10 from Panel H (fixed seed), runs all three
segmentation methods, computes quality metrics via sc_tools.bm, and generates
an HTML report.

Input: *_Probabilities.tiff files (H, W, 3) from Ilastik pixel classification.
Channels: 0=background, 1=nucleus, 2=cytoplasm.

Strategy: Run StarDist first (TF), then Cellpose (torch) to avoid torch/TF
initialization conflicts on macOS ARM.

Usage:
    cd projects/imc/ggo-imc
    python scripts/benchmark_segmentation.py [--n-per-panel 10] [--seed 42] [--gpu]
"""

from __future__ import annotations

import argparse
import gc
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent


def discover_rois(panel: str) -> list[dict]:
    """Find all ROIs with both Probabilities TIFF and DeepCell mask."""
    panel_dir = PROJECT_DIR / "processed" / f"PANEL_{panel}"
    if not panel_dir.exists():
        logger.warning("Panel directory not found: %s", panel_dir)
        return []

    rois = []
    for sample_dir in sorted(panel_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        tiff_dir = sample_dir / "tiffs"
        if not tiff_dir.exists():
            continue
        for prob_path in sorted(tiff_dir.glob("*_Probabilities.tiff")):
            # Derive mask path: *_Probabilities.tiff -> *_full_mask.tiff
            stem = prob_path.stem.replace("_Probabilities", "")
            mask_path = prob_path.parent / f"{stem}_full_mask.tiff"
            if mask_path.exists():
                rois.append(
                    {
                        "panel": panel,
                        "sample": sample_dir.name,
                        "roi_name": stem,
                        "prob_path": str(prob_path),
                        "mask_path": str(mask_path),
                    }
                )
    return rois


def run_benchmark(n_per_panel: int = 10, seed: int = 42, use_gpu: bool = False):
    """Run the segmentation benchmark."""
    rng = np.random.default_rng(seed)

    # Discover and sample ROIs
    selected_rois = []
    for panel in ["G", "H"]:
        rois = discover_rois(panel)
        logger.info("Panel %s: found %d ROIs with Probabilities + mask", panel, len(rois))
        if len(rois) == 0:
            continue
        n_select = min(n_per_panel, len(rois))
        selected_idx = rng.choice(len(rois), size=n_select, replace=False)
        for i in sorted(selected_idx):
            selected_rois.append(rois[i])

    logger.info("Selected %d ROIs total", len(selected_rois))

    # Store ROI metadata (load images lazily to save memory)
    roi_names = [r["roi_name"] for r in selected_rois]
    roi_meta = {r["roi_name"]: r for r in selected_rois}

    # --- Pass 1: Run StarDist on all ROIs (TF only, no torch loaded) ---
    logger.info("=== Pass 1: StarDist (TensorFlow) ===")
    stardist_masks = {}
    try:
        from sc_tools.bm.segment import run_stardist

        for name in roi_names:
            logger.info("  StarDist on %s...", name)
            t0 = time.time()
            try:
                prob_img = tifffile.imread(roi_meta[name]["prob_path"])
                mask = run_stardist(prob_img, nuclear_idx=1)
                stardist_masks[name] = mask
                n_cells = len(np.unique(mask)) - 1
                logger.info("  StarDist %s: %d cells (%.1fs)", name, n_cells, time.time() - t0)
                del prob_img
            except Exception:
                logger.exception("  StarDist failed on %s", name)
    except ImportError:
        logger.warning("StarDist not available (needs tensorflow). Skipping.")

    gc.collect()

    # --- Pass 2: Run Cellpose on all ROIs (torch) ---
    logger.info("=== Pass 2: Cellpose (PyTorch) ===")
    cellpose_masks = {}
    try:
        from sc_tools.bm.segment import run_cellpose

        for name in roi_names:
            logger.info("  Cellpose on %s...", name)
            t0 = time.time()
            try:
                prob_img = tifffile.imread(roi_meta[name]["prob_path"])
                mask = run_cellpose(prob_img, nuclear_idx=1, cytoplasm_idx=2, gpu=use_gpu)
                cellpose_masks[name] = mask
                n_cells = len(np.unique(mask)) - 1
                logger.info("  Cellpose %s: %d cells (%.1fs)", name, n_cells, time.time() - t0)
                del prob_img
            except Exception:
                logger.exception("  Cellpose failed on %s", name)
    except ImportError:
        logger.warning("Cellpose not available. Skipping.")

    # --- Compare segmentations per ROI ---
    logger.info("=== Computing metrics ===")
    from sc_tools.bm.segmentation import compare_segmentations
    from sc_tools.bm.report import generate_segmentation_report

    all_results = []
    for name in roi_names:
        meta = roi_meta[name]
        deepcell_mask = tifffile.imread(meta["mask_path"])
        if deepcell_mask.ndim > 2:
            deepcell_mask = deepcell_mask[0]

        # Nuclear probability channel for marker quality metric
        prob_img = tifffile.imread(meta["prob_path"])
        nuclear_img = prob_img[:, :, 1].astype(np.float32)

        masks = {"DeepCell": deepcell_mask.astype(np.int32)}
        if name in cellpose_masks:
            # Cellpose mask may be at prob resolution (2x); resize to match DeepCell
            cp_mask = cellpose_masks[name]
            if cp_mask.shape != deepcell_mask.shape:
                from skimage.transform import resize

                cp_mask = resize(
                    cp_mask,
                    deepcell_mask.shape,
                    order=0,
                    preserve_range=True,
                    anti_aliasing=False,
                ).astype(np.int32)
            masks["Cellpose"] = cp_mask.astype(np.int32)
        if name in stardist_masks:
            sd_mask = stardist_masks[name]
            if sd_mask.shape != deepcell_mask.shape:
                from skimage.transform import resize

                sd_mask = resize(
                    sd_mask,
                    deepcell_mask.shape,
                    order=0,
                    preserve_range=True,
                    anti_aliasing=False,
                ).astype(np.int32)
            masks["StarDist"] = sd_mask.astype(np.int32)

        if len(masks) < 2:
            logger.warning("Only DeepCell for %s, skipping", name)
            continue

        # Resize nuclear image to mask resolution for marker quality
        if nuclear_img.shape != deepcell_mask.shape:
            from skimage.transform import resize

            nuclear_img = resize(nuclear_img, deepcell_mask.shape, preserve_range=True).astype(
                np.float32
            )

        try:
            df = compare_segmentations(masks, intensity_image=nuclear_img)
            df["roi"] = name
            df["panel"] = meta["panel"]
            all_results.append(df)
            logger.info(
                "  %s: %s",
                name,
                ", ".join(f"{r['method']}={r['n_cells']}" for _, r in df.iterrows()),
            )
        except Exception:
            logger.exception("compare_segmentations failed on %s", name)

        del prob_img, nuclear_img

    if not all_results:
        logger.error("No results collected. Check data paths and dependencies.")
        sys.exit(1)

    # Aggregate results
    results_df = pd.concat(all_results, ignore_index=True)

    numeric_cols = results_df.select_dtypes(include=[np.number]).columns.tolist()
    exclude = {"composite_score"}
    agg_cols = [c for c in numeric_cols if c not in exclude]

    summary_rows = []
    for method in results_df["method"].unique():
        method_df = results_df[results_df["method"] == method]
        row = {"method": method, "n_rois": len(method_df)}
        for col in agg_cols:
            row[f"{col}_mean"] = method_df[col].mean()
            row[f"{col}_std"] = method_df[col].std()
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)

    report_df = pd.DataFrame(summary_rows)
    report_df = report_df.rename(columns={f"{c}_mean": c for c in agg_cols})
    if "composite_score" not in report_df.columns:
        report_df["composite_score"] = 50.0
    if "median_circularity" in report_df.columns:
        report_df = report_df.sort_values("median_circularity", ascending=False).reset_index(
            drop=True
        )

    # Save results
    output_dir = PROJECT_DIR / "figures" / "benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    results_df.to_csv(output_dir / "segmentation_per_roi.csv", index=False)
    summary_df.to_csv(output_dir / "segmentation_summary.csv", index=False)
    logger.info("Per-ROI results: %s", output_dir / "segmentation_per_roi.csv")
    logger.info("Summary: %s", output_dir / "segmentation_summary.csv")

    # Generate HTML report
    try:
        report_path = generate_segmentation_report(
            comparison_df=report_df,
            output_path=output_dir / "segmentation_benchmark_report.html",
            title="IMC Segmentation Benchmark: DeepCell vs Cellpose vs StarDist",
        )
        logger.info("HTML report: %s", report_path)
    except Exception:
        logger.exception("Report generation failed (CSVs are still saved)")

    # Print summary
    print("\n=== Segmentation Benchmark Summary ===")
    print(f"ROIs evaluated: {len(results_df['roi'].unique())}")
    print()
    for _, row in summary_df.iterrows():
        print(f"  {row['method']}:")
        print(
            f"    Cells (mean +/- std): {row.get('n_cells_mean', 0):.0f}"
            f" +/- {row.get('n_cells_std', 0):.0f}"
        )
        print(
            f"    Circularity: {row.get('median_circularity_mean', 0):.3f}"
            f" +/- {row.get('median_circularity_std', 0):.3f}"
        )
        print(
            f"    Solidity: {row.get('median_solidity_mean', 0):.3f}"
            f" +/- {row.get('median_solidity_std', 0):.3f}"
        )
        if "mean_snr_mean" in row.index:
            print(
                f"    SNR: {row.get('mean_snr_mean', 0):.3f} +/- {row.get('mean_snr_std', 0):.3f}"
            )
        print()


def main():
    parser = argparse.ArgumentParser(description="Run IMC segmentation benchmark")
    parser.add_argument("--n-per-panel", type=int, default=10, help="ROIs per panel (default: 10)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for Cellpose")
    args = parser.parse_args()

    run_benchmark(n_per_panel=args.n_per_panel, seed=args.seed, use_gpu=args.gpu)


if __name__ == "__main__":
    main()
