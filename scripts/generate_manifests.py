#!/usr/bin/env python3
"""Generate Phase 0b batch manifests from existing processed directories.

Scans processed/PANEL_G/ and processed/PANEL_H/ for sample directories
and writes per-panel batch TSVs to metadata/phase0/.

Usage (from project root):
    python scripts/generate_manifests.py
"""

from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(__file__).parent.parent


def main():
    processed_dir = PROJECT_DIR / "processed"
    phase0_dir = PROJECT_DIR / "metadata" / "phase0"
    phase0_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    for panel_dir in sorted(processed_dir.glob("PANEL_*")):
        if not panel_dir.is_dir():
            continue
        panel_name = panel_dir.name.lower()  # panel_g, panel_h
        rows = []

        for sample_dir in sorted(panel_dir.iterdir()):
            if not sample_dir.is_dir():
                continue
            rows.append({
                "sample_id": sample_dir.name,
                "processed_dir": str(sample_dir.relative_to(PROJECT_DIR)),
                "panel": panel_dir.name,
                "batch": panel_name,
            })

        if rows:
            df = pd.DataFrame(rows)
            out_path = phase0_dir / f"{panel_name}_samples.tsv"
            df.to_csv(out_path, sep="\t", index=False)
            print(f"Wrote {len(df)} samples to {out_path}")
            all_rows.extend(rows)

    if all_rows:
        all_df = pd.DataFrame(all_rows)
        all_path = phase0_dir / "all_samples.tsv"
        all_df.to_csv(all_path, sep="\t", index=False)
        print(f"Wrote {len(all_df)} total samples to {all_path}")


if __name__ == "__main__":
    main()
