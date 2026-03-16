#!/usr/bin/env python3
"""
Small Molecule Drug Design Pipeline
====================================
AI-guided generation, filtering, scoring, and iterative optimization
of small molecule candidates targeting the IL-6 Site II binding pocket.

Usage:
    python3 smallmol_pipeline.py                     # full run
    python3 smallmol_pipeline.py --max-molecules 100 # quick demo
    python3 smallmol_pipeline.py --skip-phase 3      # resume from phase 4
    python3 smallmol_pipeline.py --no-figures         # skip visualization
"""

import argparse
import logging
import os
import sys
import time

import pandas as pd

import smallmol_config as cfg


def banner():
    lines = [
        "",
        "=" * 65,
        "  SMALL MOLECULE DRUG DESIGN PIPELINE",
        "  Target: IL-6 Site II | Modality: Small Molecules",
        "  Method: Fragment-based generation + MC optimization",
        "=" * 65,
        "",
    ]
    return "\n".join(lines)


def setup_logging():
    os.makedirs(cfg.RESULTS_DIR, exist_ok=True)
    log_path = os.path.join(cfg.RESULTS_DIR, "smallmol_pipeline.log")

    root = logging.getLogger()
    root.setLevel(getattr(logging, cfg.LOG_LEVEL))

    fmt = logging.Formatter("%(asctime)s  %(name)-30s  %(levelname)-7s  %(message)s",
                            datefmt="%H:%M:%S")

    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(fmt)
    root.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    root.addHandler(sh)

    return log_path


def parse_args():
    parser = argparse.ArgumentParser(description="Small Molecule Drug Design Pipeline")
    parser.add_argument("--max-molecules", type=int, default=None,
                        help="Override number of molecules to generate")
    parser.add_argument("--skip-phase", type=int, default=0,
                        help="Skip phases up to this number (resume from next)")
    parser.add_argument("--no-figures", action="store_true",
                        help="Skip figure generation")
    return parser.parse_args()


def main():
    args = parse_args()
    log_path = setup_logging()
    logger = logging.getLogger("smallmol_pipeline")

    print(banner())
    logger.info("Small Molecule Pipeline starting")
    logger.info("Results directory: %s", cfg.RESULTS_DIR)

    if args.max_molecules:
        cfg.N_MOLECULES_GENERATE = args.max_molecules
        logger.info("Override: generating %d molecules", args.max_molecules)

    skip = args.skip_phase
    t0 = time.time()

    os.makedirs(cfg.RESULTS_DIR, exist_ok=True)
    os.makedirs(cfg.FIGURES_DIR, exist_ok=True)

    # Lazy imports to avoid loading RDKit until needed
    from src.smallmol.phase1_pocket import run_phase1
    from src.smallmol.phase2_generate import run_phase2
    from src.smallmol.phase3_filter import run_phase3
    from src.smallmol.phase4_score import run_phase4
    from src.smallmol.phase5_optimize import run_phase5
    from src.smallmol.phase6_analyze import run_phase6
    from src.smallmol.visualize import generate_all_figures

    # CSV paths
    csv_gen = os.path.join(cfg.RESULTS_DIR, "01_molecules_generated.csv")
    csv_filt = os.path.join(cfg.RESULTS_DIR, "02_molecules_filtered.csv")
    csv_scored = os.path.join(cfg.RESULTS_DIR, "03_molecules_scored.csv")
    csv_history = os.path.join(cfg.RESULTS_DIR, "04_optimization_history.csv")
    csv_variants = os.path.join(cfg.RESULTS_DIR, "05_optimized_variants.csv")
    csv_leads = os.path.join(cfg.RESULTS_DIR, "06_lead_analysis.csv")
    pocket_json = os.path.join(cfg.RESULTS_DIR, "pocket_profile.json")

    # ── Phase 1: Pocket Analysis ─────────────────────────────────────────────
    if skip >= 1:
        import json
        with open(pocket_json, "r") as f:
            pocket_profile = json.load(f)
        logger.info("Phase 1 skipped — loaded pocket profile from cache")
    else:
        pocket_profile = run_phase1(cfg.INTERFACE_JSON, output_json=pocket_json)

    # ── Phase 2: Generation ──────────────────────────────────────────────────
    if skip >= 2 and os.path.exists(csv_gen):
        df_gen = pd.read_csv(csv_gen)
        logger.info("Phase 2 skipped — loaded %d molecules from cache", len(df_gen))
    else:
        df_gen = run_phase2(pocket_profile, output_csv=csv_gen)

    # ── Phase 3: Filtering ───────────────────────────────────────────────────
    if skip >= 3 and os.path.exists(csv_filt):
        df_filt = pd.read_csv(csv_filt)
        logger.info("Phase 3 skipped — loaded %d molecules from cache", len(df_filt))
    else:
        df_filt = run_phase3(df_gen, output_csv=csv_filt)

    # Store full df (with pass/fail flags) for figures
    df_all = df_gen.copy()
    if "passes_all" not in df_all.columns:
        # Re-apply filters just for the flags
        from src.smallmol.phase3_filter import apply_filters
        df_all = apply_filters(df_all)

    # ── Phase 4: Scoring ─────────────────────────────────────────────────────
    if skip >= 4 and os.path.exists(csv_scored):
        df_scored = pd.read_csv(csv_scored)
        logger.info("Phase 4 skipped — loaded %d scored molecules from cache", len(df_scored))
    else:
        df_scored = run_phase4(df_filt, pocket_profile, output_csv=csv_scored)

    # ── Phase 5: Optimization ────────────────────────────────────────────────
    if skip >= 5 and os.path.exists(csv_history):
        history_df = pd.read_csv(csv_history)
        variants_df = pd.read_csv(csv_variants) if os.path.exists(csv_variants) else pd.DataFrame()
        logger.info("Phase 5 skipped — loaded optimization history from cache")
    else:
        history_df, variants_df = run_phase5(
            df_scored, pocket_profile,
            output_history_csv=csv_history,
            output_variants_csv=csv_variants,
        )

    # ── Phase 6: Lead Analysis ───────────────────────────────────────────────
    if skip >= 6 and os.path.exists(csv_leads):
        lead_df = pd.read_csv(csv_leads)
        logger.info("Phase 6 skipped — loaded lead analysis from cache")
    else:
        lead_df = run_phase6(df_scored, variants_df, pocket_profile, output_csv=csv_leads)

    # ── Figures ──────────────────────────────────────────────────────────────
    if not args.no_figures:
        generate_all_figures(df_all, df_scored, history_df, lead_df, cfg.FIGURES_DIR)

    # ── Summary ──────────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    logger.info("")
    logger.info("=" * 65)
    logger.info("  PIPELINE COMPLETE")
    logger.info("=" * 65)
    logger.info("  Molecules generated:  %d", len(df_gen))
    logger.info("  Passed filters:       %d", len(df_filt))
    logger.info("  Scored:               %d", len(df_scored))
    logger.info("  Leads analyzed:       %d", len(lead_df))
    logger.info("  Elapsed time:         %.1f seconds", elapsed)
    logger.info("  Results:              %s", cfg.RESULTS_DIR)
    logger.info("  Log:                  %s", log_path)
    logger.info("=" * 65)


if __name__ == "__main__":
    main()
