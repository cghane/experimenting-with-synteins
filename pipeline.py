#!/usr/bin/env python3
"""
AI-Guided Mini-Protein Binder Design Pipeline
Target: IL-6 Site II | Scaffold: Affibody Z-domain | Lab: Abiologics

Usage:
    python pipeline.py                     # Full run, mock structure prediction
    python pipeline.py --real              # Real ESMFold API calls
    python pipeline.py --max-candidates 50 # Limit for quick test
    python pipeline.py --skip-phase 4      # Resume from Phase 5

All outputs saved to results/ directory.
"""

import argparse
import logging
import os
import sys
import time

import pandas as pd

# ── Project imports ───────────────────────────────────────────────────────────
from config import (
    RECEPTOR_PDB, INTERFACE_JSON, RESULTS_DIR, TARGET_PDB_ID,
    RECEPTOR_CHAIN, LOG_LEVEL
)

# ── Logging setup ─────────────────────────────────────────────────────────────
os.makedirs(RESULTS_DIR, exist_ok=True)
logging.basicConfig(
    level=getattr(logging, LOG_LEVEL, logging.INFO),
    format='%(asctime)s  %(levelname)-7s  %(message)s',
    datefmt='%H:%M:%S',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(RESULTS_DIR, 'pipeline.log'), mode='w'),
    ]
)
logger = logging.getLogger(__name__)

# ── Result file paths ─────────────────────────────────────────────────────────
P2_CSV  = os.path.join(RESULTS_DIR, "01_candidates_generated.csv")
P3_CSV  = os.path.join(RESULTS_DIR, "02_candidates_filtered.csv")
P4_CSV  = os.path.join(RESULTS_DIR, "03_structure_metrics.csv")
P5_CSV  = os.path.join(RESULTS_DIR, "04_candidates_ranked.csv")
P6H_CSV = os.path.join(RESULTS_DIR, "05_optimization_history.csv")
P6V_CSV = os.path.join(RESULTS_DIR, "06_optimized_variants.csv")
P7_CSV  = os.path.join(RESULTS_DIR, "07_synteins_analysis.csv")


def banner():
    print("""
╔═══════════════════════════════════════════════════════════════════╗
║      AI-Guided Mini-Protein Binder Design Pipeline               ║
║      Target : IL-6 Site II (gp130 binding interface)             ║
║      Scaffold: Affibody Z-domain (55 aa, 3-helix bundle)         ║
║      Lab    : Abiologics                                          ║
╚═══════════════════════════════════════════════════════════════════╝
""")


def parse_args():
    p = argparse.ArgumentParser(description="Abiologics binder design pipeline")
    p.add_argument('--real', action='store_true',
                   help='Use ESMFold API for real structure prediction (slower)')
    p.add_argument('--max-candidates', type=int, default=None,
                   help='Limit Phase 4 to N candidates (for quick demo)')
    p.add_argument('--skip-phase', type=int, default=0,
                   help='Resume from this phase (loads cached CSVs)')
    p.add_argument('--top-n-optimize', type=int, default=10,
                   help='Number of candidates to optimise in Phase 6')
    p.add_argument('--no-figures', action='store_true',
                   help='Skip figure generation')
    return p.parse_args()


# ── Phase runners with caching ────────────────────────────────────────────────

def phase1(skip: bool = False):
    from src.phase1_target import run_phase1
    logger.info("\n" + "═" * 60)
    logger.info("PHASE 1 — Target Selection & Interface Definition")
    logger.info("═" * 60)
    if skip:
        logger.info("  [SKIP] Phase 1 already completed.")
        return {}
    t0 = time.time()
    result = run_phase1(RECEPTOR_PDB, INTERFACE_JSON, TARGET_PDB_ID, RECEPTOR_CHAIN)
    logger.info(f"  Phase 1 completed in {time.time()-t0:.1f}s")
    return result


def phase2(skip: bool = False) -> pd.DataFrame:
    from src.phase2_generate import run_phase2
    logger.info("\n" + "═" * 60)
    logger.info("PHASE 2 — Candidate Generation (Affibody Library)")
    logger.info("═" * 60)
    if skip and os.path.exists(P2_CSV):
        logger.info(f"  [CACHE] Loading {P2_CSV}")
        return pd.read_csv(P2_CSV)
    t0 = time.time()
    df = run_phase2(output_csv=P2_CSV)
    logger.info(f"  Phase 2 completed in {time.time()-t0:.1f}s → {len(df)} candidates")
    return df


def phase3(df_gen: pd.DataFrame, skip: bool = False) -> pd.DataFrame:
    from src.phase3_filter import run_phase3
    logger.info("\n" + "═" * 60)
    logger.info("PHASE 3 — Developability Filtering")
    logger.info("═" * 60)
    if skip and os.path.exists(P3_CSV):
        logger.info(f"  [CACHE] Loading {P3_CSV}")
        return pd.read_csv(P3_CSV)
    t0 = time.time()
    df = run_phase3(df_gen, output_csv=P3_CSV)
    n_pass = df['passes_filter'].sum()
    logger.info(f"  Phase 3 completed in {time.time()-t0:.1f}s → {n_pass} candidates passed")
    return df


def phase4(df_filt: pd.DataFrame, mock: bool, max_cands,
           skip: bool = False) -> pd.DataFrame:
    from src.phase4_predict import run_phase4
    logger.info("\n" + "═" * 60)
    logger.info(f"PHASE 4 — Structure Prediction ({'mock' if mock else 'ESMFold'})")
    logger.info("═" * 60)
    if skip and os.path.exists(P4_CSV):
        logger.info(f"  [CACHE] Loading {P4_CSV}")
        return pd.read_csv(P4_CSV)
    t0 = time.time()
    df = run_phase4(df_filt, mock=mock, output_csv=P4_CSV,
                    receptor_pdb=RECEPTOR_PDB, max_candidates=max_cands)
    logger.info(f"  Phase 4 completed in {time.time()-t0:.1f}s → {len(df)} structures")
    return df


def phase5(df_struct: pd.DataFrame, skip: bool = False) -> pd.DataFrame:
    from src.phase5_score import run_phase5
    logger.info("\n" + "═" * 60)
    logger.info("PHASE 5 — Composite Scoring & Ranking")
    logger.info("═" * 60)
    if skip and os.path.exists(P5_CSV):
        logger.info(f"  [CACHE] Loading {P5_CSV}")
        return pd.read_csv(P5_CSV)
    t0 = time.time()
    df = run_phase5(df_struct, output_csv=P5_CSV)
    logger.info(f"  Phase 5 completed in {time.time()-t0:.1f}s → top score: {df['final_score'].max():.1f}")
    return df


def phase6(df_ranked: pd.DataFrame, top_n: int, mock: bool,
           skip: bool = False):
    from src.phase6_optimize import run_phase6
    logger.info("\n" + "═" * 60)
    logger.info(f"PHASE 6 — Iterative Optimisation (top {top_n})")
    logger.info("═" * 60)
    if skip and os.path.exists(P6H_CSV) and os.path.exists(P6V_CSV):
        logger.info(f"  [CACHE] Loading optimisation results")
        return pd.read_csv(P6H_CSV), pd.read_csv(P6V_CSV)
    t0 = time.time()
    history, variants = run_phase6(
        df_ranked, top_n=top_n, mock=mock,
        output_history_csv=P6H_CSV,
        output_variants_csv=P6V_CSV,
    )
    logger.info(f"  Phase 6 completed in {time.time()-t0:.1f}s")
    return history, variants


def phase7(df_ranked: pd.DataFrame, skip: bool = False) -> pd.DataFrame:
    from src.phase7_synteins import run_phase7
    logger.info("\n" + "═" * 60)
    logger.info("PHASE 7 — Synteins-Style Analysis (top 5)")
    logger.info("═" * 60)
    if skip and os.path.exists(P7_CSV):
        logger.info(f"  [CACHE] Loading {P7_CSV}")
        return pd.read_csv(P7_CSV)
    t0 = time.time()
    df = run_phase7(df_ranked, top_n=5, output_csv=P7_CSV)
    logger.info(f"  Phase 7 completed in {time.time()-t0:.1f}s")
    return df


def figures(df_filt, df_scored, history_df, synteins_df):
    from src.visualize import generate_all_figures
    fig_dir = os.path.join(RESULTS_DIR, "figures")
    logger.info("\n" + "═" * 60)
    logger.info(f"FIGURES — Generating plots → {fig_dir}")
    logger.info("═" * 60)
    generate_all_figures(
        df_filtered=df_filt,
        df_scored=df_scored,
        history_df=history_df,
        synteins_df=synteins_df,
        output_dir=fig_dir,
    )


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    banner()
    args = parse_args()
    mock = not args.real

    if mock:
        logger.info("Running in MOCK mode (synthetic structure metrics).")
        logger.info("Use --real for ESMFold API predictions (requires internet).\n")
    else:
        logger.info("Running in REAL mode (ESMFold API). This may take 20–60 minutes.\n")

    skip = args.skip_phase
    pipeline_start = time.time()

    # ── Execute phases ────────────────────────────────────────────────────────
    p1_result = phase1(skip=skip > 1)
    df_gen    = phase2(skip=skip > 2)
    df_filt   = phase3(df_gen, skip=skip > 3)
    df_struct = phase4(df_filt, mock=mock, max_cands=args.max_candidates, skip=skip > 4)
    df_ranked = phase5(df_struct, skip=skip > 5)
    history, variants = phase6(df_ranked, top_n=args.top_n_optimize, mock=mock, skip=skip > 6)
    synteins  = phase7(df_ranked, skip=skip > 7)

    # ── Figures ───────────────────────────────────────────────────────────────
    if not args.no_figures:
        figures(df_filt, df_ranked, history, synteins)

    # ── Final summary ─────────────────────────────────────────────────────────
    total_time = time.time() - pipeline_start
    logger.info("\n" + "╔" + "═" * 60 + "╗")
    logger.info("║  PIPELINE COMPLETE                                         ║")
    logger.info("╠" + "═" * 60 + "╣")
    logger.info(f"║  Total runtime  : {total_time:.1f}s")
    logger.info(f"║  Candidates     : {len(df_gen)} generated → "
                f"{df_filt['passes_filter'].sum()} filtered → "
                f"{len(df_struct)} scored")
    logger.info(f"║  Top binder     : {df_ranked['id'].iloc[0]} "
                f"(score {df_ranked['final_score'].iloc[0]:.1f})")
    logger.info(f"║  Top sequence   : {df_ranked['sequence'].iloc[0]}")

    if not history.empty:
        init = df_ranked.head(args.top_n_optimize)['final_score'].mean()
        final_best = history.groupby('candidate_id')['best_score'].max().mean()
        logger.info(f"║  Opt improvement: +{final_best - init:.1f} mean score over "
                    f"{history['round'].max()} rounds")

    logger.info("╠" + "═" * 60 + "╣")
    logger.info(f"║  Results → {RESULTS_DIR}")
    logger.info("╚" + "═" * 60 + "╝\n")

    return df_ranked


if __name__ == "__main__":
    main()
