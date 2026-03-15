"""
Phase 5 — Composite Scoring Function

Combines structural and sequence-based metrics into a single
interpretable score (0–100) for candidate ranking.

Score formula:
    S = w1 × norm(interface_pLDDT)
      + w2 × norm(BSA_proxy)
      + w3 × norm(H-bond_estimate)
      - w4 × norm(aggregation_penalty)
      - w5 × norm(protease_penalty)

Weights defined in config.py. All terms normalised to [0, 1] before
weighting to ensure commensurability.

Developability bonus: small bonus for charge window compliance,
  helix propensity, and sequence entropy (robustness score).
"""

import logging
import os
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from config import SCORE_WEIGHTS, RESULTS_DIR
from src.utils import (
    total_protease_sites, max_aggregation_window, immunogenicity_score
)

logger = logging.getLogger(__name__)


# ── Normalisation helpers ─────────────────────────────────────────────────────

def _minmax(series: pd.Series, lo: float, hi: float) -> pd.Series:
    """Clip and normalise a series to [0, 1] using provided lo/hi bounds."""
    clipped = series.clip(lower=lo, upper=hi)
    if hi == lo:
        return pd.Series(0.5, index=series.index)
    return (clipped - lo) / (hi - lo)


# Reference ranges for normalisation (from literature + pilot experiments)
NORM_RANGES = {
    'interface_plddt':   (40.0, 92.0),   # ESMFold pLDDT range for mini-proteins
    'bsa_proxy_A2':      (0.0,  700.0),  # Å²; typical Affibody interface ~400–600 Å²
    'hbond_estimate':    (0.0,  10.0),   # estimated H-bonds
    'max_aggregation':   (0.0,   3.0),   # PASTA scale
    'protease_risk':     (0.0,  12.0),   # total cleavage site count
    'immunogen':         (0.0,   1.0),   # 0–1 score
}


# ── Per-candidate score computation ──────────────────────────────────────────

def compute_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute composite binding score for each candidate.
    Assumes df has columns from Phase 4 (interface_plddt, bsa_proxy_A2, etc.)
    and Phase 3 (max_aggregation, sequence, …).

    Adds columns:
        norm_*, score_components, composite_score, robustness_score,
        final_score, rank
    """
    df = df.copy()

    # ── Compute missing columns if needed ────────────────────────────────────
    if 'max_aggregation' not in df.columns:
        df['max_aggregation'] = df['sequence'].apply(lambda s: max_aggregation_window(s, 5))

    df['protease_risk']   = df['sequence'].apply(total_protease_sites)
    df['immunogen_score'] = df['sequence'].apply(immunogenicity_score)

    # ── Normalise each term ───────────────────────────────────────────────────
    df['norm_interface_plddt'] = _minmax(
        df['interface_plddt'], *NORM_RANGES['interface_plddt'])
    df['norm_bsa']             = _minmax(
        df['bsa_proxy_A2'],    *NORM_RANGES['bsa_proxy_A2'])
    df['norm_hbond']           = _minmax(
        df['hbond_estimate'],  *NORM_RANGES['hbond_estimate'])
    df['norm_aggregation']     = _minmax(
        df['max_aggregation'], *NORM_RANGES['max_aggregation'])
    df['norm_protease']        = _minmax(
        df['protease_risk'],   *NORM_RANGES['protease_risk'])
    df['norm_immunogen']       = _minmax(
        df['immunogen_score'], *NORM_RANGES['immunogen'])

    # ── Composite binding score ───────────────────────────────────────────────
    w = SCORE_WEIGHTS
    df['composite_score'] = (
        w['interface_plddt'] * df['norm_interface_plddt'] +
        w['bsa_score']       * df['norm_bsa'] +
        w['hbond_score']     * df['norm_hbond'] +
        w['aggregation_pen'] * df['norm_aggregation'] +    # negative weight in config
        w['protease_pen']    * df['norm_protease']         # negative weight in config
    )

    # Scale to 0–100
    df['composite_score'] = (df['composite_score'] * 100).clip(0, 100)

    # ── Robustness bonus (developability) ────────────────────────────────────
    # Bonus for:  good helix propensity, charge in window, low immunogenicity
    charge_ok    = ((df['net_charge_calc'] >= -3) & (df['net_charge_calc'] <= 6)).astype(float)
    helix_ok     = (df['mean_helix_prop'] >= 0.70).astype(float)
    immun_ok     = (df['norm_immunogen'] <= 0.4).astype(float)

    df['robustness_score'] = (charge_ok + helix_ok + immun_ok) / 3.0 * 10.0  # 0–10 bonus

    # ── Final score ───────────────────────────────────────────────────────────
    df['final_score'] = (df['composite_score'] + df['robustness_score']).clip(0, 100)

    # ── Rank ─────────────────────────────────────────────────────────────────
    df['rank'] = df['final_score'].rank(method='dense', ascending=False).astype(int)
    df.sort_values('rank', inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


def summarise_scores(df: pd.DataFrame) -> str:
    """Return a formatted leaderboard string of top candidates."""
    top = df.head(20)
    lines = [
        "\n" + "=" * 80,
        f"{'Rank':>4}  {'ID':>10}  {'Final':>6}  {'Comp':>6}  "
        f"{'iPLDDT':>7}  {'BSA Å²':>7}  {'nCont':>6}  {'Protease':>8}  Sequence",
        "-" * 80,
    ]
    for _, row in top.iterrows():
        seq_display = row['sequence'][:20] + "…"
        lines.append(
            f"{row['rank']:>4}  {row['id']:>10}  {row['final_score']:>6.1f}  "
            f"{row['composite_score']:>6.1f}  "
            f"{row['interface_plddt']:>7.1f}  {row['bsa_proxy_A2']:>7.0f}  "
            f"{int(row['n_contacts']):>6}  {int(row['protease_risk']):>8}  {seq_display}"
        )
    lines.append("=" * 80)
    return '\n'.join(lines)


def run_phase5(df: pd.DataFrame, output_csv: Optional[str] = None) -> pd.DataFrame:
    """Execute Phase 5: scoring and ranking."""
    logger.info("─" * 50)
    logger.info("PHASE 5 — Composite scoring & ranking")
    logger.info(f"Score = {SCORE_WEIGHTS['interface_plddt']} × iPLDDT "
                f"+ {SCORE_WEIGHTS['bsa_score']} × BSA "
                f"+ {SCORE_WEIGHTS['hbond_score']} × Hbond "
                f"+ {SCORE_WEIGHTS['aggregation_pen']} × Aggr "
                f"+ {SCORE_WEIGHTS['protease_pen']} × Protease")
    logger.info("─" * 50)

    scored = compute_scores(df)

    logger.info(f"\nScoring complete:")
    logger.info(f"  Top score   : {scored['final_score'].max():.1f}")
    logger.info(f"  Mean score  : {scored['final_score'].mean():.1f} ± "
                f"{scored['final_score'].std():.1f}")
    logger.info(f"  Median score: {scored['final_score'].median():.1f}")
    logger.info(summarise_scores(scored))

    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        scored.to_csv(output_csv, index=False)
        logger.info(f"Saved ranked candidates → {output_csv}")

    return scored
