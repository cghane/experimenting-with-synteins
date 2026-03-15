"""
Phase 3 — First-Pass Developability Filtering

Filters 1000 → ~200 candidates based on physicochemical properties
relevant to therapeutic development. Each filter is independently
motivated by known challenges in biologics manufacturing and safety.

Filter criteria:
  1. Hydrophobicity     : No window > 2.5 (KD scale) — prevents aggregation
  2. Net charge         : −3 to +6 at pH 7.4 — solubility window
  3. Cysteines          : ≤1 — avoids off-target disulfides
  4. Low complexity     : Min Shannon entropy ≥ 2.0 — removes repetitive motifs
  5. Aggregation        : Max PASTA window ≤ 2.0 — reduces amyloid risk
  6. Helix propensity   : Mean ≥ 0.65 — structural stability of 3-helix bundle
  7. Molecular weight   : 4–12 kDa — within target size range
"""

import logging
import os
from typing import Optional, Dict

import pandas as pd

from config import (
    MAX_HYDROPHOBICITY_WINDOW, MIN_NET_CHARGE, MAX_NET_CHARGE,
    MAX_CYSTEINES, MAX_AGGREGATION_WINDOW, MIN_COMPLEXITY_ENTROPY, RESULTS_DIR
)
from src.utils import (
    max_hydrophobic_window, net_charge, count_cysteines,
    min_complexity, max_aggregation_window, mean_helix_propensity,
    molecular_weight
)

logger = logging.getLogger(__name__)

# ── Individual filter functions ───────────────────────────────────────────────

def filter_hydrophobicity(seq: str) -> bool:
    """Pass if no 5-aa window exceeds KD threshold (avoids membrane insertion)."""
    return max_hydrophobic_window(seq, window=5) <= MAX_HYDROPHOBICITY_WINDOW


def filter_charge(seq: str) -> bool:
    """Pass if net charge is within the therapeutic solubility window."""
    q = net_charge(seq)
    return MIN_NET_CHARGE <= q <= MAX_NET_CHARGE


def filter_cysteines(seq: str) -> bool:
    """Pass if cysteine count ≤ threshold (≤1 avoids scrambled disulfides)."""
    return count_cysteines(seq) <= MAX_CYSTEINES


def filter_complexity(seq: str) -> bool:
    """Pass if minimum local Shannon entropy ≥ threshold (removes repetitive seqs)."""
    return min_complexity(seq, window=12) >= MIN_COMPLEXITY_ENTROPY


def filter_aggregation(seq: str) -> bool:
    """Pass if maximum PASTA aggregation window ≤ threshold."""
    return max_aggregation_window(seq, window=5) <= MAX_AGGREGATION_WINDOW


def filter_helix_propensity(seq: str, min_propensity: float = 0.65) -> bool:
    """Pass if mean helix propensity is sufficient for 3-helix bundle stability."""
    return mean_helix_propensity(seq) >= min_propensity


def filter_molecular_weight(seq: str, min_kda: float = 4.0, max_kda: float = 12.0) -> bool:
    """Pass if MW is within typical mini-protein range."""
    mw = molecular_weight(seq) / 1000  # convert to kDa
    return min_kda <= mw <= max_kda


# ── Composite filter ──────────────────────────────────────────────────────────

FILTERS: Dict[str, callable] = {
    'hydrophobicity':    filter_hydrophobicity,
    'charge':            filter_charge,
    'cysteines':         filter_cysteines,
    'low_complexity':    filter_complexity,
    'aggregation':       filter_aggregation,
    'helix_propensity':  filter_helix_propensity,
    'molecular_weight':  filter_molecular_weight,
}


def apply_filters(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply all developability filters to candidate DataFrame.
    Adds per-filter pass/fail columns and a summary 'passes_filter' column.
    Also computes additional physicochemical metrics for downstream use.
    """
    logger.info(f"Applying {len(FILTERS)} developability filters to {len(df)} candidates …")

    # Compute metrics columns
    df = df.copy()
    df['max_hydro_window']   = df['sequence'].apply(lambda s: max_hydrophobic_window(s, 5))
    df['net_charge_calc']    = df['sequence'].apply(net_charge)
    df['n_cysteines']        = df['sequence'].apply(count_cysteines)
    df['min_entropy']        = df['sequence'].apply(lambda s: min_complexity(s, 12))
    df['max_aggregation']    = df['sequence'].apply(lambda s: max_aggregation_window(s, 5))
    df['mean_helix_prop']    = df['sequence'].apply(mean_helix_propensity)
    df['mol_weight_kda']     = df['sequence'].apply(lambda s: molecular_weight(s) / 1000)

    # Apply each filter
    filter_results = {}
    for name, fn in FILTERS.items():
        filter_results[f'pass_{name}'] = df['sequence'].apply(fn)

    filter_df = pd.DataFrame(filter_results, index=df.index)
    df = pd.concat([df, filter_df], axis=1)

    # Combined pass
    filter_cols = [c for c in df.columns if c.startswith('pass_')]
    df['passes_filter'] = df[filter_cols].all(axis=1)

    # Log breakdown
    total = len(df)
    passed = df['passes_filter'].sum()
    logger.info(f"Filter results ({total} → {passed} candidates):")
    for name in FILTERS:
        col = f'pass_{name}'
        n_pass = df[col].sum()
        logger.info(f"  {name:20s}: {n_pass:4d}/{total} pass ({100*n_pass/total:.1f}%)")

    logger.info(f"  {'COMBINED':20s}: {passed:4d}/{total} pass ({100*passed/total:.1f}%)")

    return df


def run_phase3(df: pd.DataFrame, output_csv: Optional[str] = None) -> pd.DataFrame:
    """Execute Phase 3 filtering and optionally save to CSV."""
    logger.info("─" * 50)
    logger.info("PHASE 3 — Developability filtering")
    logger.info("─" * 50)

    df_filtered = apply_filters(df)
    passed = df_filtered[df_filtered['passes_filter']]

    logger.info(f"\nPhase 3 complete: {len(passed)} candidates passed all filters.")

    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        df_filtered.to_csv(output_csv, index=False)
        logger.info(f"Saved filtered results → {output_csv}")

    return df_filtered
