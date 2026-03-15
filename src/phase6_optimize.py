"""
Phase 6 — Iterative Optimization (Closed-Loop Design)

Monte Carlo simulated annealing over variable positions:
  - Start from top-10 candidates
  - At each step: propose 1–2 point mutations at variable positions
  - Accept if score improves; accept with Boltzmann probability if worse
  - Temperature annealed from T_start → T_end over rounds
  - Track score trajectory → demonstrates AI-driven design loop

This closes the design loop:
    generate → score → mutate → re-score → accept/reject → repeat

Reference: Metropolis–Hastings MCMC; SA in protein engineering:
  Bloom et al. PNAS 2005; Goldsmith & Tawfik Curr Opin Struct Biol 2012
"""

import copy
import json
import logging
import math
import os
import random
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

from config import (
    VARIABLE_POSITIONS, ALLOWED_AAS, OPTIM_ROUNDS, OPTIM_CANDIDATES,
    OPTIM_MUTATIONS, OPTIM_T_START, OPTIM_T_END, RANDOM_SEED, RESULTS_DIR
)
from src.phase4_predict import predict_candidate
from src.phase5_score import compute_scores

logger = logging.getLogger(__name__)


# ── Mutation proposal ─────────────────────────────────────────────────────────

def propose_mutation(
    sequence: str,
    variable_positions: List[int],
    allowed_aas: List[str],
    n_mutations: int = 1,
    rng: Optional[random.Random] = None,
) -> Tuple[str, List[Tuple[int, str, str]]]:
    """
    Propose 1–n_mutations point mutations at variable positions.
    Returns (new_sequence, [(position, old_aa, new_aa), ...]).
    """
    if rng is None:
        rng = random.Random()

    seq = list(sequence)
    changes = []
    positions = rng.sample(variable_positions, k=min(n_mutations, len(variable_positions)))

    for pos in positions:
        old_aa = seq[pos]
        # Exclude current amino acid from choices
        choices = [aa for aa in allowed_aas if aa != old_aa]
        if not choices:
            continue
        new_aa = rng.choice(choices)
        seq[pos] = new_aa
        changes.append((pos, old_aa, new_aa))

    return ''.join(seq), changes


# ── Single-candidate optimisation ────────────────────────────────────────────

def optimise_candidate(
    candidate_row: pd.Series,
    n_rounds: int = OPTIM_ROUNDS,
    mutations_per_round: int = OPTIM_MUTATIONS,
    t_start: float = OPTIM_T_START,
    t_end: float = OPTIM_T_END,
    mock: bool = True,
    seed: int = RANDOM_SEED,
) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Run simulated annealing for a single candidate.

    Returns:
        history_df : DataFrame with round-by-round best score
        final_variants : List of top variant records
    """
    rng = random.Random(seed + int(abs(hash(candidate_row['sequence'])) % 1000))

    current_seq   = candidate_row['sequence']
    current_score = candidate_row['final_score']
    best_seq      = current_seq
    best_score    = current_score
    cand_id       = candidate_row['id']

    history = [{
        'round': 0,
        'candidate_id': cand_id,
        'sequence': current_seq,
        'score': current_score,
        'best_score': best_score,
        'accepted': True,
        'mutations': '[]',
        'temperature': t_start,
    }]

    # Top variants store
    top_variants = [{'id': cand_id, 'sequence': current_seq, 'score': current_score}]

    for round_idx in range(1, n_rounds + 1):
        # Anneal temperature
        progress = (round_idx - 1) / max(1, n_rounds - 1)
        T = t_start * (t_end / t_start) ** progress

        round_best_score = current_score
        round_best_seq   = current_seq

        for mut_idx in range(mutations_per_round):
            # Propose 1 or 2 mutations (more exploration early, less late)
            n_muts = 2 if round_idx <= n_rounds // 2 else 1
            new_seq, changes = propose_mutation(
                current_seq, VARIABLE_POSITIONS, ALLOWED_AAS, n_muts, rng
            )

            # Score the mutant
            struct_metrics = predict_candidate(new_seq, n_seeds=1, mock=mock)

            # Build a mini DataFrame for score computation
            row_dict = candidate_row.to_dict()
            row_dict.update(struct_metrics)
            row_dict['sequence']   = new_seq
            row_dict['max_aggregation'] = max(
                [row_dict.get('max_aggregation', 1.0)] +
                [struct_metrics.get('bsa_proxy_A2', 0) / 200]  # rough proxy
            )
            mini_df = pd.DataFrame([row_dict])

            # Re-use Phase 5 scoring
            try:
                scored_mini = compute_scores(mini_df)
                new_score = float(scored_mini['final_score'].iloc[0])
            except Exception:
                new_score = current_score * 0.95  # fallback

            # Metropolis acceptance
            delta = new_score - current_score
            accept = (delta > 0) or (rng.random() < math.exp(delta / T))

            if accept:
                current_seq   = new_seq
                current_score = new_score

            if new_score > round_best_score:
                round_best_score = new_score
                round_best_seq   = new_seq

            if new_score > best_score:
                best_score = new_score
                best_seq   = new_seq
                top_variants.append({
                    'id': f"{cand_id}_R{round_idx}M{mut_idx+1}",
                    'sequence': new_seq,
                    'score': new_score,
                })

        history.append({
            'round': round_idx,
            'candidate_id': cand_id,
            'sequence': round_best_seq,
            'score': round_best_score,
            'best_score': best_score,
            'accepted': True,
            'mutations': json.dumps([c[1] + str(c[0]) + c[2] for c in []]),
            'temperature': T,
        })

    history_df = pd.DataFrame(history)

    # Keep top-5 variants per candidate
    top_variants.sort(key=lambda x: x['score'], reverse=True)
    top_variants = top_variants[:5]

    improvement = best_score - candidate_row['final_score']
    logger.info(f"  {cand_id}: {candidate_row['final_score']:.1f} → {best_score:.1f} "
                f"(+{improvement:.1f} over {n_rounds} rounds)")

    return history_df, top_variants


# ── Phase 6 runner ────────────────────────────────────────────────────────────

def run_phase6(
    df: pd.DataFrame,
    top_n: int = OPTIM_CANDIDATES,
    mock: bool = True,
    output_history_csv: Optional[str] = None,
    output_variants_csv: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Execute Phase 6: iterative optimisation for the top-N candidates.

    Returns:
        all_history  : DataFrame of per-round scores for all candidates
        top_variants : DataFrame of best variants discovered
    """
    logger.info("─" * 50)
    logger.info(f"PHASE 6 — Iterative optimisation (MC/SA)")
    logger.info(f"  Candidates  : top {top_n}")
    logger.info(f"  Rounds      : {OPTIM_ROUNDS}")
    logger.info(f"  Mutations/round : {OPTIM_MUTATIONS}")
    logger.info(f"  Temperature : {OPTIM_T_START} → {OPTIM_T_END} (annealing)")
    logger.info("─" * 50)

    top_candidates = df.head(top_n)
    all_histories  = []
    all_variants   = []

    for i, (_, row) in enumerate(tqdm(top_candidates.iterrows(),
                                      total=top_n, desc="Optimising")):
        seed = RANDOM_SEED + i * 100
        history_df, variants = optimise_candidate(
            row,
            n_rounds=OPTIM_ROUNDS,
            mutations_per_round=OPTIM_MUTATIONS,
            mock=mock,
            seed=seed,
        )
        all_histories.append(history_df)
        all_variants.extend(variants)

    all_history_df  = pd.concat(all_histories, ignore_index=True)
    all_variants_df = pd.DataFrame(all_variants).drop_duplicates('sequence')
    all_variants_df.sort_values('score', ascending=False, inplace=True)
    all_variants_df.reset_index(drop=True, inplace=True)

    # Summary
    initial_mean = top_candidates['final_score'].mean()
    final_mean   = all_history_df.groupby('candidate_id')['best_score'].max().mean()
    logger.info(f"\nOptimisation complete:")
    logger.info(f"  Mean initial score : {initial_mean:.1f}")
    logger.info(f"  Mean final score   : {final_mean:.1f}")
    logger.info(f"  Mean improvement   : +{final_mean - initial_mean:.1f}")
    logger.info(f"  Total variants explored: {len(all_variants_df)}")

    if output_history_csv:
        os.makedirs(os.path.dirname(output_history_csv), exist_ok=True)
        all_history_df.to_csv(output_history_csv, index=False)
        logger.info(f"Saved optimisation history → {output_history_csv}")

    if output_variants_csv:
        os.makedirs(os.path.dirname(output_variants_csv), exist_ok=True)
        all_variants_df.to_csv(output_variants_csv, index=False)
        logger.info(f"Saved optimised variants → {output_variants_csv}")

    return all_history_df, all_variants_df
