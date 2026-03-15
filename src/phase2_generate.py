"""
Phase 2 — Candidate Mini-Protein Generation

Strategy: Affibody Z-domain scaffold (55 aa, 3-helix bundle).
  - 13 variable positions on the binding face of helices 1 & 2
  - Combinatorial sampling biased toward IL-6 Site II complementarity
  - Two generation modes:
      (a) Biased design — residue probabilities from interface_residues.json
      (b) Random exploration — uniform sampling for diversity

Reference: Löfblom et al. FEBS Letters 2010; Nygren lab affibody library.
"""

import json
import logging
import random
import os
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd

from config import (
    AFFIBODY_SCAFFOLD, VARIABLE_POSITIONS, ALLOWED_AAS,
    N_CANDIDATES_GENERATE, RANDOM_SEED, INTERFACE_JSON, RESULTS_DIR
)
from src.utils import (
    mean_hydrophobicity, net_charge, mean_helix_propensity,
    sequence_similarity, sequence_entropy
)

logger = logging.getLogger(__name__)

# ── Amino acid design probabilities ──────────────────────────────────────────

def build_position_distributions(interface_json: str) -> Dict[int, Dict[str, float]]:
    """
    Build position-specific amino acid probability distributions
    informed by IL-6 Site II complementarity analysis.
    Positions where we have design rationale receive biased probabilities;
    others receive uniform distribution over ALLOWED_AAS.
    """
    with open(interface_json) as f:
        data = json.load(f)

    bias = data.get('scaffold_variable_positions', {}).get('position_design_bias', {})

    distributions: Dict[int, Dict[str, float]] = {}
    for pos in VARIABLE_POSITIONS:
        pos_str = str(pos)
        if pos_str in bias:
            preferred = bias[pos_str]['preferred']
            # Preferred AAs get 10x weight, others get 1x
            weights = {aa: (10.0 if aa in preferred else 1.0) for aa in ALLOWED_AAS}
        else:
            weights = {aa: 1.0 for aa in ALLOWED_AAS}

        total = sum(weights.values())
        distributions[pos] = {aa: w / total for aa, w in weights.items()}

    return distributions


def sample_sequence(
    scaffold: str,
    variable_positions: List[int],
    distributions: Dict[int, Dict[str, float]],
    rng: random.Random
) -> Tuple[str, Dict[int, str]]:
    """
    Sample a candidate sequence by substituting variable positions
    according to the given probability distributions.
    Returns (sequence, {position: amino_acid}).
    """
    seq = list(scaffold)
    substitutions: Dict[int, str] = {}

    for pos in variable_positions:
        dist = distributions[pos]
        aas = list(dist.keys())
        probs = list(dist.values())
        chosen = rng.choices(aas, weights=probs, k=1)[0]
        seq[pos] = chosen
        substitutions[pos] = chosen

    return ''.join(seq), substitutions


# ── Sequence diversity ────────────────────────────────────────────────────────

def pairwise_diversity(sequences: List[str]) -> float:
    """
    Mean pairwise Hamming distance normalised by sequence length.
    Measures library diversity (0 = identical, 1 = maximally diverse).
    """
    if len(sequences) < 2:
        return 0.0
    n = len(sequences[0])
    total = 0
    count = 0
    # Sample up to 500 pairs for speed
    sample_size = min(len(sequences), 100)
    sampled = sequences[:sample_size]
    for i in range(len(sampled)):
        for j in range(i + 1, len(sampled)):
            total += sum(a != b for a, b in zip(sampled[i], sampled[j]))
            count += 1
    return (total / count / n) if count else 0.0


# ── Candidate generation ──────────────────────────────────────────────────────

def generate_candidates(
    n: int = N_CANDIDATES_GENERATE,
    biased_fraction: float = 0.7,
    seed: int = RANDOM_SEED,
    interface_json: str = INTERFACE_JSON,
) -> pd.DataFrame:
    """
    Generate n candidate Affibody sequences for IL-6 binding.

    Args:
        n: Number of candidates to generate.
        biased_fraction: Fraction sampled with interface-informed bias (vs uniform).
        seed: Random seed for reproducibility.
        interface_json: Path to interface residues JSON.

    Returns:
        DataFrame with columns: [id, sequence, substitutions,
                                  mean_hydrophobicity, net_charge,
                                  helix_propensity, entropy, is_biased]
    """
    rng = random.Random(seed)
    np.random.seed(seed)

    # Build biased and uniform distributions
    biased_dist  = build_position_distributions(interface_json)
    uniform_dist = {pos: {aa: 1.0 / len(ALLOWED_AAS) for aa in ALLOWED_AAS}
                    for pos in VARIABLE_POSITIONS}

    records = []
    seen_sequences = set()
    n_biased  = int(n * biased_fraction)
    n_uniform = n - n_biased

    logger.info(f"Generating {n_biased} biased + {n_uniform} random candidates …")

    for idx in range(n):
        is_biased = (idx < n_biased)
        dist = biased_dist if is_biased else uniform_dist

        # Rejection-sample for uniqueness (max 100 attempts per candidate)
        for _ in range(100):
            seq, subs = sample_sequence(AFFIBODY_SCAFFOLD, VARIABLE_POSITIONS, dist, rng)
            if seq not in seen_sequences:
                seen_sequences.add(seq)
                break
        else:
            logger.debug(f"Duplicate sequence at candidate {idx}, accepting anyway")

        records.append({
            'id':               f"ABIO-{idx+1:04d}",
            'sequence':         seq,
            'substitutions':    json.dumps(subs),
            'scaffold':         AFFIBODY_SCAFFOLD,
            'mean_hydrophobicity': mean_hydrophobicity(seq),
            'net_charge':       net_charge(seq),
            'helix_propensity': mean_helix_propensity(seq),
            'entropy':          sequence_entropy(seq),
            'is_biased':        is_biased,
            'length':           len(seq),
        })

    df = pd.DataFrame(records)

    # Summary stats
    diversity = pairwise_diversity(df['sequence'].tolist())
    logger.info(f"Generated {len(df)} candidates (library diversity: {diversity:.3f})")
    logger.info(f"  Hydrophobicity range : {df['mean_hydrophobicity'].min():.2f} – "
                f"{df['mean_hydrophobicity'].max():.2f}")
    logger.info(f"  Net charge range     : {df['net_charge'].min():.1f} – "
                f"{df['net_charge'].max():.1f}")

    return df


def run_phase2(output_csv: Optional[str] = None) -> pd.DataFrame:
    """Execute Phase 2 and optionally save to CSV."""
    logger.info("─" * 50)
    logger.info("PHASE 2 — Candidate generation (Affibody scaffold)")
    logger.info(f"  Scaffold : {AFFIBODY_SCAFFOLD}")
    logger.info(f"  Length   : {len(AFFIBODY_SCAFFOLD)} aa")
    logger.info(f"  Variable : {len(VARIABLE_POSITIONS)} positions → {VARIABLE_POSITIONS}")
    logger.info("─" * 50)

    df = generate_candidates()

    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        df.to_csv(output_csv, index=False)
        logger.info(f"Saved {len(df)} candidates → {output_csv}")

    return df
