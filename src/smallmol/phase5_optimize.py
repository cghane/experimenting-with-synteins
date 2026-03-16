"""
Phase 5 — Monte Carlo Optimization
Simulated annealing with SMARTS-based molecular mutations.
Iteratively improves top candidates using Metropolis acceptance.
"""

import logging
import math
import os
import random
from typing import Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from src.smallmol.utils import smiles_to_mol, mol_to_canonical_smiles
from src.smallmol.phase4_score import score_molecule
import smallmol_config as cfg

logger = logging.getLogger(__name__)

# SMARTS-based mutation transforms
# Each is (reaction_smarts, name, description)
MUTATION_TRANSFORMS = [
    # Bioisosteric replacements
    ("[CH3:1]>>[NH2:1]", "methyl_to_amine", "CH3 -> NH2"),
    ("[NH2:1]>>[CH3:1]", "amine_to_methyl", "NH2 -> CH3"),
    ("[OH:1]>>[F:1]", "hydroxyl_to_F", "OH -> F"),
    ("[F:1]>>[OH:1]", "F_to_hydroxyl", "F -> OH"),
    ("[F:1]>>[Cl:1]", "F_to_Cl", "F -> Cl"),
    ("[Cl:1]>>[F:1]", "Cl_to_F", "Cl -> F"),
    ("[OH:1]>>[NH2:1]", "hydroxyl_to_amine", "OH -> NH2"),
    ("[NH2:1]>>[OH:1]", "amine_to_hydroxyl", "NH2 -> OH"),
    # Ring modifications
    ("[cH:1]>>[nH0:1]", "CH_to_N_in_ring", "aromatic CH -> N"),
    ("[nH0:1]>>[cH:1]", "N_to_CH_in_ring", "aromatic N -> CH"),
    # Functional group additions/removals
    ("[cH:1]>>[c:1]O", "add_hydroxyl", "add -OH to aromatic"),
    ("[cH:1]>>[c:1]F", "add_fluorine", "add -F to aromatic"),
    ("[cH:1]>>[c:1]C", "add_methyl", "add -CH3 to aromatic"),
    ("[cH:1]>>[c:1]OC", "add_methoxy", "add -OCH3 to aromatic"),
    ("[cH:1]>>[c:1]N", "add_amine", "add -NH2 to aromatic"),
]

# Pre-compile reactions
COMPILED_MUTATIONS = []
for smarts, name, desc in MUTATION_TRANSFORMS:
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is not None:
            COMPILED_MUTATIONS.append((rxn, name, desc))
    except Exception:
        pass


def propose_mutation(smiles, rng):
    # type: (str, random.Random) -> Tuple[Optional[str], str]
    """
    Propose a single mutation to a molecule.
    Returns (new_smiles, mutation_name) or (None, "") on failure.
    """
    mol = smiles_to_mol(smiles)
    if mol is None or not COMPILED_MUTATIONS:
        return None, ""

    # Shuffle and try mutations until one works
    candidates = list(COMPILED_MUTATIONS)
    rng.shuffle(candidates)

    for rxn, name, desc in candidates:
        try:
            products = rxn.RunReactants((mol,))
            if products:
                # Pick a random product
                product_set = list(products)
                rng.shuffle(product_set)
                for prod_tuple in product_set:
                    prod = prod_tuple[0]
                    try:
                        Chem.SanitizeMol(prod)
                        new_smi = mol_to_canonical_smiles(prod)
                        if new_smi and new_smi != smiles:
                            return new_smi, name
                    except Exception:
                        continue
        except Exception:
            continue

    return None, ""


def optimize_candidate(candidate_id, smiles, pocket_profile, rng,
                       n_rounds, mutations_per_round, t_start, t_end):
    # type: (str, str, Dict, random.Random, int, int, float, float) -> Tuple[List[Dict], List[Dict]]
    """
    Optimize a single molecule via simulated annealing.
    Returns (history, top_variants).
    """
    current_smi = smiles
    current_scores = score_molecule(current_smi, pocket_profile)
    current_score = current_scores["final_score"]
    best_score = current_score
    best_smi = current_smi

    history = []
    variants = [(current_score, current_smi)]

    for rnd in range(1, n_rounds + 1):
        progress = (rnd - 1) / max(n_rounds - 1, 1)
        temp = t_start * (t_end / max(t_start, 1e-9)) ** progress

        accepted_this_round = 0
        for _ in range(mutations_per_round):
            new_smi, mut_name = propose_mutation(current_smi, rng)
            if new_smi is None:
                continue

            new_scores = score_molecule(new_smi, pocket_profile)
            new_score = new_scores["final_score"]
            delta = new_score - current_score

            # Metropolis acceptance
            accept = False
            if delta > 0:
                accept = True
            elif temp > 0:
                prob = math.exp(delta / max(temp, 1e-9))
                accept = rng.random() < prob

            if accept:
                current_smi = new_smi
                current_score = new_score
                accepted_this_round += 1

                if current_score > best_score:
                    best_score = current_score
                    best_smi = current_smi

                variants.append((current_score, current_smi))

        history.append({
            "molecule_id": candidate_id,
            "round": rnd,
            "temperature": round(temp, 4),
            "best_score": round(best_score, 2),
            "current_score": round(current_score, 2),
            "accepted": accepted_this_round,
            "best_smiles": best_smi,
        })

    # Deduplicate and get top-5 variants
    seen = set()
    unique_variants = []
    for score, smi in sorted(variants, key=lambda x: -x[0]):
        if smi not in seen:
            seen.add(smi)
            unique_variants.append({"molecule_id": candidate_id, "smiles": smi,
                                     "final_score": score})
        if len(unique_variants) >= 5:
            break

    return history, unique_variants


def run_phase5(df, pocket_profile, output_history_csv=None, output_variants_csv=None):
    # type: (pd.DataFrame, Dict, Optional[str], Optional[str]) -> Tuple[pd.DataFrame, pd.DataFrame]
    """Execute Phase 5: Monte Carlo optimization."""
    logger.info("=" * 60)
    logger.info("PHASE 5 — Monte Carlo Optimization")
    logger.info("=" * 60)

    top_n = min(cfg.OPTIM_TOP_N, len(df))
    df_top = df.head(top_n).copy()

    logger.info("Optimizing top %d candidates", top_n)
    logger.info("Rounds: %d, Mutations/round: %d, T: %.2f -> %.2f",
                cfg.OPTIM_ROUNDS, cfg.OPTIM_MUTATIONS_PER_ROUND,
                cfg.OPTIM_T_START, cfg.OPTIM_T_END)

    rng = random.Random(cfg.RANDOM_SEED + 100)
    all_history = []
    all_variants = []

    for i, row in df_top.iterrows():
        cid = row["id"]
        smi = row["smiles"]
        initial_score = row["final_score"]

        history, variants = optimize_candidate(
            cid, smi, pocket_profile, rng,
            cfg.OPTIM_ROUNDS, cfg.OPTIM_MUTATIONS_PER_ROUND,
            cfg.OPTIM_T_START, cfg.OPTIM_T_END,
        )
        all_history.extend(history)
        all_variants.extend(variants)

        best = max(v["final_score"] for v in variants)
        delta = best - initial_score
        logger.info("  %s: %.1f -> %.1f (%+.1f)", cid, initial_score, best, delta)

    history_df = pd.DataFrame(all_history)
    variants_df = pd.DataFrame(all_variants)

    mean_improvement = 0.0
    if not history_df.empty:
        initial_scores = df_top.set_index("id")["final_score"]
        final_scores = history_df.groupby("molecule_id")["best_score"].max()
        improvements = final_scores - initial_scores.reindex(final_scores.index)
        mean_improvement = improvements.mean()

    logger.info("Mean improvement: %+.1f score units", mean_improvement)

    if output_history_csv:
        os.makedirs(os.path.dirname(output_history_csv), exist_ok=True)
        history_df.to_csv(output_history_csv, index=False)
        logger.info("Saved history -> %s", output_history_csv)

    if output_variants_csv:
        os.makedirs(os.path.dirname(output_variants_csv), exist_ok=True)
        variants_df.to_csv(output_variants_csv, index=False)
        logger.info("Saved variants -> %s", output_variants_csv)

    return history_df, variants_df
