"""
Phase 4 — Composite Scoring
Score each molecule on drug-likeness (QED), synthetic accessibility,
pharmacophore complementarity to the IL-6 Site II pocket, and
ADMET-related properties.
"""

import logging
import os
from typing import Dict, Optional

import pandas as pd

from src.smallmol.utils import (
    smiles_to_mol, compute_qed, compute_sa_score,
    count_pharmacophore_features,
)
import smallmol_config as cfg

logger = logging.getLogger(__name__)


def _minmax(value, lo, hi):
    # type: (float, float, float) -> float
    """Clip and normalize a value to [0, 1]."""
    return max(0.0, min(1.0, (value - lo) / max(hi - lo, 1e-9)))


def compute_binding_proxy(mol, pocket_profile):
    # type: (...) -> float
    """
    Pharmacophore complementarity score (0–1).
    Measures how well the molecule's features match the pocket's hotspots.
    """
    features = count_pharmacophore_features(mol)

    pocket_hydro = pocket_profile["n_hydrophobic_hotspots"]
    pocket_charged = pocket_profile["n_charged_hotspots"]
    pocket_hbond = pocket_profile["n_hbond_hotspots"]

    # Feature matching: fraction of pocket hotspots covered
    hydro_match = min(features["hydrophobic"], pocket_hydro) / max(pocket_hydro, 1)
    charge_match = min(features["charged"], pocket_charged) / max(pocket_charged, 1)

    hbond_total = features["hbd"] + features["hba"]
    hbond_match = min(hbond_total, pocket_hbond * 2) / max(pocket_hbond * 2, 1)

    # Weighted combination favoring charged contacts (pocket is electrostatic-dominant)
    binding = 0.25 * hydro_match + 0.45 * charge_match + 0.30 * hbond_match

    return min(binding, 1.0)


def score_molecule(smiles, pocket_profile):
    # type: (str, Dict) -> Dict[str, float]
    """Compute all score components for a single molecule."""
    mol = smiles_to_mol(smiles)
    if mol is None:
        return {"qed": 0.0, "sa_score": 10.0, "sa_inv": 0.0,
                "binding_proxy": 0.0, "logp_dev": 4.0, "tpsa_score": 0.0,
                "composite_score": 0.0, "robustness_bonus": 0.0, "final_score": 0.0}

    qed = compute_qed(mol)
    sa = compute_sa_score(mol)
    sa_inv = 1.0 - sa / 10.0  # invert: lower SA = better -> higher sa_inv
    binding = compute_binding_proxy(mol, pocket_profile)

    from rdkit.Chem import Descriptors
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    mw = Descriptors.MolWt(mol)

    # LogP penalty: deviation from ideal ~2.5
    logp_dev = abs(logp - 2.5)

    # TPSA bonus: 60–120 is oral bioavailability sweet spot
    if 60 <= tpsa <= 120:
        tpsa_score = 1.0
    elif 40 <= tpsa < 60 or 120 < tpsa <= 140:
        tpsa_score = 0.5
    else:
        tpsa_score = 0.0

    # Normalize
    norm_qed = _minmax(qed, *cfg.NORM_RANGES["qed"])
    norm_sa = _minmax(sa_inv, *cfg.NORM_RANGES["sa_inv"])
    norm_binding = _minmax(binding, *cfg.NORM_RANGES["binding_proxy"])
    norm_logp = _minmax(logp_dev, *cfg.NORM_RANGES["logp_dev"])
    norm_tpsa = _minmax(tpsa_score, *cfg.NORM_RANGES["tpsa_score"])

    # Composite score (weights from config)
    w = cfg.SCORE_WEIGHTS
    composite = (w["qed"] * norm_qed
                 + w["sa_score"] * norm_sa
                 + w["binding_proxy"] * norm_binding
                 + w["logp_penalty"] * norm_logp
                 + w["tpsa_bonus"] * norm_tpsa)

    composite_100 = max(0.0, min(100.0, composite * 100))

    # Robustness bonus (up to 10 points)
    bonus = 0.0
    if qed > 0.5:
        bonus += 3.3
    if sa < 4.0:
        bonus += 3.3
    if 250 < mw < 450:
        bonus += 3.4

    final = min(100.0, composite_100 + bonus)

    return {
        "qed": round(qed, 4),
        "sa_score": round(sa, 2),
        "sa_inv": round(sa_inv, 4),
        "binding_proxy": round(binding, 4),
        "logp_dev": round(logp_dev, 4),
        "tpsa_score": round(tpsa_score, 2),
        "norm_qed": round(norm_qed, 4),
        "norm_sa": round(norm_sa, 4),
        "norm_binding": round(norm_binding, 4),
        "norm_logp": round(norm_logp, 4),
        "norm_tpsa": round(norm_tpsa, 4),
        "composite_score": round(composite_100, 2),
        "robustness_bonus": round(bonus, 2),
        "final_score": round(final, 2),
    }


def run_phase4(df, pocket_profile, output_csv=None):
    # type: (pd.DataFrame, Dict, Optional[str]) -> pd.DataFrame
    """Execute Phase 4: Composite scoring."""
    logger.info("=" * 60)
    logger.info("PHASE 4 — Composite Scoring")
    logger.info("=" * 60)

    scores = []
    for i, row in df.iterrows():
        s = score_molecule(row["smiles"], pocket_profile)
        scores.append(s)

    score_df = pd.DataFrame(scores)
    df = pd.concat([df.reset_index(drop=True), score_df], axis=1)

    # Rank
    df = df.sort_values("final_score", ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)

    logger.info("Scored %d molecules", len(df))
    logger.info("Score range: %.1f – %.1f", df["final_score"].min(), df["final_score"].max())
    logger.info("Mean score: %.1f", df["final_score"].mean())

    top5 = df.head(5)
    logger.info("Top 5 candidates:")
    for _, r in top5.iterrows():
        logger.info("  %s  score=%.1f  QED=%.2f  SA=%.1f  binding=%.2f  %s",
                     r["id"], r["final_score"], r["qed"], r["sa_score"],
                     r["binding_proxy"], r["smiles"][:60])

    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        df.to_csv(output_csv, index=False)
        logger.info("Saved -> %s", output_csv)

    return df
