"""
Phase 6 — Lead Analysis
ADMET profiling, structural alerts, diversity analysis, and
pharmacophore summary for top lead compounds.
"""

import logging
import os
from typing import Dict, List, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from src.smallmol.utils import (
    smiles_to_mol, compute_descriptors, compute_qed, compute_sa_score,
    tanimoto_similarity, count_pharmacophore_features,
)
import smallmol_config as cfg

logger = logging.getLogger(__name__)

# Build Brenk filter catalog (structural alerts / toxicophores)
_brenk_params = FilterCatalogParams()
_brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
BRENK_CATALOG = FilterCatalog(_brenk_params)

# Common CYP inhibitor substructures (simplified)
CYP_PATTERNS = [
    Chem.MolFromSmarts("[#7]1~[#6]~[#6]~[#7]~[#6]~1"),  # imidazole-like
    Chem.MolFromSmarts("c1ccncc1"),                        # pyridine
    Chem.MolFromSmarts("[#7]1~[#6]~[#6]~[#6]~[#7]~1"),   # piperazine-like
]
CYP_PATTERNS = [p for p in CYP_PATTERNS if p is not None]


def analyze_lead(smiles, pocket_profile):
    # type: (str, Dict) -> Dict
    """Comprehensive analysis of a single lead compound."""
    mol = smiles_to_mol(smiles)
    if mol is None:
        return {"error": "invalid SMILES"}

    desc = compute_descriptors(mol)
    qed = compute_qed(mol)
    sa = compute_sa_score(mol)
    pharma = count_pharmacophore_features(mol)

    # ADMET estimates
    tpsa = desc["tpsa"]
    mw = desc["mw"]
    logp = desc["logp"]
    rot = desc["rotatable_bonds"]

    oral_bioavail = "likely" if (rot <= 10 and tpsa <= 140) else "unlikely"
    bbb_permeable = "likely" if (mw < 400 and tpsa < 90 and 1 < logp < 3) else "unlikely"

    # CYP liability
    cyp_flags = 0
    for pat in CYP_PATTERNS:
        if mol.HasSubstructMatch(pat):
            cyp_flags += 1
    cyp_risk = "high" if cyp_flags >= 2 else ("moderate" if cyp_flags == 1 else "low")

    # Structural alerts (Brenk)
    brenk_alerts = []
    entry = BRENK_CATALOG.GetFirstMatch(mol)
    if entry is not None:
        brenk_alerts.append(entry.GetDescription())

    # Solubility estimate (simplified: based on LogP and MW)
    if logp < 2 and mw < 300:
        solubility = "high"
    elif logp < 4 and mw < 500:
        solubility = "moderate"
    else:
        solubility = "low"

    return {
        "smiles": smiles,
        "mw": round(mw, 1),
        "logp": round(logp, 2),
        "hbd": desc["hbd"],
        "hba": desc["hba"],
        "tpsa": round(tpsa, 1),
        "rotatable_bonds": rot,
        "ring_count": desc["ring_count"],
        "aromatic_rings": desc["aromatic_rings"],
        "heavy_atoms": desc["heavy_atoms"],
        "qed": round(qed, 3),
        "sa_score": round(sa, 2),
        "oral_bioavailability": oral_bioavail,
        "bbb_permeability": bbb_permeable,
        "cyp_risk": cyp_risk,
        "solubility_estimate": solubility,
        "structural_alerts": "; ".join(brenk_alerts) if brenk_alerts else "none",
        "n_pharmacophore_hbd": pharma["hbd"],
        "n_pharmacophore_hba": pharma["hba"],
        "n_pharmacophore_hydrophobic": pharma["hydrophobic"],
        "n_pharmacophore_charged": pharma["charged"],
    }


def compute_diversity_matrix(smiles_list):
    # type: (List[str]) -> List[List[float]]
    """Compute pairwise Tanimoto distance matrix."""
    mols = [smiles_to_mol(s) for s in smiles_list]
    n = len(mols)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            if mols[i] and mols[j]:
                sim = tanimoto_similarity(mols[i], mols[j])
                dist = 1.0 - sim
                matrix[i][j] = round(dist, 3)
                matrix[j][i] = round(dist, 3)

    return matrix


def run_phase6(df_scored, variants_df, pocket_profile, output_csv=None):
    # type: (pd.DataFrame, pd.DataFrame, Dict, Optional[str]) -> pd.DataFrame
    """Execute Phase 6: Lead analysis."""
    logger.info("=" * 60)
    logger.info("PHASE 6 — Lead Analysis")
    logger.info("=" * 60)

    # Get top leads: best variant per molecule_id, then overall top-N
    if not variants_df.empty:
        best_variants = (variants_df
                         .sort_values("final_score", ascending=False)
                         .drop_duplicates("molecule_id")
                         .head(cfg.LEAD_TOP_N))
        lead_smiles = best_variants["smiles"].tolist()
        lead_ids = best_variants["molecule_id"].tolist()
    else:
        top = df_scored.head(cfg.LEAD_TOP_N)
        lead_smiles = top["smiles"].tolist()
        lead_ids = top["id"].tolist()

    logger.info("Analyzing %d lead compounds", len(lead_smiles))

    analyses = []
    for lid, smi in zip(lead_ids, lead_smiles):
        result = analyze_lead(smi, pocket_profile)
        result["lead_id"] = lid
        analyses.append(result)

    lead_df = pd.DataFrame(analyses)

    # Diversity
    if len(lead_smiles) > 1:
        dist_matrix = compute_diversity_matrix(lead_smiles)
        mean_dist = 0.0
        count = 0
        for i in range(len(dist_matrix)):
            for j in range(i + 1, len(dist_matrix)):
                mean_dist += dist_matrix[i][j]
                count += 1
        mean_dist /= max(count, 1)
        logger.info("Mean pairwise Tanimoto distance: %.3f", mean_dist)
        lead_df["mean_tanimoto_distance"] = round(mean_dist, 3)

    # Log summary
    for _, row in lead_df.iterrows():
        logger.info("  %s: MW=%.0f LogP=%.1f QED=%.2f SA=%.1f oral=%s cyp=%s alerts=%s",
                     row["lead_id"], row["mw"], row["logp"], row["qed"],
                     row["sa_score"], row["oral_bioavailability"],
                     row["cyp_risk"], row["structural_alerts"])

    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        lead_df.to_csv(output_csv, index=False)
        logger.info("Saved -> %s", output_csv)

    return lead_df
