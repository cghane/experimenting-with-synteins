"""
Phase 3 — Molecular Filtering
Apply Lipinski Rule of 5, Veber rules, PAINS, and drug-likeness filters
to reduce the candidate library to viable leads.
"""

import logging
import os
from typing import Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from src.smallmol.utils import smiles_to_mol
import smallmol_config as cfg

logger = logging.getLogger(__name__)

# Build PAINS filter catalog once at module load
_pains_params = FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
PAINS_CATALOG = FilterCatalog(_pains_params)


def _check_lipinski(row):
    # type: (pd.Series) -> bool
    """Lipinski Rule of 5: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10."""
    return (row["mw"] <= cfg.LIPINSKI_MW_MAX
            and row["logp"] <= cfg.LIPINSKI_LOGP_MAX
            and row["hbd"] <= cfg.LIPINSKI_HBD_MAX
            and row["hba"] <= cfg.LIPINSKI_HBA_MAX)


def _check_veber(row):
    # type: (pd.Series) -> bool
    """Veber rules: TPSA <= 140, rotatable bonds <= 10."""
    return (row["tpsa"] <= cfg.VEBER_TPSA_MAX
            and row["rotatable_bonds"] <= cfg.VEBER_ROTBONDS_MAX)


def _check_complexity(row):
    # type: (pd.Series) -> bool
    """Minimum molecular complexity: at least 1 ring, MW >= 150."""
    return row["ring_count"] >= cfg.MIN_RINGS and row["mw"] >= cfg.MIN_MW


def _check_aggregation(row):
    # type: (pd.Series) -> bool
    """Reject likely colloidal aggregators (high LogP + high MW)."""
    if row["logp"] > cfg.AGGREGATION_LOGP_MIN and row["mw"] > cfg.AGGREGATION_MW_MIN:
        return False
    return True


def _check_pains(smiles):
    # type: (str) -> bool
    """Check for PAINS (pan-assay interference) substructures. Returns True if clean."""
    mol = smiles_to_mol(smiles)
    if mol is None:
        return False
    return not PAINS_CATALOG.HasMatch(mol)


FILTERS = {
    "lipinski":    _check_lipinski,
    "veber":       _check_veber,
    "complexity":  _check_complexity,
    "aggregation": _check_aggregation,
}


def apply_filters(df):
    # type: (pd.DataFrame) -> pd.DataFrame
    """Apply all filters and add pass/fail columns."""
    df = df.copy()

    # Row-based filters
    for name, func in FILTERS.items():
        col = "pass_%s" % name
        df[col] = df.apply(func, axis=1)

    # PAINS filter (needs SMILES directly)
    df["pass_pains"] = df["smiles"].apply(_check_pains)

    # Overall pass
    filter_cols = ["pass_%s" % n for n in FILTERS] + ["pass_pains"]
    df["passes_all"] = df[filter_cols].all(axis=1)

    return df


def run_phase3(df, output_csv=None):
    # type: (pd.DataFrame, Optional[str]) -> pd.DataFrame
    """Execute Phase 3: Molecular filtering."""
    logger.info("=" * 60)
    logger.info("PHASE 3 — Molecular Filtering")
    logger.info("=" * 60)

    df = apply_filters(df)

    n_total = len(df)
    n_pass = df["passes_all"].sum()
    logger.info("Input: %d molecules", n_total)

    for name in list(FILTERS.keys()) + ["pains"]:
        col = "pass_%s" % name
        n = df[col].sum()
        logger.info("  %-15s  %d / %d pass (%.0f%%)", name, n, n_total, 100 * n / max(n_total, 1))

    logger.info("Overall: %d / %d pass all filters (%.0f%%)",
                n_pass, n_total, 100 * n_pass / max(n_total, 1))

    df_filtered = df[df["passes_all"]].copy().reset_index(drop=True)

    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        df.to_csv(output_csv, index=False)
        logger.info("Saved (all with pass flags) -> %s", output_csv)

    return df_filtered
