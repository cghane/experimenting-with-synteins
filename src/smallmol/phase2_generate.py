"""
Phase 2 — Molecule Generation
Fragment-based enumeration: attach decoration fragments to pharmacologically
relevant scaffolds to build a diverse candidate library.
"""

import logging
import os
import random
from typing import Dict, List, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from src.smallmol.utils import smiles_to_mol, compute_descriptors, mol_to_canonical_smiles
import smallmol_config as cfg

logger = logging.getLogger(__name__)

# Enumerate by combining scaffolds + substituent SMILES
# using direct SMILES string construction at open positions

# Known IL-6 inhibitor-inspired decorated scaffolds
DECORATED_TEMPLATES = [
    # Indole derivatives
    "c1ccc2[nH]c({R1})c({R2})c2c1",
    "c1cc({R1})c2[nH]c({R2})cc2c1",
    # Phenylpyridine derivatives
    "c1cc({R1})c(-c2cc({R2})ccn2)cc1",
    "c1ccc(-c2cccc({R1})n2)c({R2})c1",
    # Quinazoline derivatives
    "c1cc({R2})c2nc({R1})ncc2c1",
    "c1cnc2cc({R1})c({R2})cc2n1",
    # Benzimidazole derivatives
    "c1cc({R1})c2[nH]c({R2})nc2c1",
    "c1ccc2[nH]c({R1})nc2c({R2})1",
    # Pyrazole derivatives
    "c1c({R1})n({R2})nc1",
    "c1c({R1})c({R2})n[nH]1",
]

# Fragment SMILES for R-group substitution
POLAR_FRAGMENTS = [
    "O", "N", "C(=O)N", "C(=O)O", "S(=O)(=O)N", "C(N)=O",
    "C1COCCN1", "C1CNCCN1", "c1nn[nH]n1", "C(=O)NC",
    "NC(=O)C", "OC", "NCC",
]

HYDROPHOBIC_FRAGMENTS = [
    "C", "C(F)(F)F", "F", "Cl", "c1ccccc1", "C1CC1", "CC",
    "C(C)C", "CC(C)C",
]

ALL_FRAGMENTS = POLAR_FRAGMENTS + HYDROPHOBIC_FRAGMENTS


def _build_molecule(template, r1_frag, r2_frag):
    # type: (str, str, str) -> Optional[str]
    """Build a molecule by substituting R-groups into a template."""
    smiles = template.replace("{R1}", r1_frag).replace("{R2}", r2_frag)
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    return mol_to_canonical_smiles(mol)


def _get_scaffold_name(template):
    # type: (str) -> str
    """Extract scaffold name from template."""
    if "nH]c" in template and "c2c1" in template:
        return "indole"
    elif "ccn2" in template or "ncc2" in template:
        return "phenylpyridine"
    elif "nc(" in template.lower() and "ncc2" in template:
        return "quinazoline"
    elif "[nH]c" in template and "nc2c" in template:
        return "benzimidazole"
    elif "n[nH]" in template or "n({R2})nc" in template:
        return "pyrazole"
    # Fallback: check for specific patterns
    for name in ["indole", "phenylpyridine", "quinazoline", "benzimidazole", "pyrazole"]:
        if name in template.lower():
            return name
    return "other"


def generate_library(pocket_profile, n_molecules, seed, biased_fraction):
    # type: (Dict, int, int, float) -> pd.DataFrame
    """
    Generate a library of small molecules by combining scaffolds with fragments.
    Biased molecules use fragments that complement the pocket pharmacophore.
    """
    rng = random.Random(seed)
    n_biased = int(n_molecules * biased_fraction)
    n_random = n_molecules - n_biased

    # Bias fragment selection toward pocket complementarity
    # Pocket is electrostatic/polar dominant -> prefer polar fragments
    biased_weights_r1 = []
    for f in ALL_FRAGMENTS:
        if f in POLAR_FRAGMENTS:
            biased_weights_r1.append(3.0)  # 3x weight for polar
        else:
            biased_weights_r1.append(1.0)

    biased_weights_r2 = list(biased_weights_r1)  # same weighting

    molecules = []
    seen_smiles = set()
    attempts = 0
    max_attempts = n_molecules * 20

    # Generate biased molecules
    while len(molecules) < n_biased and attempts < max_attempts:
        attempts += 1
        template = rng.choice(DECORATED_TEMPLATES)
        r1 = rng.choices(ALL_FRAGMENTS, weights=biased_weights_r1, k=1)[0]
        r2 = rng.choices(ALL_FRAGMENTS, weights=biased_weights_r2, k=1)[0]
        smi = _build_molecule(template, r1, r2)
        if smi and smi not in seen_smiles:
            seen_smiles.add(smi)
            scaffold = _get_scaffold_name(template)
            molecules.append({
                "smiles": smi,
                "scaffold": scaffold,
                "r1_fragment": r1,
                "r2_fragment": r2,
                "is_biased": True,
            })

    # Generate random molecules
    attempts = 0
    while len(molecules) < n_biased + n_random and attempts < max_attempts:
        attempts += 1
        template = rng.choice(DECORATED_TEMPLATES)
        r1 = rng.choice(ALL_FRAGMENTS)
        r2 = rng.choice(ALL_FRAGMENTS)
        smi = _build_molecule(template, r1, r2)
        if smi and smi not in seen_smiles:
            seen_smiles.add(smi)
            scaffold = _get_scaffold_name(template)
            molecules.append({
                "smiles": smi,
                "scaffold": scaffold,
                "r1_fragment": r1,
                "r2_fragment": r2,
                "is_biased": False,
            })

    df = pd.DataFrame(molecules)
    df["id"] = ["SMOL-%04d" % (i + 1) for i in range(len(df))]

    # Compute descriptors
    desc_rows = []
    for smi in df["smiles"]:
        mol = smiles_to_mol(smi)
        if mol:
            desc_rows.append(compute_descriptors(mol))
        else:
            desc_rows.append({})

    desc_df = pd.DataFrame(desc_rows)
    df = pd.concat([df, desc_df], axis=1)

    return df


def run_phase2(pocket_profile, output_csv=None):
    # type: (Dict, Optional[str]) -> pd.DataFrame
    """Execute Phase 2: Molecule generation."""
    logger.info("=" * 60)
    logger.info("PHASE 2 — Molecule Generation")
    logger.info("=" * 60)

    df = generate_library(
        pocket_profile,
        n_molecules=cfg.N_MOLECULES_GENERATE,
        seed=cfg.RANDOM_SEED,
        biased_fraction=cfg.BIASED_FRACTION,
    )

    n_biased = df["is_biased"].sum()
    n_random = len(df) - n_biased
    logger.info("Generated %d molecules (%d biased, %d random)", len(df), n_biased, n_random)
    logger.info("Unique scaffolds: %s", sorted(df["scaffold"].unique()))
    logger.info("MW range: %.0f – %.0f", df["mw"].min(), df["mw"].max())

    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        df.to_csv(output_csv, index=False)
        logger.info("Saved -> %s", output_csv)

    return df
