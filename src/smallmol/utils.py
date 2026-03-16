"""
Shared RDKit utility functions for the small molecule pipeline.
"""

import logging
from typing import Optional, Dict, List

from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, AllChem, DataStructs

logger = logging.getLogger(__name__)


def smiles_to_mol(smiles):
    # type: (str) -> Optional[Chem.Mol]
    """Safely parse SMILES to RDKit Mol, returns None on failure."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def compute_descriptors(mol):
    # type: (Chem.Mol) -> Dict[str, float]
    """Compute standard molecular descriptors."""
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "tpsa": Descriptors.TPSA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "ring_count": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "num_heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
    }


def compute_qed(mol):
    # type: (Chem.Mol) -> float
    """Quantitative Estimate of Drug-likeness (0–1)."""
    try:
        return QED.qed(mol)
    except Exception:
        return 0.0


def compute_sa_score(mol):
    # type: (Chem.Mol) -> float
    """
    Simplified synthetic accessibility score (1–10, lower = easier).
    Based on fragment complexity heuristics.
    """
    try:
        ring_count = rdMolDescriptors.CalcNumRings(mol)
        stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        hetero_ratio = rdMolDescriptors.CalcNumHeteroatoms(mol) / max(mol.GetNumHeavyAtoms(), 1)
        spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        bridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)

        score = 2.0
        score += ring_count * 0.5
        score += stereo * 0.8
        score += spiro * 1.0
        score += bridgehead * 1.2
        score += max(0, hetero_ratio - 0.3) * 3.0
        score += max(0, mol.GetNumHeavyAtoms() - 30) * 0.1

        return min(max(score, 1.0), 10.0)
    except Exception:
        return 5.0


def compute_morgan_fp(mol, radius=2, n_bits=2048):
    # type: (Chem.Mol, int, int) -> DataStructs.ExplicitBitVect
    """Morgan (circular) fingerprint."""
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def tanimoto_similarity(mol1, mol2):
    # type: (Chem.Mol, Chem.Mol) -> float
    """Tanimoto similarity between two molecules (Morgan FP)."""
    fp1 = compute_morgan_fp(mol1)
    fp2 = compute_morgan_fp(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def count_pharmacophore_features(mol):
    # type: (Chem.Mol) -> Dict[str, int]
    """Count pharmacophore-relevant features."""
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)

    # Hydrophobic: aromatic rings + aliphatic carbons in rings
    aromatic = rdMolDescriptors.CalcNumAromaticRings(mol)
    aliphatic_rings = rdMolDescriptors.CalcNumAliphaticRings(mol)
    hydrophobic = aromatic + aliphatic_rings

    # Charged groups: carboxylic acids (neg) + basic amines (pos)
    n_neg = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]"))) if Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]") else 0
    n_pos = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]"))) if Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]") else 0

    # Sulfonamide as negative-charge mimic
    n_sulfonamide = len(mol.GetSubstructMatches(Chem.MolFromSmarts("S(=O)(=O)N"))) if Chem.MolFromSmarts("S(=O)(=O)N") else 0
    # Tetrazole as carboxylic acid bioisostere
    n_tetrazole = len(mol.GetSubstructMatches(Chem.MolFromSmarts("c1nn[nH]n1"))) if Chem.MolFromSmarts("c1nn[nH]n1") else 0

    charged = n_neg + n_pos + n_sulfonamide + n_tetrazole

    return {
        "hbd": hbd,
        "hba": hba,
        "hydrophobic": hydrophobic,
        "charged": charged,
        "n_positive": n_pos,
        "n_negative": n_neg + n_sulfonamide + n_tetrazole,
    }


def mol_to_canonical_smiles(mol):
    # type: (Chem.Mol) -> Optional[str]
    """Return canonical SMILES or None."""
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return None
