"""
Central configuration for the Small Molecule Drug Design Pipeline.
Target: IL-6 Site II binding pocket | Modality: Small molecules
"""

import os

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results", "smallmol")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")

# ── Target ───────────────────────────────────────────────────────────────────
TARGET_NAME    = "IL-6"
TARGET_PDB_ID  = "4CNI"
RECEPTOR_CHAIN = "A"
RECEPTOR_PDB   = os.path.join(DATA_DIR, "receptor.pdb")
INTERFACE_JSON = os.path.join(DATA_DIR, "interface_residues.json")

# Pocket pharmacophore profile (from interface_residues.json analysis)
POCKET_PROFILE = {
    "n_hydrophobic_hotspots": 3,   # F73, I82, W157
    "n_charged_hotspots": 7,       # K77, E110, R113, K120, R160, R168, K171
    "n_hbond_hotspots": 4,         # S38, Q75, Q118, Q175
    "dominant_character": "electrostatic/polar",
    "net_charge": "+3 to +4 at pH 7.4",
}

# ── Scaffolds ────────────────────────────────────────────────────────────────
# Pharmacologically relevant scaffolds for cytokine inhibitors
SCAFFOLD_SMILES = {
    "indole":          "c1ccc2[nH]ccc2c1",
    "phenylpyridine":  "c1ccc(-c2ccccn2)cc1",
    "quinazoline":     "c1cnc2ccccc2n1",
    "benzimidazole":   "c1ccc2[nH]cnc2c1",
    "pyrazole":        "c1cn[nH]c1",
}

# Decoration fragments — functional groups attached to scaffolds
FRAGMENTS = {
    # H-bond donors/acceptors (complement pocket polar residues)
    "hydroxyl":      "[OH]",
    "amine":         "[NH2]",
    "amide":         "C(=O)N",
    "carboxamide":   "C(=O)[NH2]",
    "methoxy":       "OC",
    "morpholine":    "C1COCCN1",
    "piperazine":    "C1CNCCN1",
    "sulfonamide":   "S(=O)(=O)N",
    # Hydrophobic (complement F73, I82, W157)
    "methyl":        "C",
    "trifluoromethyl": "C(F)(F)F",
    "fluorine":      "F",
    "chlorine":      "Cl",
    "phenyl":        "c1ccccc1",
    "cyclopropyl":   "C1CC1",
    # Charged (complement K77, R113, R168)
    "carboxylic":    "C(=O)O",
    "tetrazole":     "c1nn[nH]n1",
    "acylsulfonamide": "S(=O)(=O)NC(=O)",
}

# ── Generation ───────────────────────────────────────────────────────────────
N_MOLECULES_GENERATE = 500
RANDOM_SEED          = 42
BIASED_FRACTION      = 0.70   # 70% biased toward pocket complementarity

# ── Lipinski / ADMET thresholds ──────────────────────────────────────────────
LIPINSKI_MW_MAX       = 500.0
LIPINSKI_LOGP_MAX     = 5.0
LIPINSKI_HBD_MAX      = 5
LIPINSKI_HBA_MAX      = 10
VEBER_TPSA_MAX        = 140.0
VEBER_ROTBONDS_MAX    = 10
MIN_RINGS             = 1
MIN_MW                = 150.0
AGGREGATION_LOGP_MIN  = 4.0
AGGREGATION_MW_MIN    = 350.0

# ── Scoring weights ──────────────────────────────────────────────────────────
SCORE_WEIGHTS = {
    "qed":           0.25,
    "sa_score":      0.20,   # inverted: lower SA = better
    "binding_proxy": 0.30,
    "logp_penalty": -0.10,
    "tpsa_bonus":    0.15,
}

# Normalization ranges for scoring
NORM_RANGES = {
    "qed":           (0.0, 1.0),
    "sa_inv":        (0.0, 1.0),
    "binding_proxy": (0.0, 1.0),
    "logp_dev":      (0.0, 4.0),     # |logP - 2.5|
    "tpsa_score":    (0.0, 1.0),
}

# ── Optimization ─────────────────────────────────────────────────────────────
OPTIM_ROUNDS           = 8
OPTIM_TOP_N            = 10
OPTIM_MUTATIONS_PER_ROUND = 15
OPTIM_T_START          = 1.0
OPTIM_T_END            = 0.05

# ── Lead analysis ────────────────────────────────────────────────────────────
LEAD_TOP_N = 5

# ── Logging ──────────────────────────────────────────────────────────────────
LOG_LEVEL = "INFO"
