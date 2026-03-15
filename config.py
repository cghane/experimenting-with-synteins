"""
Central configuration for the AI-Guided Mini-Protein Binder Design Pipeline.
Target: IL-6 (Interleukin-6) | Scaffold: Affibody Z-domain
"""

import os

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(BASE_DIR, "data")
RESULTS_DIR= os.path.join(BASE_DIR, "results")
REPORT_DIR = os.path.join(BASE_DIR, "report")

# ── Target ───────────────────────────────────────────────────────────────────
TARGET_NAME   = "IL-6"
TARGET_PDB_ID = "4CNI"          # IL-6 / IL-6Rα / gp130 ternary complex (2.4Å)
RECEPTOR_CHAIN= "A"             # IL-6 chain in 4CNI
RECEPTOR_PDB  = os.path.join(DATA_DIR, "receptor.pdb")
INTERFACE_JSON= os.path.join(DATA_DIR, "interface_residues.json")

# ── Scaffold ─────────────────────────────────────────────────────────────────
# Affibody Z-domain (Protein A, 55 aa). Clinically validated scaffold.
# PDB refs: 2B97, 2SPZ. Used in izokibep (anti-IL-17A, Phase 2).
AFFIBODY_SCAFFOLD = "VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDA"
SCAFFOLD_NAME     = "Z-domain (Affibody)"

# 13 variable positions on the binding face (0-indexed), helices 1 & 2.
# Nygren lab standard: positions 9–18 (H1 face) + 24–28, 32 (H2 face).
VARIABLE_POSITIONS = [8, 9, 10, 12, 13, 15, 17, 23, 24, 25, 26, 27, 31]

# Amino acids allowed at variable positions (excludes C, P, G for stability)
ALLOWED_AAS = list("ADEFHIKLMNQRSTVWY")

# ── Candidate generation ─────────────────────────────────────────────────────
N_CANDIDATES_GENERATE = 1000
RANDOM_SEED           = 42

# ── Developability filters ───────────────────────────────────────────────────
MAX_HYDROPHOBICITY_WINDOW = 2.5    # KD scale, 5-aa window mean
MIN_NET_CHARGE            = -3
MAX_NET_CHARGE            = 6
MAX_CYSTEINES             = 1
MAX_AGGREGATION_WINDOW    = 2.0    # PASTA scale, 5-aa window mean
MIN_COMPLEXITY_ENTROPY    = 2.0    # Shannon entropy, 12-aa window

# ── Phase 4: Structure prediction ────────────────────────────────────────────
ESMATLAS_URL  = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMATLAS_TIMEOUT = 60             # seconds per request
N_SEEDS          = 3             # predictions per candidate
CONTACT_CUTOFF_A = 5.0           # Å — interface contact distance

# ── Phase 5: Scoring weights ─────────────────────────────────────────────────
SCORE_WEIGHTS = {
    "interface_plddt":  0.35,
    "bsa_score":        0.30,
    "hbond_score":      0.15,
    "aggregation_pen": -0.10,
    "protease_pen":    -0.10,
}

# ── Phase 6: Optimization ────────────────────────────────────────────────────
OPTIM_ROUNDS     = 6
OPTIM_CANDIDATES = 10
OPTIM_MUTATIONS  = 20   # candidates per round per binder
OPTIM_T_START    = 1.0
OPTIM_T_END      = 0.05

# ── Phase 7: Synteins thresholds ─────────────────────────────────────────────
MAX_PROTEASE_SITES   = 3
MAX_IMMUNOGEN_SCORE  = 0.4
HELIX_STAPLE_SPACING = 4     # i, i+4 for lactam bridge

# ── Logging ──────────────────────────────────────────────────────────────────
LOG_LEVEL = "INFO"
