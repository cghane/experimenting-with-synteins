"""
Shared amino acid property tables and utility functions.

All scales from peer-reviewed literature:
- KD hydrophobicity: Kyte & Doolittle, J Mol Biol 1982
- Helix propensity: Pace & Scholtz, Biophys J 1998
- Aggregation: PASTA scale, Trovato et al., PLoS Comput Biol 2006
- Protease cleavage: MEROPS database heuristics
"""

import math
import logging
from typing import List, Dict, Tuple, Optional

logger = logging.getLogger(__name__)

# ── Amino acid property scales ────────────────────────────────────────────────

# Kyte-Doolittle hydrophobicity (normalized, higher = more hydrophobic)
KD_HYDROPHOBICITY: Dict[str, float] = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}

# Net charge at pH 7.4 (simplified Henderson-Hasselbalch)
CHARGE_PH74: Dict[str, float] = {
    'R': +1.0, 'K': +1.0, 'H': +0.1,
    'D': -1.0, 'E': -1.0,
    'A':  0.0, 'C':  0.0, 'F':  0.0, 'G':  0.0, 'I':  0.0,
    'L':  0.0, 'M':  0.0, 'N':  0.0, 'P':  0.0, 'Q':  0.0,
    'S':  0.0, 'T':  0.0, 'V':  0.0, 'W':  0.0, 'Y':  0.0,
}

# N-terminus (+1) and C-terminus (-1) charges included separately
NTERM_CHARGE = +1.0
CTERM_CHARGE = -1.0

# Helix propensity (Pace & Scholtz 1998; A=1.0 reference)
HELIX_PROPENSITY: Dict[str, float] = {
    'A': 1.00, 'R': 0.79, 'N': 0.65, 'D': 0.72, 'C': 0.68,
    'Q': 0.85, 'E': 0.99, 'G': 0.16, 'H': 0.79, 'I': 0.79,
    'L': 0.98, 'K': 0.87, 'M': 0.93, 'F': 0.81, 'P': 0.00,
    'S': 0.72, 'T': 0.72, 'W': 0.81, 'Y': 0.72, 'V': 0.83,
}

# PASTA aggregation propensity (higher = more aggregation-prone)
AGGREGATION_PROPENSITY: Dict[str, float] = {
    'A':  0.06, 'R': -0.89, 'N': -0.82, 'D': -0.72, 'C':  0.84,
    'Q': -0.58, 'E': -0.64, 'G':  0.00, 'H': -0.28, 'I':  1.51,
    'L':  1.03, 'K': -1.14, 'M':  0.64, 'F':  2.14, 'P': -0.16,
    'S': -0.35, 'T': -0.31, 'W':  1.48, 'Y':  0.95, 'V':  1.08,
}

# Molecular weight (Da)
MW: Dict[str, float] = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
    'Q': 146.2, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
    'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
    'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1,
}

# One-letter to three-letter code
AA1TO3: Dict[str, str] = {
    'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU',
    'G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE',
    'P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'
}

AA3TO1: Dict[str, str] = {v: k for k, v in AA1TO3.items()}

# Known trypsin-like protease cleavage sites (C-terminal to K/R, not before P)
TRYPSIN_PATTERN_RE = r'[KR](?!P)'

# Chymotrypsin cleavage: after F, Y, W (not before P)
CHYMO_PATTERN_RE = r'[FYW](?!P)'


# ── Sequence property calculators ────────────────────────────────────────────

def net_charge(seq: str) -> float:
    """Net charge at pH 7.4 including termini."""
    return NTERM_CHARGE + CTERM_CHARGE + sum(CHARGE_PH74.get(aa, 0) for aa in seq)


def molecular_weight(seq: str) -> float:
    """Approximate molecular weight in Da (subtract water for each peptide bond)."""
    return sum(MW.get(aa, 110.0) for aa in seq) - (len(seq) - 1) * 18.02


def mean_hydrophobicity(seq: str) -> float:
    """Mean Kyte-Doolittle hydrophobicity across whole sequence."""
    return sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in seq) / len(seq)


def windowed_hydrophobicity(seq: str, window: int = 5) -> List[float]:
    """KD hydrophobicity averaged over a sliding window."""
    scores = []
    for i in range(len(seq) - window + 1):
        window_seq = seq[i:i + window]
        scores.append(sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in window_seq) / window)
    return scores


def max_hydrophobic_window(seq: str, window: int = 5) -> float:
    """Maximum hydrophobicity in any sliding window (key filter metric)."""
    windows = windowed_hydrophobicity(seq, window)
    return max(windows) if windows else 0.0


def windowed_aggregation(seq: str, window: int = 5) -> List[float]:
    """PASTA aggregation propensity in sliding window."""
    scores = []
    for i in range(len(seq) - window + 1):
        window_seq = seq[i:i + window]
        scores.append(sum(AGGREGATION_PROPENSITY.get(aa, 0) for aa in window_seq) / window)
    return scores


def max_aggregation_window(seq: str, window: int = 5) -> float:
    """Peak aggregation propensity (key filter metric)."""
    windows = windowed_aggregation(seq, window)
    return max(windows) if windows else 0.0


def sequence_entropy(seq: str) -> float:
    """Shannon entropy over amino acid composition (bits)."""
    from collections import Counter
    counts = Counter(seq)
    total = len(seq)
    return -sum((c / total) * math.log2(c / total) for c in counts.values())


def windowed_entropy(seq: str, window: int = 12) -> List[float]:
    """Shannon entropy in sliding window — low entropy = low complexity."""
    scores = []
    for i in range(len(seq) - window + 1):
        scores.append(sequence_entropy(seq[i:i + window]))
    return scores


def min_complexity(seq: str, window: int = 12) -> float:
    """Minimum local Shannon entropy; low values flag repetitive regions."""
    if len(seq) < window:
        return sequence_entropy(seq)
    return min(windowed_entropy(seq, window))


def mean_helix_propensity(seq: str) -> float:
    """Mean helix propensity across sequence (0–1)."""
    return sum(HELIX_PROPENSITY.get(aa, 0.5) for aa in seq) / len(seq)


def count_cysteines(seq: str) -> int:
    return seq.count('C')


def count_prolines(seq: str) -> int:
    return seq.count('P')


def sequence_similarity(seq1: str, seq2: str) -> float:
    """Fraction of identical residues at same positions (global, no gaps)."""
    if len(seq1) != len(seq2):
        return 0.0
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches / len(seq1)


def hamming_distance(seq1: str, seq2: str) -> int:
    """Number of positions where sequences differ."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be same length for Hamming distance")
    return sum(a != b for a, b in zip(seq1, seq2))


def find_protease_sites(seq: str) -> Dict[str, List[int]]:
    """
    Locate protease recognition motifs (positions are 0-indexed).
    Returns dict mapping protease name to list of cleavage positions.
    """
    import re
    sites = {}

    # Trypsin: cleaves after K or R (not before P)
    trypsin = [m.start() for m in re.finditer(r'[KR](?!P)', seq)]
    if trypsin:
        sites['trypsin'] = trypsin

    # Chymotrypsin: cleaves after F, Y, W (not before P)
    chymo = [m.start() for m in re.finditer(r'[FYW](?!P)', seq)]
    if chymo:
        sites['chymotrypsin'] = chymo

    # Asp-N: cleaves before D
    aspn = [m.start() for m in re.finditer(r'(?=D)', seq[1:])]
    if aspn:
        sites['asp_n'] = [p + 1 for p in aspn]

    # Thrombin: LVPR↓GS
    thrombin = [m.start() for m in re.finditer(r'LVRP', seq)]
    if thrombin:
        sites['thrombin'] = thrombin

    return sites


def total_protease_sites(seq: str) -> int:
    """Total count of protease cleavage sites across all enzymes."""
    return sum(len(v) for v in find_protease_sites(seq).values())


def immunogenicity_score(seq: str) -> float:
    """
    Simplified immunogenicity heuristic (0–1 scale, lower = better).
    Penalises:
      - Low sequence entropy (repetitive = T-cell epitopes)
      - Extreme hydrophobic stretches (MHC-II binding)
      - Long runs of a single amino acid
    """
    score = 0.0

    # Low-entropy penalty
    ent = sequence_entropy(seq)
    # Max entropy for 20 AAs = log2(20) ≈ 4.32
    score += max(0.0, (3.0 - ent) / 3.0) * 0.4

    # Hydrophobic stretch penalty (MHC-II binds 9-mer hydrophobic cores)
    hydro_windows = windowed_hydrophobicity(seq, window=9)
    if hydro_windows:
        max_hydro = max(hydro_windows)
        score += max(0.0, (max_hydro - 1.5) / 3.0) * 0.4

    # Long runs penalty
    max_run = _longest_run(seq)
    score += min(1.0, max_run / 5.0) * 0.2

    return min(1.0, score)


def _longest_run(seq: str) -> int:
    """Length of longest single-amino-acid run in sequence."""
    if not seq:
        return 0
    max_run = 1
    current_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


# ── Structural annotation helpers ────────────────────────────────────────────

def identify_helix_staple_sites(seq: str, spacing: int = 4) -> List[Tuple[int, int]]:
    """
    Find (i, i+spacing) pairs suitable for lactam bridge stapling.
    Ideal: E at position i, K at position i+spacing (or vice versa).
    """
    pairs = []
    for i in range(len(seq) - spacing):
        j = i + spacing
        if (seq[i] in 'EDA' and seq[j] in 'KR') or (seq[i] in 'KR' and seq[j] in 'EDA'):
            pairs.append((i, j))
    return pairs


def identify_cyclization_sites(seq: str) -> Dict[str, str]:
    """
    Annotate N-C terminus cyclization potential.
    Head-to-tail cyclization is feasible when:
    - N-terminus: G, A (flexible)
    - C-terminus: G, A, or any (native chemical ligation with Cys)
    """
    suggestions = {}
    if seq[0] in 'GA':
        suggestions['n_terminus'] = f"Position 1 ({seq[0]}): suitable for head-to-tail cyclization"
    if seq[-1] == 'C':
        suggestions['c_terminus'] = "C-terminal Cys: native chemical ligation (NCL) compatible"
    elif seq[-1] in 'GA':
        suggestions['c_terminus'] = f"C-terminal {seq[-1]}: flexible linker for ring closure"
    return suggestions


def identify_d_amino_acid_sites(seq: str) -> List[Dict]:
    """
    Identify positions where D-amino acid substitution would improve
    protease resistance without disrupting binding.
    Targets: exposed Gly (no side chain penalty) or loop residues.
    Heuristic: positions 0 and len-1, and positions adjacent to protease sites.
    """
    suggestions = []
    protease = find_protease_sites(seq)
    cleavage_positions = set()
    for sites in protease.values():
        for pos in sites:
            cleavage_positions.update([pos, pos + 1])

    for pos in sorted(cleavage_positions):
        if 0 <= pos < len(seq):
            aa = seq[pos]
            suggestions.append({
                'position': pos + 1,  # 1-indexed
                'residue': aa,
                'rationale': f"Adjacent to protease cleavage site; D-{AA1TO3.get(aa, aa)} swap disrupts recognition",
            })

    # Also suggest Gly positions as D-Gly has no stereochemical penalty
    for i, aa in enumerate(seq):
        if aa == 'G':
            suggestions.append({
                'position': i + 1,
                'residue': 'G',
                'rationale': "D-Gly substitution: no side chain, minimal structural perturbation",
            })

    return suggestions


# ── PDB utilities ─────────────────────────────────────────────────────────────

def parse_plddt_from_pdb(pdb_string: str) -> List[float]:
    """
    Extract per-residue pLDDT values from ESMFold/AlphaFold PDB output.
    pLDDT is stored in the B-factor (temperature factor) column.
    """
    plddts = []
    seen_residues = set()
    for line in pdb_string.split('\n'):
        if not line.startswith('ATOM'):
            continue
        # Only take CA atoms to get one value per residue
        atom_name = line[12:16].strip()
        if atom_name != 'CA':
            continue
        chain = line[21]
        resnum = int(line[22:26].strip())
        key = (chain, resnum)
        if key in seen_residues:
            continue
        seen_residues.add(key)
        try:
            bfactor = float(line[60:66].strip())
            plddts.append(bfactor)
        except (ValueError, IndexError):
            continue
    return plddts


def calc_contact_map(
    coords1: List[Tuple[float, float, float]],
    coords2: List[Tuple[float, float, float]],
    cutoff: float = 5.0
) -> Tuple[int, float]:
    """
    Count inter-chain contacts within cutoff (Å).
    Returns (n_contacts, mean_min_distance).
    """
    contacts = 0
    min_dists = []
    for c1 in coords1:
        dists = [
            math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)
            for c2 in coords2
        ]
        min_d = min(dists)
        min_dists.append(min_d)
        if min_d <= cutoff:
            contacts += 1
    mean_min = sum(min_dists) / len(min_dists) if min_dists else 999.0
    return contacts, mean_min


def extract_ca_coords(pdb_string: str, chain: Optional[str] = None) -> List[Tuple[float, float, float]]:
    """Extract Cα coordinates from PDB string, optionally filtering by chain."""
    coords = []
    for line in pdb_string.split('\n'):
        if not line.startswith('ATOM'):
            continue
        atom_name = line[12:16].strip()
        if atom_name != 'CA':
            continue
        if chain and line[21] != chain:
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append((x, y, z))
        except (ValueError, IndexError):
            continue
    return coords
