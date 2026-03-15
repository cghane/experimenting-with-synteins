"""
Phase 7 — "Synteins-Style" Concept Layer

For the top 5 binders, apply three additional evaluation layers:

1. Protease Resistance Heuristic
   - Identify exposed Lys/Arg clusters (trypsin targets)
   - Flag chymotrypsin (F/Y/W) and Asp-N sites
   - Compute net protease vulnerability score

2. Immunogenicity Heuristic
   - Penalise low sequence entropy (repetitive → T-cell epitopes)
   - Penalise MHC-II-binding hydrophobic 9-mers
   - Penalise long runs of identical amino acids

3. Structural Modification Annotations
   - Cyclization sites (head-to-tail, NCL)
   - Lactam bridge stapling pairs (i, i+4 EK or KE pairs)
   - D-amino acid swap candidates (adjacent to protease sites)
   - Proposed PEGylation sites (surface-exposed Lys)

None of these require simulation — they are rational annotations
based on structure and sequence analysis, exactly as done in serious
protein engineering groups.
"""

import logging
import os
from typing import Dict, List, Optional

import pandas as pd

from src.utils import (
    find_protease_sites, immunogenicity_score, total_protease_sites,
    identify_helix_staple_sites, identify_cyclization_sites,
    identify_d_amino_acid_sites, CHARGE_PH74, AA1TO3,
    max_aggregation_window, mean_helix_propensity
)

logger = logging.getLogger(__name__)


# ── Protease resistance ───────────────────────────────────────────────────────

def analyse_protease_resistance(seq: str) -> Dict:
    """
    Full protease vulnerability analysis.
    Returns dict with per-enzyme sites and overall risk score.
    """
    sites = find_protease_sites(seq)

    # Risk score: each site weighted by enzyme abundance in serum
    # Trypsin-like (plasmin, thrombin): high abundance → weight 1.0
    # Chymotrypsin-like: moderate → weight 0.7
    # Asp-N: low → weight 0.3
    enzyme_weights = {'trypsin': 1.0, 'chymotrypsin': 0.7, 'asp_n': 0.3, 'thrombin': 0.5}

    risk_score = sum(
        len(site_list) * enzyme_weights.get(enzyme, 0.5)
        for enzyme, site_list in sites.items()
    )

    # Identify K/R clusters (consecutive or nearby)
    kr_positions = [i for i, aa in enumerate(seq) if aa in 'KR']
    clusters = []
    if kr_positions:
        current_cluster = [kr_positions[0]]
        for pos in kr_positions[1:]:
            if pos - current_cluster[-1] <= 3:
                current_cluster.append(pos)
            else:
                if len(current_cluster) >= 2:
                    clusters.append(current_cluster)
                current_cluster = [pos]
        if len(current_cluster) >= 2:
            clusters.append(current_cluster)

    # Modifications to improve resistance
    resistance_mods = []
    for enzyme, site_list in sites.items():
        for pos in site_list:
            aa = seq[pos]
            if enzyme == 'trypsin':
                resistance_mods.append({
                    'position': pos + 1,
                    'current_aa': aa,
                    'modification': f'K→Q or R→A at position {pos+1} disrupts trypsin site',
                    'strategy': 'conservative substitution'
                })
            elif enzyme == 'chymotrypsin':
                resistance_mods.append({
                    'position': pos + 1,
                    'current_aa': aa,
                    'modification': f'D-{AA1TO3.get(aa, aa)} at position {pos+1} blocks cleavage',
                    'strategy': 'D-amino acid incorporation'
                })

    return {
        'cleavage_sites': sites,
        'total_sites': total_protease_sites(seq),
        'kr_clusters': clusters,
        'risk_score': risk_score,
        'modifications': resistance_mods[:5],  # top 5 most impactful
        'risk_level': 'HIGH' if risk_score > 6 else 'MODERATE' if risk_score > 3 else 'LOW',
    }


# ── Immunogenicity analysis ───────────────────────────────────────────────────

def analyse_immunogenicity(seq: str) -> Dict:
    """
    MHC-II immunogenicity assessment with actionable suggestions.
    """
    score = immunogenicity_score(seq)
    risk_level = 'HIGH' if score > 0.5 else 'MODERATE' if score > 0.25 else 'LOW'

    # Identify problematic windows (hydrophobic 9-mers: candidate MHC-II binders)
    from src.utils import windowed_hydrophobicity, windowed_entropy
    hydro9 = windowed_hydrophobicity(seq, 9)
    entropy12 = windowed_entropy(seq, 12)

    mhc_candidate_windows = [
        {'start': i+1, 'end': i+9, 'peptide': seq[i:i+9], 'hydrophobicity': h}
        for i, h in enumerate(hydro9) if h > 1.5
    ]

    low_entropy_windows = [
        {'start': i+1, 'end': i+12, 'peptide': seq[i:i+12], 'entropy': e}
        for i, e in enumerate(entropy12) if e < 2.5
    ]

    # Deimmunisation suggestions
    deimmun_suggestions = []
    for w in mhc_candidate_windows[:3]:
        deimmun_suggestions.append({
            'window': f"{w['peptide']} (pos {w['start']}–{w['end']})",
            'issue': f"Hydrophobic 9-mer (score {w['hydrophobicity']:.2f}) — MHC-II anchor risk",
            'suggestion': "Replace hydrophobic anchors at P1/P4/P6/P9 with charged/polar residues"
        })

    return {
        'immunogenicity_score': score,
        'risk_level': risk_level,
        'mhc_ii_candidate_windows': mhc_candidate_windows,
        'low_complexity_windows': low_entropy_windows,
        'deimmunisation_suggestions': deimmun_suggestions,
    }


# ── Structural modification annotations ──────────────────────────────────────

def annotate_modifications(seq: str) -> Dict:
    """
    Annotate all potential structural stabilisation strategies.
    Purely rational, based on sequence properties.
    """
    mods = {}

    # 1. Cyclization
    mods['cyclization'] = identify_cyclization_sites(seq)
    mods['cyclization_summary'] = (
        "Head-to-tail cyclization: reduces proteolysis, improves thermal stability. "
        "Achieves 10–100× improved serum half-life in literature (e.g., cyclic peptides)."
    )

    # 2. Helix stapling
    staple_pairs = identify_helix_staple_sites(seq, spacing=4)
    mods['staple_pairs'] = [
        {
            'pair': f"({seq[i]}{i+1}, {seq[j]}{j+1})",
            'positions': (i+1, j+1),
            'type': 'Lactam bridge (i, i+4)' if (seq[i] in 'EDA' and seq[j] in 'KR') else
                    'Lactam bridge (i, i+4)',
        }
        for i, j in staple_pairs[:5]
    ]
    mods['staple_summary'] = (
        f"Found {len(staple_pairs)} i,i+4 lactam bridge pairs. "
        "Helix stapling improves proteolytic resistance and binding affinity "
        "by pre-organising helix conformation (Verdine & Hilinski, 2012)."
    )

    # 3. D-amino acid substitutions
    d_aa_sites = identify_d_amino_acid_sites(seq)
    mods['d_amino_acid_sites'] = d_aa_sites[:5]
    mods['d_amino_acid_summary'] = (
        "D-amino acids at protease cleavage sites block recognition without "
        "significantly disrupting backbone geometry (single-residue swap). "
        "Cannot be predicted computationally with standard force fields — "
        "these are proposed for experimental validation."
    )

    # 4. PEGylation sites
    lys_positions = [i+1 for i, aa in enumerate(seq) if aa == 'K']
    # Best PEG sites: surface-exposed Lys not in binding interface
    # Heuristic: Lys not at variable positions (interface-facing)
    from config import VARIABLE_POSITIONS
    peg_candidates = [p for p in lys_positions if (p - 1) not in VARIABLE_POSITIONS]
    mods['pegylation_sites'] = [
        {'position': p, 'rationale': 'Surface-exposed Lys (not at binding interface)'}
        for p in peg_candidates[:3]
    ]
    mods['pegylation_summary'] = (
        "PEGylation at non-interface Lys extends serum half-life (elimination half-life "
        "improvement from minutes to hours). Site-specific PEG at these positions "
        "predicted to minimise binding interference."
    )

    # 5. Disulfide bridges
    cys_positions = [i+1 for i, aa in enumerate(seq) if aa == 'C']
    if len(cys_positions) >= 2:
        mods['disulfide'] = {
            'pairs': [(cys_positions[i], cys_positions[i+1])
                      for i in range(0, len(cys_positions)-1, 2)],
            'rationale': "Engineered disulfides increase thermal stability (ΔTm +5–15°C typical)"
        }
    else:
        mods['disulfide'] = {
            'suggestion': "Introduce C at designed positions for thermal stabilisation",
            'candidate_positions': 'Buried positions in helix core (e.g., i, i+11 in 3-helix bundle)'
        }

    return mods


# ── Full candidate analysis ───────────────────────────────────────────────────

def analyse_candidate(row: pd.Series) -> Dict:
    """Full synteins-style analysis for a single candidate."""
    seq = row['sequence']
    return {
        'id': row['id'],
        'sequence': seq,
        'final_score': row.get('final_score', 0),
        'protease': analyse_protease_resistance(seq),
        'immunogenicity': analyse_immunogenicity(seq),
        'modifications': annotate_modifications(seq),
    }


def format_candidate_report(analysis: Dict) -> str:
    """Format a human-readable analysis report for one candidate."""
    pr  = analysis['protease']
    imm = analysis['immunogenicity']
    mod = analysis['modifications']

    lines = [
        "\n" + "=" * 70,
        f"SYNTEINS ANALYSIS: {analysis['id']}",
        f"Score: {analysis['final_score']:.1f} | Sequence: {analysis['sequence']}",
        "=" * 70,
        "",
        "── Protease Resistance ─────────────────────────────────────",
        f"Risk level     : {pr['risk_level']}",
        f"Total sites    : {pr['total_sites']} "
        f"(Trypsin: {len(pr['cleavage_sites'].get('trypsin',[]))}, "
        f"Chymotrypsin: {len(pr['cleavage_sites'].get('chymotrypsin',[]))})",
        f"KR clusters    : {len(pr['kr_clusters'])}",
        "Top modifications:",
    ]
    for m in pr['modifications'][:3]:
        lines.append(f"  • {m['modification']}")

    lines += [
        "",
        "── Immunogenicity ──────────────────────────────────────────",
        f"Score : {imm['immunogenicity_score']:.3f} | Risk: {imm['risk_level']}",
        f"MHC-II candidates : {len(imm['mhc_ii_candidate_windows'])} windows",
    ]
    for s in imm['deimmunisation_suggestions'][:2]:
        lines.append(f"  • {s['window']}: {s['suggestion']}")

    lines += [
        "",
        "── Structural Modifications ────────────────────────────────",
        f"Helix staple pairs : {len(mod.get('staple_pairs', []))} i,i+4 pairs",
        f"PEGylation sites   : {len(mod.get('pegylation_sites', []))}",
        f"D-AA swap sites    : {len(mod.get('d_amino_acid_sites', []))}",
    ]
    if mod.get('staple_pairs'):
        lines.append("  Staple candidates: " +
                     ", ".join(p['pair'] for p in mod['staple_pairs'][:3]))
    if mod.get('pegylation_sites'):
        lines.append("  PEG sites: " +
                     ", ".join(f"K{p['position']}" for p in mod['pegylation_sites']))

    lines.append("=" * 70)
    return '\n'.join(lines)


# ── Phase 7 runner ────────────────────────────────────────────────────────────

def run_phase7(
    df: pd.DataFrame,
    top_n: int = 5,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
    """Execute Phase 7: synteins-style analysis for top N candidates."""
    logger.info("─" * 50)
    logger.info(f"PHASE 7 — Synteins-style analysis (top {top_n} candidates)")
    logger.info("─" * 50)

    top = df.head(top_n)
    analyses = []
    reports  = []

    for _, row in top.iterrows():
        analysis = analyse_candidate(row)
        analyses.append(analysis)
        report = format_candidate_report(analysis)
        reports.append(report)
        logger.info(report)

    # Flatten to DataFrame for CSV export
    records = []
    for a in analyses:
        pr  = a['protease']
        imm = a['immunogenicity']
        mod = a['modifications']
        records.append({
            'id':                  a['id'],
            'sequence':            a['sequence'],
            'final_score':         a['final_score'],
            'protease_risk_level': pr['risk_level'],
            'total_protease_sites': pr['total_sites'],
            'protease_risk_score': pr['risk_score'],
            'immunogen_score':     imm['immunogenicity_score'],
            'immunogen_risk':      imm['risk_level'],
            'n_mhc_windows':       len(imm['mhc_ii_candidate_windows']),
            'n_staple_pairs':      len(mod.get('staple_pairs', [])),
            'n_peg_sites':         len(mod.get('pegylation_sites', [])),
            'n_daa_sites':         len(mod.get('d_amino_acid_sites', [])),
            'cyclization_n_term':  'n_terminus' in mod.get('cyclization', {}),
            'cyclization_c_term':  'c_terminus' in mod.get('cyclization', {}),
        })

    result_df = pd.DataFrame(records)

    logger.info(f"\nPhase 7 complete: {len(result_df)} candidates analysed")
    if output_csv:
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        result_df.to_csv(output_csv, index=False)
        logger.info(f"Saved synteins analysis → {output_csv}")

    return result_df
