"""
Phase 1 — Target Selection & Interface Definition

Downloads the IL-6 PDB structure, validates the binding site,
and produces a short biophysical summary of the target interface.
"""

import os
import json
import logging
import requests
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


# ── PDB download ─────────────────────────────────────────────────────────────

def fetch_pdb(pdb_id: str, output_path: str, chain: Optional[str] = None) -> bool:
    """
    Download a PDB file from RCSB. Optionally filter to a single chain.
    Returns True on success.
    """
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    logger.info(f"Fetching PDB {pdb_id} from RCSB …")
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        content = resp.text
    except requests.RequestException as e:
        logger.error(f"Failed to fetch {pdb_id}: {e}")
        return False

    if chain:
        content = _filter_chain(content, chain)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(content)

    logger.info(f"Saved {pdb_id} → {output_path} ({len(content):,} bytes)")
    return True


def _filter_chain(pdb_content: str, chain: str) -> str:
    """Keep only ATOM/HETATM records for the specified chain."""
    lines = []
    for line in pdb_content.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            if line[21] == chain:
                lines.append(line)
        elif line.startswith(('HEADER', 'TITLE', 'REMARK', 'SEQRES', 'END')):
            lines.append(line)
    return '\n'.join(lines)


# ── Interface analysis ────────────────────────────────────────────────────────

def load_interface_residues(json_path: str) -> Dict:
    """Load interface residue definitions from JSON file."""
    with open(json_path) as f:
        return json.load(f)


def extract_receptor_sequence(pdb_path: str, chain: str = 'A') -> str:
    """
    Parse receptor sequence from PDB ATOM records (CA atoms).
    Returns one-letter code string.
    """
    from src.utils import AA3TO1
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            if line[21] != chain:
                continue
            resnum = int(line[22:26].strip())
            resname = line[17:20].strip()
            residues[resnum] = AA3TO1.get(resname, 'X')

    return ''.join(residues[k] for k in sorted(residues))


def summarise_interface(interface_data: Dict) -> str:
    """
    Return a formatted string summarising the target interface properties.
    """
    site = interface_data.get('site_ii_residues', {})
    residues = site.get('residues', [])
    design = interface_data.get('design_rationale', {})
    surface = design.get('target_surface_properties', {})

    hydrophobic = [r['one_letter'] + str(r['resnum'])
                   for r in residues if r['role'] == 'hydrophobic anchor' or 'hydrophobic' in r['role']]
    electrostatic = [r['one_letter'] + str(r['resnum'])
                     for r in residues if 'electrostatic' in r['role'] or 'charge' in r['role']]
    hbond = [r['one_letter'] + str(r['resnum'])
             for r in residues if 'H-bond' in r['role'] or 'polar' in r['role']]

    lines = [
        "=" * 60,
        f"Target: {interface_data['target']} (PDB: {interface_data['pdb_id']})",
        f"Binding site: {interface_data['binding_site']}",
        "=" * 60,
        f"Interface residues: {len(residues)}",
        f"  Hydrophobic anchors : {', '.join(hydrophobic)}",
        f"  Electrostatic       : {', '.join(electrostatic)}",
        f"  H-bond formers      : {', '.join(hbond)}",
        "",
        "Surface properties (Site II):",
        f"  Net charge          : {surface.get('net_charge_site_ii', 'N/A')}",
        f"  Hydrophobic patch   : {surface.get('hydrophobic_patch', 'N/A')}",
        f"  Polar shell         : {surface.get('polar_shell', 'N/A')}",
        "=" * 60,
    ]
    return '\n'.join(lines)


def run_phase1(receptor_pdb_path: str, interface_json_path: str,
               pdb_id: str = "4CNI", chain: str = "A") -> Dict:
    """
    Execute Phase 1: download receptor PDB and summarise interface.
    Returns a results dict for pipeline integration.
    """
    logger.info("─" * 50)
    logger.info("PHASE 1 — Target selection & interface definition")
    logger.info("─" * 50)

    # Download PDB if not cached
    if not os.path.exists(receptor_pdb_path):
        success = fetch_pdb(pdb_id, receptor_pdb_path, chain=chain)
        if not success:
            logger.warning("PDB download failed — some downstream steps may use cached/mock data.")
    else:
        logger.info(f"PDB already cached at {receptor_pdb_path}")

    # Load interface definitions
    interface_data = load_interface_residues(interface_json_path)
    summary = summarise_interface(interface_data)
    logger.info("\n" + summary)

    # Try to parse receptor sequence
    receptor_seq = ""
    if os.path.exists(receptor_pdb_path):
        try:
            receptor_seq = extract_receptor_sequence(receptor_pdb_path, chain=chain)
            logger.info(f"Receptor sequence ({len(receptor_seq)} aa): {receptor_seq[:30]}…")
        except Exception as e:
            logger.warning(f"Could not parse receptor sequence: {e}")

    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "receptor_pdb": receptor_pdb_path,
        "interface_data": interface_data,
        "receptor_sequence": receptor_seq,
        "n_interface_residues": len(interface_data.get('site_ii_residues', {}).get('residues', [])),
    }
