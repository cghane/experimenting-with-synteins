"""
Phase 1 — Pocket Analysis
Parse IL-6 receptor structure and characterize the Site II binding pocket
for pharmacophore-based molecule design.
"""

import json
import logging
import os
from typing import Dict

logger = logging.getLogger(__name__)


def analyze_pocket(interface_json):
    # type: (str) -> Dict
    """
    Extract pocket pharmacophore profile from interface_residues.json.
    Returns a dict with hotspot counts and pocket character.
    """
    with open(interface_json, "r") as f:
        data = json.load(f)

    site = data["site_ii_residues"]
    residues = site["residues"]

    hydrophobic_res = [r for r in residues if r["role"].startswith("hydrophobic")]
    charged_res = [r for r in residues if "electrostatic" in r["role"]
                   or "charge" in r["role"] or "salt" in r["role"]]
    hbond_res = [r for r in residues if "H-bond" in r["role"]
                 or "polar" in r["role"] or r["role"] == "backbone H-bond"]
    aromatic_res = [r for r in residues if "aromatic" in r["role"]]

    profile = {
        "n_hydrophobic_hotspots": len(hydrophobic_res) + len(aromatic_res),
        "n_charged_hotspots": len(charged_res),
        "n_hbond_hotspots": len(hbond_res),
        "total_residues": len(residues),
        "hydrophobic_residues": [r["resnum"] for r in hydrophobic_res + aromatic_res],
        "charged_residues": [r["resnum"] for r in charged_res],
        "hbond_residues": [r["resnum"] for r in hbond_res],
        "dominant_character": "electrostatic/polar" if len(charged_res) > len(hydrophobic_res)
                             else "hydrophobic",
    }

    return profile


def run_phase1(interface_json, output_json=None):
    # type: (str, str) -> Dict
    """Execute Phase 1: Pocket analysis."""
    logger.info("=" * 60)
    logger.info("PHASE 1 — Pocket Analysis")
    logger.info("=" * 60)

    profile = analyze_pocket(interface_json)

    logger.info("IL-6 Site II pocket profile:")
    logger.info("  Hydrophobic hotspots: %d", profile["n_hydrophobic_hotspots"])
    logger.info("  Charged hotspots:     %d", profile["n_charged_hotspots"])
    logger.info("  H-bond hotspots:      %d", profile["n_hbond_hotspots"])
    logger.info("  Dominant character:    %s", profile["dominant_character"])
    logger.info("  Total interface res:   %d", profile["total_residues"])

    if output_json:
        os.makedirs(os.path.dirname(output_json), exist_ok=True)
        with open(output_json, "w") as f:
            json.dump(profile, f, indent=2)
        logger.info("Saved pocket profile -> %s", output_json)

    return profile
