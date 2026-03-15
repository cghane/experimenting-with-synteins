# AI-Guided Mini-Protein Binder Design Pipeline

A closed-loop computational pipeline for designing mini-proteins that block **IL-6** (a key driver of autoimmune disease). This project is my attempt to demonstrate the design-to-chemistry workflow at the core of next-generation Syntein-style therapeutics.

---

## Problem Statement

Existing IL-6 blockers (e.g. tocilizumab) are 150 kDa antibodies costing ~$20k/year per patient, made in living cells. A 55 amino acid mini-protein could do the same job, being small enough to chemically synthesize, cheaper to manufacture, and engineerable in ways biology alone can't achieve (D-amino acids, cyclization, helix stapling).

My pipeline automates the full design process: generate candidates → filter for drug-like properties → score by predicted binding quality → iteratively improve the best ones → annotate for chemical modification.

---

## Pipeline Overview

```
Phase 1  Target setup      Download IL-6 structure (PDB 4CNI), define 18-residue binding site
Phase 2  Generation        1,000 Affibody Z-domain candidates (interface-biased + random)
Phase 3  Filtering         1,000 → ~380 (hydrophobicity, charge, aggregation, protease sites)
Phase 4  Structure scoring ESMFold prediction, extract interface confidence + buried surface area
Phase 5  Ranking           Composite score 0–100 across all candidates
Phase 6  Optimization      Monte Carlo loop: mutate top-10, re-score, keep improvements
Phase 7  Syntein layer      Annotate top-5 for D-amino acid swaps, cyclization, helix stapling
```

**Demo run:** ~75 seconds (mock mode) · **Full run:** ~45 min (ESMFold API)

---

## Quick Start

```bash
pip install -r requirements.txt

python3 pipeline.py                     # full pipeline, mock mode (~75s)
python3 pipeline.py --max-candidates 50 # quick demo
python3 pipeline.py --real              # real ESMFold structure prediction
```

---

## Key Results

| Metric | Value |
|---|---|
| Candidates generated | 1,000 |
| Passed developability filters | ~380 |
| Top candidate score (pre-optimization) | 53.1 / 100 |
| Score improvement after Monte Carlo | +5.3 pts average across top-10 |
| Optimization rounds | 6 × 20 mutation trials per candidate |

---

## Scoring Function

Each candidate is scored on binding quality, stability, and therapeutic viability:

```
Score = 0.35 × interface confidence   (ESMFold predicted alignment at binding face)
      + 0.30 × buried surface area    (proxy for binding affinity)
      + 0.15 × H-bond count
      − 0.10 × aggregation risk       (PASTA scale)
      − 0.10 × protease cleavage risk
      +  0–10 robustness bonus        (charge, helix propensity, low immunogenicity)
```

All terms normalized to [0,1]. Final score scaled to 0–100.

---

## Syntein Annotation Layer (Phase 7)

For the top 5 candidates, the pipeline identifies concrete chemical modification sites:

| Modification | Purpose |
|---|---|
| D-amino acid swaps | Block trypsin/chymotrypsin cleavage at exposed K/R/F/Y sites |
| Head-to-tail cyclization | Lock binding conformation, resist exoprotease degradation |
| Helix stapling (i, i+4 pairs) | Pre-organize helix structure to improve binding affinity |
| PEGylation sites | Extend half-life without disrupting binding interface |

These are rational design annotations — i.e., the same decisions I believe a medicinal chemist would make before sending a sequence to synthesis.

---

## Why the Affibody Scaffold?

- **55 amino acids** - within practical range of solid-phase chemical synthesis
- **Clinically validated** - same scaffold class as izokibep (anti-IL-17A, Phase 2 trials)
- **Stable 3-helix bundle** - well-characterized fold, predictable behavior
- **Designable interface** - 13 positions on helices 1 & 2 face the binding site

```
VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDA
         ★★★★★ ★★    ★★★★★★    ★
                     ↑ 13 variable positions (binding face)
```

---

## Output Figures

| Figure | Description |
|---|---|
| Fig 1 | Library property distributions |
| Fig 2 | Filter attrition funnel (1000 → 380) |
| Fig 3 | Score distribution + component breakdown |
| Fig 4 | Top-20 scorecard heatmap |
| Fig 5 | Optimization trajectories per candidate |
| Fig 6 | Protease risk vs. immunogenicity (top-5) |
| Fig 7 | Sequence logo — variable positions in top-20 binders |

All figures auto-generated to `results/figures/`.

---

## Repository Structure

```
├── pipeline.py              Main runner
├── config.py                All parameters (weights, thresholds, API settings)
├── requirements.txt
├── data/
│   ├── receptor.pdb         IL-6 structure
│   └── interface_residues.json
├── src/
│   ├── phase1_target.py
│   ├── phase2_generate.py
│   ├── phase3_filter.py
│   ├── phase4_predict.py
│   ├── phase5_score.py
│   ├── phase6_optimize.py
│   ├── phase7_synteins.py
│   ├── visualize.py
│   └── utils.py
├── results/                 CSVs + figures (auto-generated)
└── report/
    └── technical_report.md  Full 8-page writeup
```

---

## References

- Lin et al. *Science* 2023 — ESMFold
- Löfblom et al. *FEBS Letters* 2010 — Affibody library design
- Boulanger et al. *Science* 2003 — IL-6/gp130 complex structure
- Kyte & Doolittle *J Mol Biol* 1982 — Hydrophobicity scale
- Trovato et al. *PLoS Comput Biol* 2006 — PASTA aggregation scale
