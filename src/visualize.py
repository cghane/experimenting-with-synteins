"""
Visualisation module — generates publication-quality figures for all phases.

Figures produced:
  Fig 1 — Library distributions (Phase 3): hydrophobicity, charge, aggregation
  Fig 2 — Filter funnel: candidates passing each filter
  Fig 3 — Score distributions (Phase 5): composite score histogram
  Fig 4 — Top-20 candidate scorecard: heatmap of normalised metrics
  Fig 5 — Optimisation trajectories (Phase 6): score vs round
  Fig 6 — Synteins summary (Phase 7): protease vs immunogenicity scatter
  Fig 7 — Sequence logo of top-10 binders at variable positions
"""

import logging
import os
from typing import List, Optional

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

from config import VARIABLE_POSITIONS, AFFIBODY_SCAFFOLD, RESULTS_DIR

logger = logging.getLogger(__name__)

# ── Style ─────────────────────────────────────────────────────────────────────
PALETTE  = sns.color_palette("coolwarm_r", 10)
ACCENT   = "#2563EB"   # Abiologics blue
WARN     = "#DC2626"   # Warning red
PASS_CLR = "#16A34A"   # Pass green
FAIL_CLR = "#DC2626"   # Fail red

plt.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
    'axes.spines.top':  False,
    'axes.spines.right': False,
    'font.family':      'sans-serif',
    'font.size':        10,
    'axes.titlesize':   12,
    'axes.labelsize':   10,
})


def _savefig(fig, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"  Saved → {path}")


# ── Figure 1: Library property distributions ──────────────────────────────────

def plot_library_distributions(df_filtered: pd.DataFrame, output_dir: str):
    """
    Phase 3: Violin/box plots of key physicochemical properties,
    coloured by pass/fail status.
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("Phase 3 — Candidate Library Distributions",
                 fontsize=14, fontweight='bold', y=1.01)

    df = df_filtered.copy()
    pass_mask = df['passes_filter']

    props = [
        ('max_hydro_window',  'Max Hydrophobicity Window\n(KD scale, 5-aa)',
         'max_hydrophobic_window', 2.5),
        ('net_charge_calc',   'Net Charge at pH 7.4',
         'net_charge', None),
        ('n_cysteines',       'Cysteine Count',
         'cysteines', 1),
        ('min_entropy',       'Min Local Sequence Entropy\n(12-aa window)',
         'entropy', 2.0),
        ('max_aggregation',   'Max Aggregation Window\n(PASTA, 5-aa)',
         'aggregation', 2.0),
        ('mean_helix_prop',   'Mean Helix Propensity\n(Pace & Scholtz)',
         'helix_prop', 0.65),
    ]

    for ax, (col, label, _, threshold) in zip(axes.flat, props):
        if col not in df.columns:
            ax.set_visible(False)
            continue

        pass_data = df.loc[pass_mask, col].dropna()
        fail_data = df.loc[~pass_mask, col].dropna()

        ax.hist(fail_data, bins=30, alpha=0.6, color=FAIL_CLR, label='Fail', density=True)
        ax.hist(pass_data, bins=30, alpha=0.6, color=PASS_CLR, label='Pass', density=True)

        if threshold is not None:
            ax.axvline(threshold, color='k', linestyle='--', linewidth=1.5,
                       label=f'Threshold = {threshold}')

        ax.set_xlabel(label)
        ax.set_ylabel('Density')
        ax.legend(fontsize=8, framealpha=0.5)
        ax.set_title(label.split('\n')[0])

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig1_library_distributions.png"))


# ── Figure 2: Filter funnel ───────────────────────────────────────────────────

def plot_filter_funnel(df_filtered: pd.DataFrame, output_dir: str):
    """Waterfall / funnel chart showing candidate attrition at each filter."""
    filter_cols = [c for c in df_filtered.columns if c.startswith('pass_')]
    filter_names = [c.replace('pass_', '').replace('_', ' ').title() for c in filter_cols]
    filter_names = ['All Candidates'] + filter_names + ['All Filters']

    # Cumulative pass (each applied independently then combined)
    n_total = len(df_filtered)
    counts = [n_total]
    for col in filter_cols:
        counts.append(int(df_filtered[col].sum()))
    counts.append(int(df_filtered['passes_filter'].sum()))

    fig, ax = plt.subplots(figsize=(12, 5))
    colors = [ACCENT] + [PASS_CLR] * len(filter_cols) + [WARN]
    bars = ax.barh(filter_names, counts, color=colors, alpha=0.85, edgecolor='white')

    for bar, count in zip(bars, counts):
        ax.text(count + 5, bar.get_y() + bar.get_height()/2,
                f'{count} ({100*count/n_total:.0f}%)',
                va='center', fontsize=9)

    ax.set_xlabel('Candidates Remaining')
    ax.set_title('Phase 3 — Filter Attrition Funnel', fontweight='bold')
    ax.set_xlim(0, n_total * 1.18)
    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig2_filter_funnel.png"))


# ── Figure 3: Score distribution ─────────────────────────────────────────────

def plot_score_distribution(df_scored: pd.DataFrame, output_dir: str):
    """Phase 5: Final score distribution with top-20 highlighted."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Phase 5 — Composite Score Distribution", fontsize=13, fontweight='bold')

    scores = df_scored['final_score'].dropna()
    top20_thresh = df_scored['final_score'].nlargest(20).min()

    # Left: histogram
    ax = axes[0]
    ax.hist(scores, bins=40, color=ACCENT, alpha=0.75, edgecolor='white')
    ax.axvline(top20_thresh, color=WARN, linestyle='--', linewidth=2,
               label=f'Top-20 threshold ({top20_thresh:.1f})')
    ax.set_xlabel('Final Score (0–100)')
    ax.set_ylabel('Count')
    ax.set_title('Score Distribution')
    ax.legend()

    # Right: breakdown by score component
    ax2 = axes[1]
    score_cols = ['norm_interface_plddt', 'norm_bsa', 'norm_hbond',
                  'norm_aggregation', 'norm_protease']
    score_labels = ['Interface pLDDT', 'BSA', 'H-bonds', 'Aggregation\n(penalty)', 'Protease\n(penalty)']
    avail = [c for c in score_cols if c in df_scored.columns]

    if avail:
        means = [df_scored[c].mean() for c in avail]
        lbls  = [score_labels[score_cols.index(c)] for c in avail]
        bars  = ax2.bar(lbls, means, color=PALETTE[:len(avail)], alpha=0.85, edgecolor='white')
        for bar, m in zip(bars, means):
            ax2.text(bar.get_x() + bar.get_width()/2, m + 0.01,
                     f'{m:.2f}', ha='center', va='bottom', fontsize=9)
        ax2.set_ylim(0, 1.1)
        ax2.set_ylabel('Normalised Mean Score Component')
        ax2.set_title('Score Component Breakdown\n(Population Mean)')

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig3_score_distribution.png"))


# ── Figure 4: Top-20 scorecard heatmap ───────────────────────────────────────

def plot_top20_scorecard(df_scored: pd.DataFrame, output_dir: str):
    """Heatmap of all normalised metrics for the top-20 candidates."""
    top20 = df_scored.head(20).copy()
    metric_cols = ['interface_plddt', 'bsa_proxy_A2', 'hbond_estimate',
                   'n_contacts', 'max_aggregation', 'protease_risk',
                   'plddt_mean', 'final_score']
    metric_labels = ['Interface\npLDDT', 'BSA (Å²)', 'H-bonds',
                     'Contacts', 'Aggregation\n(↓ better)', 'Protease\n(↓ better)',
                     'Mean pLDDT', 'Final\nScore']

    avail_cols   = [c for c in metric_cols if c in top20.columns]
    avail_labels = [metric_labels[metric_cols.index(c)] for c in avail_cols]

    # Normalise each column to [0, 1] for display
    hmap = top20[avail_cols].copy()
    for c in avail_cols:
        lo, hi = hmap[c].min(), hmap[c].max()
        if hi > lo:
            hmap[c] = (hmap[c] - lo) / (hi - lo)

    fig, ax = plt.subplots(figsize=(len(avail_cols) * 1.4 + 2, 8))
    sns.heatmap(
        hmap.values,
        ax=ax,
        cmap='RdYlGn',
        xticklabels=avail_labels,
        yticklabels=top20['id'].tolist(),
        annot=True,
        fmt='.2f',
        linewidths=0.4,
        linecolor='white',
        cbar_kws={'label': 'Normalised score (0–1)'},
    )
    ax.set_title("Top-20 Candidates — Metric Scorecard", fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=7)
    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig4_top20_scorecard.png"))


# ── Figure 5: Optimisation trajectories ──────────────────────────────────────

def plot_optimisation_trajectories(history_df: pd.DataFrame, output_dir: str):
    """Phase 6: Best score per round for all optimised candidates."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Phase 6 — Iterative Optimisation (Monte Carlo / SA)",
                 fontsize=13, fontweight='bold')

    # Left: individual trajectories
    ax = axes[0]
    cmap = plt.cm.get_cmap('tab10')
    candidates = history_df['candidate_id'].unique()
    for i, cand in enumerate(candidates):
        sub = history_df[history_df['candidate_id'] == cand].sort_values('round')
        color = cmap(i % 10)
        ax.plot(sub['round'], sub['best_score'], '-o', markersize=4,
                color=color, linewidth=1.8, label=cand, alpha=0.85)

    ax.set_xlabel('Optimisation Round')
    ax.set_ylabel('Best Score (0–100)')
    ax.set_title('Score Trajectory per Candidate')
    if len(candidates) <= 10:
        ax.legend(fontsize=7, framealpha=0.5)
    ax.grid(axis='y', alpha=0.3)

    # Right: mean ± std across all candidates
    ax2 = axes[1]
    grouped = history_df.groupby('round')['best_score'].agg(['mean', 'std']).reset_index()
    ax2.plot(grouped['round'], grouped['mean'], '-o', color=ACCENT, linewidth=2, markersize=5)
    ax2.fill_between(grouped['round'],
                     grouped['mean'] - grouped['std'],
                     grouped['mean'] + grouped['std'],
                     color=ACCENT, alpha=0.2, label='±1 SD')

    initial = grouped['mean'].iloc[0]
    final   = grouped['mean'].iloc[-1]
    ax2.annotate(f'+{final - initial:.1f}',
                 xy=(grouped['round'].iloc[-1], final),
                 xytext=(-30, 10), textcoords='offset points',
                 fontsize=11, color=PASS_CLR, fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color=PASS_CLR))

    ax2.set_xlabel('Optimisation Round')
    ax2.set_ylabel('Mean Best Score (0–100)')
    ax2.set_title('Population Mean Improvement')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig5_optimisation_trajectories.png"))


# ── Figure 6: Synteins scatter ────────────────────────────────────────────────

def plot_synteins_analysis(synteins_df: pd.DataFrame, output_dir: str):
    """Phase 7: Protease risk vs immunogenicity scatter with score as bubble size."""
    if synteins_df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Phase 7 — Synteins-Style Safety & Stability Analysis",
                 fontsize=13, fontweight='bold')

    # Left: protease risk vs immunogenicity
    ax = axes[0]
    sc = ax.scatter(
        synteins_df['total_protease_sites'],
        synteins_df['immunogen_score'],
        c=synteins_df['final_score'],
        cmap='RdYlGn',
        s=200, alpha=0.85, edgecolors='gray', linewidths=0.8,
        zorder=3,
    )
    for _, row in synteins_df.iterrows():
        ax.annotate(row['id'], (row['total_protease_sites'], row['immunogen_score']),
                    textcoords='offset points', xytext=(5, 3), fontsize=7)

    cb = plt.colorbar(sc, ax=ax)
    cb.set_label('Final Score')
    ax.set_xlabel('Total Protease Cleavage Sites')
    ax.set_ylabel('Immunogenicity Score (0–1, ↓ better)')
    ax.set_title('Protease Risk vs Immunogenicity')
    ax.axvline(3, color=WARN, linestyle='--', alpha=0.5, label='Max sites = 3')
    ax.axhline(0.4, color=WARN, linestyle=':', alpha=0.5, label='Max immuno = 0.4')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)

    # Right: modification opportunity count
    ax2 = axes[1]
    mod_data = {
        'Staple Pairs\n(i,i+4)':  synteins_df['n_staple_pairs'],
        'PEGylation\nSites':       synteins_df['n_peg_sites'],
        'D-AA Swap\nSites':        synteins_df['n_daa_sites'],
    }
    x = np.arange(len(synteins_df))
    width = 0.25
    colors_mods = [ACCENT, PASS_CLR, WARN]
    for i, (label, vals) in enumerate(mod_data.items()):
        ax2.bar(x + i * width, vals.values, width, label=label, color=colors_mods[i], alpha=0.8)

    ax2.set_xticks(x + width)
    ax2.set_xticklabels(synteins_df['id'].tolist(), rotation=30, ha='right', fontsize=8)
    ax2.set_ylabel('Modification Opportunities')
    ax2.set_title('Available Stabilisation Modifications')
    ax2.legend(fontsize=8)
    ax2.grid(axis='y', alpha=0.2)

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig6_synteins_analysis.png"))


# ── Figure 7: Sequence logo ───────────────────────────────────────────────────

def plot_sequence_logo(df_scored: pd.DataFrame, output_dir: str, top_n: int = 20):
    """Sequence logo at variable positions for top-N binders."""
    try:
        import logomaker
    except ImportError:
        logger.warning("logomaker not installed — skipping sequence logo")
        return

    top = df_scored.head(top_n)
    if top.empty:
        return

    # Build position frequency matrix for variable positions only
    counts = {pos: {} for pos in VARIABLE_POSITIONS}
    for seq in top['sequence']:
        for pos in VARIABLE_POSITIONS:
            aa = seq[pos]
            counts[pos][aa] = counts[pos].get(aa, 0) + 1

    all_aas = list("ACDEFGHIKLMNPQRSTVWY")
    matrix = pd.DataFrame(
        {i: {aa: counts[pos].get(aa, 0) / top_n for aa in all_aas}
         for i, pos in enumerate(VARIABLE_POSITIONS)}
    ).T

    fig, ax = plt.subplots(figsize=(max(8, len(VARIABLE_POSITIONS) * 0.8), 3))
    logo = logomaker.Logo(matrix, ax=ax, color_scheme='chemistry')
    logo.style_xticks(rotation=0, fmt='%d')
    ax.set_xticks(range(len(VARIABLE_POSITIONS)))
    ax.set_xticklabels([f"{AFFIBODY_SCAFFOLD[p]}{p+1}" for p in VARIABLE_POSITIONS], fontsize=8)
    ax.set_xlabel('Variable Position (scaffold residue + position)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Sequence Logo — Variable Positions (Top {top_n} Binders)', fontweight='bold')
    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig7_sequence_logo.png"))


# ── Master plot function ──────────────────────────────────────────────────────

def generate_all_figures(
    df_filtered:   Optional[pd.DataFrame] = None,
    df_scored:     Optional[pd.DataFrame] = None,
    history_df:    Optional[pd.DataFrame] = None,
    synteins_df:   Optional[pd.DataFrame] = None,
    output_dir:    str = os.path.join(RESULTS_DIR, "figures"),
):
    """Generate all figures for the pipeline report."""
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"\nGenerating figures → {output_dir}")

    if df_filtered is not None:
        logger.info("  Fig 1: Library distributions …")
        plot_library_distributions(df_filtered, output_dir)
        logger.info("  Fig 2: Filter funnel …")
        plot_filter_funnel(df_filtered, output_dir)

    if df_scored is not None:
        logger.info("  Fig 3: Score distribution …")
        plot_score_distribution(df_scored, output_dir)
        logger.info("  Fig 4: Top-20 scorecard …")
        plot_top20_scorecard(df_scored, output_dir)
        logger.info("  Fig 7: Sequence logo …")
        plot_sequence_logo(df_scored, output_dir)

    if history_df is not None:
        logger.info("  Fig 5: Optimisation trajectories …")
        plot_optimisation_trajectories(history_df, output_dir)

    if synteins_df is not None:
        logger.info("  Fig 6: Synteins analysis …")
        plot_synteins_analysis(synteins_df, output_dir)

    logger.info(f"All figures saved to {output_dir}")
