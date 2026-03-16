"""
Visualization module for the small molecule pipeline.
Generates 5 publication-quality figures.
"""

import logging
import os
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw

logger = logging.getLogger(__name__)

# Style constants (match protein pipeline)
ACCENT  = "#2563EB"
WARN    = "#DC2626"
PASS_CLR = "#16A34A"
FAIL_CLR = "#9CA3AF"
BG_CLR   = "#F8FAFC"

plt.rcParams.update({
    "figure.facecolor": BG_CLR,
    "axes.facecolor": "white",
    "axes.edgecolor": "#D1D5DB",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 10,
    "figure.dpi": 150,
})


def _savefig(fig, path):
    # type: (plt.Figure, str) -> None
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("  Saved figure -> %s", path)


def plot_property_distributions(df_all, output_dir):
    # type: (pd.DataFrame, str) -> None
    """Fig 1: Property distributions colored by filter pass/fail."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("Small Molecule Library — Property Distributions", fontsize=14, fontweight="bold")

    properties = [
        ("mw", "Molecular Weight (Da)", None),
        ("logp", "LogP", None),
        ("hbd", "H-Bond Donors", None),
        ("hba", "H-Bond Acceptors", None),
        ("tpsa", "TPSA (\u00c5\u00b2)", None),
        ("rotatable_bonds", "Rotatable Bonds", None),
    ]

    pass_col = "passes_all" if "passes_all" in df_all.columns else None

    for ax, (col, label, _) in zip(axes.ravel(), properties):
        if col not in df_all.columns:
            ax.set_visible(False)
            continue

        if pass_col and pass_col in df_all.columns:
            pass_vals = df_all.loc[df_all[pass_col], col].dropna()
            fail_vals = df_all.loc[~df_all[pass_col], col].dropna()
            ax.hist([pass_vals, fail_vals], bins=25, stacked=True,
                    color=[PASS_CLR, FAIL_CLR], label=["Pass", "Fail"], alpha=0.8)
            ax.legend(fontsize=8)
        else:
            ax.hist(df_all[col].dropna(), bins=25, color=ACCENT, alpha=0.8)

        ax.set_xlabel(label)
        ax.set_ylabel("Count")

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig1_property_distributions.png"))


def plot_filter_funnel(df_all, output_dir):
    # type: (pd.DataFrame, str) -> None
    """Fig 2: Filter attrition funnel."""
    filter_names = ["lipinski", "veber", "complexity", "aggregation", "pains"]
    counts = {"Generated": len(df_all)}

    for name in filter_names:
        col = "pass_%s" % name
        if col in df_all.columns:
            counts[name.title()] = int(df_all[col].sum())

    if "passes_all" in df_all.columns:
        counts["All Filters"] = int(df_all["passes_all"].sum())

    labels = list(counts.keys())
    values = list(counts.values())

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [ACCENT] + [PASS_CLR if v > values[-1] else WARN for v in values[1:]]
    bars = ax.barh(labels[::-1], values[::-1], color=colors[::-1], height=0.6)

    for bar, val in zip(bars, values[::-1]):
        ax.text(bar.get_width() + 5, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=10, fontweight="bold")

    ax.set_xlabel("Number of Molecules")
    ax.set_title("Filter Attrition Funnel", fontsize=13, fontweight="bold")
    ax.set_xlim(0, max(values) * 1.15)

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig2_filter_funnel.png"))


def plot_score_distribution(df_scored, output_dir):
    # type: (pd.DataFrame, str) -> None
    """Fig 3: Score distribution + component breakdown."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: score histogram
    scores = df_scored["final_score"].dropna()
    ax1.hist(scores, bins=30, color=ACCENT, alpha=0.8, edgecolor="white")
    top10_thresh = df_scored["final_score"].nlargest(10).min()
    ax1.axvline(top10_thresh, color=WARN, linestyle="--", linewidth=2,
                label="Top-10 threshold (%.1f)" % top10_thresh)
    ax1.set_xlabel("Final Score")
    ax1.set_ylabel("Count")
    ax1.set_title("Score Distribution", fontweight="bold")
    ax1.legend(fontsize=9)

    # Right: score component breakdown (means)
    components = ["norm_qed", "norm_sa", "norm_binding", "norm_logp", "norm_tpsa"]
    labels = ["QED", "SA (inv)", "Binding\nProxy", "LogP\nPenalty", "TPSA\nBonus"]
    available = [c for c in components if c in df_scored.columns]
    avail_labels = [labels[components.index(c)] for c in available]

    if available:
        means = [df_scored[c].mean() for c in available]
        colors = [PASS_CLR if m > 0.5 else ACCENT for m in means]
        ax2.bar(avail_labels, means, color=colors, alpha=0.85, edgecolor="white")
        ax2.set_ylabel("Mean Normalized Value")
        ax2.set_title("Score Components (Mean)", fontweight="bold")
        ax2.set_ylim(0, 1.1)

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig3_score_distribution.png"))


def plot_optimization_trajectories(history_df, output_dir):
    # type: (pd.DataFrame, str) -> None
    """Fig 4: Optimization trajectories per molecule + population mean."""
    if history_df.empty:
        logger.warning("No optimization history — skipping Fig 4")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    mol_ids = history_df["molecule_id"].unique()
    cmap = plt.cm.get_cmap("tab10")

    # Left: individual trajectories
    for i, mid in enumerate(mol_ids):
        sub = history_df[history_df["molecule_id"] == mid]
        color = cmap(i % 10)
        ax1.plot(sub["round"], sub["best_score"], "o-", color=color,
                 markersize=4, linewidth=1.5, label=mid, alpha=0.8)

    ax1.set_xlabel("Round")
    ax1.set_ylabel("Best Score")
    ax1.set_title("Individual Trajectories", fontweight="bold")
    if len(mol_ids) <= 10:
        ax1.legend(fontsize=7, loc="lower right")

    # Right: population mean +/- SD
    grouped = history_df.groupby("round")["best_score"]
    means = grouped.mean()
    stds = grouped.std().fillna(0)

    ax2.plot(means.index, means.values, "o-", color=ACCENT, linewidth=2, markersize=5)
    ax2.fill_between(means.index, means - stds, means + stds,
                     color=ACCENT, alpha=0.15)
    ax2.set_xlabel("Round")
    ax2.set_ylabel("Mean Best Score")
    ax2.set_title("Population Mean \u00b1 SD", fontweight="bold")

    plt.tight_layout()
    _savefig(fig, os.path.join(output_dir, "fig4_optimization_trajectories.png"))


def plot_lead_dashboard(lead_df, output_dir):
    # type: (pd.DataFrame, str) -> None
    """Fig 5: Lead compound dashboard with 2D structures + radar chart."""
    n_leads = min(len(lead_df), 5)
    if n_leads == 0:
        logger.warning("No leads — skipping Fig 5")
        return

    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(2, n_leads, height_ratios=[1, 1], hspace=0.4, wspace=0.3)

    # Top row: 2D molecular structures
    for i in range(n_leads):
        ax = fig.add_subplot(gs[0, i])
        row = lead_df.iloc[i]
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol:
            img = Draw.MolToImage(mol, size=(250, 250))
            ax.imshow(img)
        ax.set_title("%s\nScore: %.1f" % (row.get("lead_id", "Lead %d" % (i + 1)),
                                            row.get("qed", 0) * 100),
                     fontsize=9, fontweight="bold")
        ax.axis("off")

    # Bottom row: radar chart comparing all leads
    ax_radar = fig.add_subplot(gs[1, :], polar=True)

    categories = ["MW\n(norm)", "LogP\n(norm)", "QED", "SA\n(inv)", "TPSA\n(norm)", "HBD\n(norm)"]
    n_cats = len(categories)
    angles = [n / float(n_cats) * 2 * np.pi for n in range(n_cats)]
    angles += angles[:1]

    for i in range(n_leads):
        row = lead_df.iloc[i]
        values = [
            min(row.get("mw", 0) / 500.0, 1.0),
            min(row.get("logp", 0) / 5.0, 1.0),
            row.get("qed", 0),
            1.0 - row.get("sa_score", 5) / 10.0,
            min(row.get("tpsa", 0) / 140.0, 1.0),
            min(row.get("hbd", 0) / 5.0, 1.0),
        ]
        values += values[:1]
        color = plt.cm.get_cmap("tab10")(i % 10)
        ax_radar.plot(angles, values, "o-", linewidth=1.5, color=color,
                      label=row.get("lead_id", "Lead %d" % (i + 1)), markersize=3)
        ax_radar.fill(angles, values, alpha=0.05, color=color)

    ax_radar.set_xticks(angles[:-1])
    ax_radar.set_xticklabels(categories, fontsize=8)
    ax_radar.set_ylim(0, 1.1)
    ax_radar.set_title("Lead Compound Property Profiles", fontweight="bold", pad=20)
    ax_radar.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.3, 1.1))

    _savefig(fig, os.path.join(output_dir, "fig5_lead_dashboard.png"))


def generate_all_figures(df_all, df_scored, history_df, lead_df, output_dir):
    # type: (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str) -> None
    """Generate all 5 figures."""
    logger.info("Generating figures...")
    os.makedirs(output_dir, exist_ok=True)

    plot_property_distributions(df_all, output_dir)
    plot_filter_funnel(df_all, output_dir)
    plot_score_distribution(df_scored, output_dir)
    plot_optimization_trajectories(history_df, output_dir)
    plot_lead_dashboard(lead_df, output_dir)

    logger.info("All 5 figures generated.")
