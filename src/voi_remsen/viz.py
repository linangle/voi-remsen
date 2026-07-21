from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import PercentFormatter

from .paths import interim_file, output_file

MONTH_ABBR = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

YEARS = list(range(2005, 2027))

TIER_LABELS = {
    0: "no data",
    1: "in-situ physical only (SST, chlorophyll)",
    2: "water hazard only (PN cells, particulate DA)",
    3: "physical + water hazard (week partial)",
    4: "physical + water hazard (full week matchup)",
}
TIER_COLORS = ["#f0f0f0", "#c6dbef", "#fdae6b", "#a1d99b", "#31a354"]
TISSUE_EDGE = "#d62728"


def _year_tier(sub):
    if len(sub) == 0:
        return 0
    haz = set(sub["hazard_truth_status"]); phy = set(sub["insitu_physical_status"])
    has_haz = "has" in haz or "partial" in haz
    has_phy = "has" in phy or "partial" in phy
    if has_haz and has_phy:
        return 4 if ("has" in haz and "has" in phy) else 3
    return 2 if has_haz else (1 if has_phy else 0)


def plot_coverage(weekly_csv=None, cdph_csv=None, out_pdf=None):
    weekly_csv = weekly_csv or interim_file("calhabmap_weekly.csv")
    cdph_csv = cdph_csv or interim_file("cdph_weekly.csv")
    out_pdf = out_pdf or output_file("coverage_heatmap.pdf")

    m = pd.read_csv(weekly_csv, parse_dates=["week_start"])
    m["year"] = m["week_start"].dt.year
    order = list(m["location_name"].dropna().unique())
    code_to_name = m.dropna(subset=["location_name"]).drop_duplicates(
        "location_code").set_index("location_code")["location_name"].to_dict()

    tis = pd.read_csv(cdph_csv, parse_dates=["week_start"])
    tis["year"] = tis["week_start"].dt.year
    tis["location_name"] = tis["location_code"].map(code_to_name)
    tissue_n = tis.groupby(["location_name", "year"]).size().to_dict()

    grid = np.zeros((len(order), len(YEARS)), int)
    counts = np.zeros_like(grid)
    for i, s in enumerate(order):
        for j, y in enumerate(YEARS):
            sub = m[(m["location_name"] == s) & (m["year"] == y)]
            counts[i, j] = len(sub); grid[i, j] = _year_tier(sub)

    cmap = ListedColormap(TIER_COLORS)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.imshow(grid, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(YEARS))); ax.set_xticklabels(YEARS, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=9)

    for i, s in enumerate(order):
        for j, y in enumerate(YEARS):
            if counts[i, j] > 0:
                ax.text(j, i, str(counts[i, j]), ha="center", va="center", fontsize=6,
                        color="white" if grid[i, j] == 4 else "#333")
            nt = tissue_n.get((s, y), 0)
            if nt > 0:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor=TISSUE_EDGE, lw=1.6))
                ax.text(j + 0.38, i + 0.36, str(nt), ha="right", va="bottom",
                        fontsize=5, color=TISSUE_EDGE)

    handles = [Patch(facecolor=TIER_COLORS[k], label=TIER_LABELS[k]) for k in TIER_LABELS]
    handles.append(Patch(facecolor="none", edgecolor=TISSUE_EDGE, linewidth=1.6,
                         label="CDPH tissue DA present (joint-model support; red = # tissue weeks)"))
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    ax.set_title("HAB data coverage by site x year\n"
                 "background = calHABMAP (centre # = weekly rows); red border = CDPH tissue DA",
                 fontsize=11)
    plt.tight_layout(); plt.savefig(out_pdf)
    print(f"saved {out_pdf}")
    return out_pdf

def plot_seasonal_bloom(joint_panel_csv=None, idata_nc=None, out_pdf=None):
    """Monthly bloom probability from the K=2 HMM forward-backward smoother.

    Each observed week gets a posterior P(S_t = bloom | ALL data) from the
    smoother (not a threshold on that week's own cell count), averaged over
    posterior draws. A month's bar is the mean of those probabilities over its
    weeks; the error bar is the 5-95% posterior interval (parameter uncertainty
    propagated through the smoother). The dashed line is the overall posterior
    mean, which equals the model's stationary bloom occupancy (~57%).
    """
    from .worldmodel import build_sequences
    from .worldmodel.selection import smoothed_bloom_prob
    import arviz as az

    joint_panel_csv = joint_panel_csv or interim_file("joint_panel.csv")
    idata_nc = idata_nc or output_file("worldmodel_K2.nc")
    out_pdf = out_pdf or output_file("seasonal_bloom.pdf")

    panel = pd.read_csv(joint_panel_csv, parse_dates=["week_start"])
    seqs = build_sequences(panel)
    idata = az.from_netcdf(idata_nc)
    bloom, weeks = smoothed_bloom_prob(idata, seqs, 2)     # (S, n_obs), (n_obs,)
    months = pd.to_datetime(weeks).month.to_numpy()

    S = bloom.shape[0]
    monthly = np.empty((S, 12))       # per-draw monthly mean prob
    counts = []
    for i, m in enumerate(range(1, 13)):
        mask = months == m
        counts.append(int(mask.sum()))
        monthly[:, i] = bloom[:, mask].mean(1)
    mean_m = monthly.mean(0)
    lo, hi = np.quantile(monthly, [0.05, 0.95], axis=0)
    overall = float(bloom.mean())

    ABOVE, BELOW = "#2a78d6", "#a9ccf3"      # two steps of one blue ramp
    colors = [ABOVE if f >= overall else BELOW for f in mean_m]

    fig, ax = plt.subplots(figsize=(10, 5.0), layout="constrained")
    x = np.arange(12)
    ax.bar(x, mean_m, width=0.66, color=colors, zorder=3)
    ax.errorbar(x, mean_m, yerr=[mean_m - lo, hi - mean_m], fmt="none",
                ecolor="#3a3a3a", elinewidth=1.0, capsize=2.5, zorder=5)

    ax.axhline(overall, color="#5b5b5b", linestyle=(0, (5, 4)), linewidth=1.1, zorder=2)
    ax.text(12.3, overall + 0.02,
            f"overall posterior mean · {overall:.0%}  (= stationary occupancy)",
            va="bottom", ha="right", fontsize=8.5, color="#5b5b5b", zorder=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1})

    for xi, (m, h) in enumerate(zip(mean_m, hi)):
        ax.text(xi, h + 0.02, f"{m:.0%}", ha="center", va="bottom",
                fontsize=8.5, color="#3a3a3a", zorder=6)

    ax.set_xticks(x)
    ax.set_xticklabels(MONTH_ABBR, fontsize=9.5, color="#3a3a3a")
    ax.tick_params(axis="x", length=0, pad=6)
    ax.set_xlim(-0.7, 12.4)

    ax.set_ylim(0, 1.0)
    ax.set_ylabel("P(week in bloom state | data)", fontsize=10, color="#3a3a3a")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.tick_params(axis="y", labelsize=8.5, color="#bbbbbb", labelcolor="#5b5b5b")

    legend = [Patch(facecolor=ABOVE, label="above average"),
              Patch(facecolor=BELOW, label="below average")]
    ax.legend(handles=legend, loc="upper left", fontsize=8.5, frameon=False,
              handlelength=1.1, handleheight=1.1, bbox_to_anchor=(0.005, 1.0))

    fig.suptitle("Monthly bloom probability", fontsize=14.5, fontweight="semibold",
                 x=0.012, ha="left")
    nlo, nhi = min(counts), max(counts)
    ax.set_title("K=2 HMM forward-backward smoother, posterior P(bloom | all data)   ·   "
                 f"error bars 5–95% posterior   ·   {len(months):,} weeks, "
                 f"{nlo}–{nhi} per month",
                 fontsize=9.5, color="#6a6a6a", loc="left", pad=8)

    ax.grid(axis="y", color="#e7e5e0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#bbbbbb")

    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"saved {out_pdf}")
    return out_pdf
