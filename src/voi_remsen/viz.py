from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle

from .paths import interim_file, output_file

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
