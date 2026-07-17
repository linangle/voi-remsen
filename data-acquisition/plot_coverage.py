from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle

from paths import output_file, ensure_dirs

INPUT_CSV = output_file("matchup_table_sitexweek.csv")
OUTPUT_PDF = output_file("coverage_heatmap.pdf")
ANCHOR = "Stearns Wharf"

SITE_ORDER = [
    "Trinidad Pier",
    "Humboldt",
    "Humboldt South Bay",
    "Bodega Marine Lab",
    "Bodega Marine Lab Buoy",
    "Tomales Bay Mouth",
    "Tomales Bay Mid-Channel Buoy",
    "Inner Tomales Bay",
    "Santa Cruz Wharf",
    "Monterey Wharf",
    "Morro Bay Back Bay",
    "Morro Bay Front Bay",
    "Cal Poly Pier",
    "Stearns Wharf",
    "Santa Monica Pier",
    "Newport Beach Pier",
    "Scripps Pier",
]

YEARS = list(range(2005, 2027))

TIER_LABELS = {
    0: "no data",
    1: "physics only",
    2: "hazard-truth only",
    3: "both (some partial)",
    4: "both (full matchup)",
}
TIER_COLORS = ["#f0f0f0", "#c6dbef", "#fdae6b", "#a1d99b", "#31a354"]


def year_tier(sub):
    if len(sub) == 0:
        return 0
    haz = set(sub["hazard_truth_status"])
    phy = set(sub["insitu_physics_status"])
    has_haz = "has" in haz or "partial" in haz
    has_phy = "has" in phy or "partial" in phy
    if has_haz and has_phy:
        both_full = ("has" in haz) and ("has" in phy)
        return 4 if both_full else 3
    if has_haz:
        return 2
    if has_phy:
        return 1
    return 0


def build_grid(m, order):
    grid = np.zeros((len(order), len(YEARS)), dtype=int)
    counts = np.zeros_like(grid)
    for i, s in enumerate(order):
        for j, y in enumerate(YEARS):
            sub = m[(m["location_name"] == s) & (m["year"] == y)]
            counts[i, j] = len(sub)
            grid[i, j] = year_tier(sub)
    return grid, counts


def plot(grid, counts, order):
    cmap = ListedColormap(TIER_COLORS)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.imshow(grid, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(YEARS)))
    ax.set_xticklabels(
        [str(year) for year in YEARS],
        rotation=45,
        ha="right",
        fontsize=8
    )
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=9)

    for i in range(len(order)):
        for j in range(len(YEARS)):
            n = counts[i, j]
            if n > 0:
                ax.text(
                    j, i, str(n),
                    ha="center", va="center", fontsize=6,
                    color="white" if grid[i, j] == 4 else "#333",
                )

    legend = [Patch(facecolor=TIER_COLORS[k], label=TIER_LABELS[k]) for k in TIER_LABELS]
    ax.legend(handles=legend, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, frameon=False)
    ax.set_title(
        "calHABMAP coverage by site x year  (cell = # weekly matchup rows)",
        fontsize=11,
    )

    plt.tight_layout()
    plt.savefig(OUTPUT_PDF)
    print(f"saved {OUTPUT_PDF}")

def main():
    ensure_dirs()
    m = pd.read_csv(INPUT_CSV, parse_dates=["week_start"])
    m["year"] = m["week_start"].dt.year
    order = [s for s in SITE_ORDER if s in m["location_name"].unique()]
    grid, counts = build_grid(m, order)
    plot(grid, counts, order)


if __name__ == "__main__":
    main()
