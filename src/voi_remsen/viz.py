from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import PercentFormatter

from .config import MODEL_SITES
from .paths import interim_file, output_file

MONTH_ABBR = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

YEARS = list(range(2005, 2027))

CELL_COL, PDA_COL = "pn_total_cells_l", "pda_water_ng_ml"

TIER_LABELS = {
    0: "no data",
    1: "PN cell counts only",
    2: "particulate DA (pDA) only",
    3: "both cells + pDA",
}
TIER_COLORS = ["#f0f0f0", "#9ec5f4", "#f6c667", "#7cc47f"]
TISSUE_EDGE = "#d62728"


def _year_tier(sub):
    """Which water-column hazard markers this site-year has (SST/chl ignored)."""
    if len(sub) == 0:
        return 0
    has_cells = sub[CELL_COL].notna().any() if CELL_COL in sub else False
    has_pda = sub[PDA_COL].notna().any() if PDA_COL in sub else False
    return 3 if (has_cells and has_pda) else (1 if has_cells else (2 if has_pda else 0))


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
    n_cell = np.zeros_like(grid); n_pda = np.zeros_like(grid)
    for i, s in enumerate(order):
        for j, y in enumerate(YEARS):
            sub = m[(m["location_name"] == s) & (m["year"] == y)]
            grid[i, j] = _year_tier(sub)
            if len(sub):
                n_cell[i, j] = int(sub[CELL_COL].notna().sum()) if CELL_COL in sub else 0
                n_pda[i, j] = int(sub[PDA_COL].notna().sum()) if PDA_COL in sub else 0

    cmap = ListedColormap(TIER_COLORS)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.imshow(grid, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(YEARS))); ax.set_xticklabels(YEARS, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=9)

    for i, s in enumerate(order):
        for j, y in enumerate(YEARS):
            nc, np_ = n_cell[i, j], n_pda[i, j]
            if nc or np_:
                # cells top-left, pDA bottom-left; single number when only one marker
                if nc and np_:
                    ax.text(j - 0.42, i - 0.12, str(nc), ha="left", va="center",
                            fontsize=5.5, color="#1a3e5c")
                    ax.text(j - 0.42, i + 0.28, str(np_), ha="left", va="center",
                            fontsize=5.5, color="#7a5200")
                else:
                    ax.text(j - 0.42, i + 0.06, str(nc or np_), ha="left",
                            va="center", fontsize=6, color="#333")
            nt = tissue_n.get((s, y), 0)
            if nt > 0:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor=TISSUE_EDGE, lw=1.6))
                ax.text(j + 0.44, i + 0.34, str(nt), ha="right", va="bottom",
                        fontsize=5, color=TISSUE_EDGE)

    handles = [Patch(facecolor=TIER_COLORS[k], label=TIER_LABELS[k]) for k in TIER_LABELS]
    handles.append(Patch(facecolor="none", edgecolor=TISSUE_EDGE, linewidth=1.6,
                         label="CDPH tissue DA matched (red no. = tissue weeks)"))
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, frameon=False)
    ax.set_title("HAB water-column data coverage by site x year\n"
                 "fill = which calHABMAP markers are present; numbers = weeks with "
                 "cells (blue) / pDA (amber); red border = matched CDPH tissue DA",
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
    mean of these smoothed probabilities -- an observation-weighted data average,
    NOT the CTMC equilibrium u/(u+d) (close here, but different quantities: the
    average depends on the sampling schedule and finite sequences).
    """
    from .worldmodel import build_sequences
    from .worldmodel.selection import smoothed_bloom_prob
    import arviz as az

    joint_panel_csv = joint_panel_csv or interim_file("joint_panel.csv")
    K = int(np.load(output_file("worldmodel.npz"))["K"])   # production K
    idata_nc = idata_nc or output_file(f"worldmodel_K{K}.nc")
    out_pdf = out_pdf or output_file("seasonal_bloom.pdf")

    panel = pd.read_csv(joint_panel_csv, parse_dates=["week_start"])
    seqs = build_sequences(panel, sites=MODEL_SITES)
    idata = az.from_netcdf(idata_nc)
    bloom, weeks = smoothed_bloom_prob(idata, seqs, K)     # top state (K-1)
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
    ax.text(12.3, overall + float(hi.max()) * 0.06,
            f"overall posterior-mean P(toxic) · {overall:.1%}",
            va="bottom", ha="right", fontsize=8.5, color="#5b5b5b", zorder=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1})

    for xi, (m, h) in enumerate(zip(mean_m, hi)):
        ax.text(xi, h + float(hi.max()) * 0.06, f"{m:.1%}", ha="center", va="bottom",
                fontsize=8.5, color="#3a3a3a", zorder=6)

    ax.set_xticks(x)
    ax.set_xticklabels(MONTH_ABBR, fontsize=9.5, color="#3a3a3a")
    ax.tick_params(axis="x", length=0, pad=6)
    ax.set_xlim(-0.7, 12.4)

    ax.set_ylim(0, max(float(hi.max()) * 1.45, 0.02))
    ax.set_ylabel("P(week in toxic state | data)", fontsize=10, color="#3a3a3a")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.tick_params(axis="y", labelsize=8.5, color="#bbbbbb", labelcolor="#5b5b5b")

    legend = [Patch(facecolor=ABOVE, label="above average"),
              Patch(facecolor=BELOW, label="below average")]
    ax.legend(handles=legend, loc="upper left", fontsize=8.5, frameon=False,
              handlelength=1.1, handleheight=1.1, bbox_to_anchor=(0.005, 1.0))

    fig.suptitle("Monthly toxic-state probability", fontsize=14.5, fontweight="semibold",
                 x=0.012, ha="left")
    nlo, nhi = min(counts), max(counts)
    ax.set_title("SCW forward-backward smoother, posterior P(top/toxic state | all data)   ·   "
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


def plot_seasonal_cycle(seasonal_npz=None, stat_nc=None, joint_panel_csv=None,
                        out_pdf=None):
    """Phase D: the fitted seasonal bloom cycle vs the empirical monthly curve.

    Line + band: the seasonal model's instantaneous bloom occupancy up(woy) /
    (up(woy)+down(woy)) across week-of-year, with a 90% posterior band from the
    harmonic-coefficient draws. Points: the empirical monthly P(bloom) from the
    *stationary* K=2 forward-backward smoother (a separate fit). Their agreement
    is the cross-check that the parametric seasonal Q captures the real cycle.
    """
    from .worldmodel import build_sequences
    from .worldmodel.seasonal import season_design, WOY, PERIOD
    from .worldmodel.selection import smoothed_bloom_prob
    import arviz as az

    seasonal_npz = seasonal_npz or output_file("seasonal_worldmodel.npz")
    stat_nc = stat_nc or output_file(
        f"worldmodel_K{int(np.load(seasonal_npz)['K'])}.nc")
    joint_panel_csv = joint_panel_csv or interim_file("joint_panel.csv")
    out_pdf = out_pdf or output_file("seasonal_cycle.pdf")

    from .worldmodel.linalg import birthdeath_Q, expm_np_batched
    from .worldmodel.seasonal import periodic_occupancy

    d = np.load(seasonal_npz)
    H = int(d["H"])
    Xd = season_design(np.arange(WOY), H)                 # (53, ncoef)
    Kk = int(d["K"])
    cu, cd = d["coef_up"], d["coef_down"]                 # (S, K-1, ncoef)
    up = np.exp(np.einsum("wc,skc->swk", Xd, cu))         # (S,53,K-1)
    down = np.exp(np.einsum("wc,skc->swk", Xd, cd))
    S = up.shape[0]
    # frozen-rate equilibrium of a birth-death chain (diagnostic only)
    ratio = np.cumprod(up / down, axis=2)
    wt = np.concatenate([np.ones((S, WOY, 1)), ratio], axis=2)
    frozen = (wt / wt.sum(2, keepdims=True))[:, :, Kk - 1]
    # TRUE propagated periodic occupancy: build P_woy per draw, then propagate
    Q = np.zeros((S, WOY, Kk, Kk))
    for i in range(Kk - 1):
        Q[:, :, i, i + 1] = up[:, :, i]
        Q[:, :, i + 1, i] = down[:, :, i]
    for k in range(Kk):
        Q[:, :, k, k] = -Q[:, :, k, :].sum(-1)
    P_woy = expm_np_batched(Q)                            # (S,53,K,K)
    occ = periodic_occupancy(P_woy)[:, :, Kk - 1]         # (S,53) top-state prevalence
    occ_m = occ.mean(0)
    lo, hi = np.quantile(occ, [0.05, 0.95], axis=0)
    frozen_m = frozen.mean(0)

    # stationary-model monthly smoother (SAME data -- consistency check only)
    panel = pd.read_csv(joint_panel_csv, parse_dates=["week_start"])
    seqs = build_sequences(panel, sites=MODEL_SITES)
    bloom, weeks = smoothed_bloom_prob(az.from_netcdf(stat_nc), seqs, Kk)
    wk = pd.to_datetime(weeks)
    woy_obs = np.minimum((wk.dayofyear.to_numpy() - 1) // 7, WOY - 1)
    months = wk.month.to_numpy()
    emp_x = np.array([woy_obs[months == m].mean() for m in range(1, 13)])
    emp_y = np.array([bloom[:, months == m].mean() for m in range(1, 13)])

    woy = np.arange(WOY)
    fig, ax = plt.subplots(figsize=(10, 5.0), layout="constrained")
    ax.fill_between(woy, lo, hi, color="#2a78d6", alpha=0.15, zorder=2,
                    label="propagated occupancy, 90% posterior")
    ax.plot(woy, occ_m, color="#2a78d6", lw=2.2, zorder=4,
            label="seasonal model, propagated toxic-state occupancy")
    ax.plot(woy, frozen_m, color="#8f8d84", lw=1.3, ls=(0, (5, 3)), zorder=3,
            label="frozen-rate ratio u/(u+d) (leads prevalence)")
    ax.plot(emp_x, emp_y, "o", ms=7, color="#eb6834", zorder=5,
            markeredgecolor="white", markeredgewidth=1.2,
            label="stationary-model monthly smoother (same data)")

    peak = int(occ_m.argmax())
    ax.axvline(peak, color="#8f8d84", ls=(0, (2, 3)), lw=1, zorder=1)
    ax.annotate(f"occupancy peak ~week {peak}\n({occ_m[peak]:.1%})",
                xy=(peak, occ_m[peak]),
                xytext=(peak + 5, occ_m[peak] * 0.42), fontsize=8.5,
                color="#5b5b5b",
                arrowprops={"arrowstyle": "->", "color": "#8f8d84", "lw": 1})

    month_woy = [int(((pd.Timestamp(2020, mth, 1).dayofyear) - 1) // 7)
                 for mth in range(1, 13)]
    ax.set_xticks(month_woy)
    ax.set_xticklabels(MONTH_ABBR, fontsize=9.5, color="#3a3a3a")
    ax.tick_params(axis="x", length=0, pad=6)
    ax.set_xlim(-0.5, WOY - 0.5)
    ax.set_ylim(0, max(float(hi.max()), float(frozen_m.max())) * 1.35)
    ax.set_ylabel("P(in toxic state)", fontsize=10, color="#3a3a3a")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.tick_params(axis="y", labelsize=8.5, color="#bbbbbb", labelcolor="#5b5b5b")

    ax.legend(loc="upper right", fontsize=8, frameon=False)
    fig.suptitle("Phase D: fitted seasonal toxic-state cycle", fontsize=14.5,
                 fontweight="semibold", x=0.012, ha="left")
    ax.set_title("SCW, non-stationary Q(day-of-year), one annual harmonic   ·   "
                 "propagated occupancy vs the same-data stationary smoother "
                 "(internal consistency, not independent validation)",
                 fontsize=9.5, color="#6a6a6a", loc="left", pad=8)
    ax.grid(axis="y", color="#e7e5e0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#bbbbbb")

    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"saved {out_pdf}")
    return out_pdf
