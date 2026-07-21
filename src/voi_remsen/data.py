from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C
from .paths import input_file, interim_file

def load_calhabmap_raw(path=None) -> pd.DataFrame:
    path = path or input_file(C.CALHABMAP_FILE)
    df = pd.read_csv(path, low_memory=False)
    df["time_utc"] = pd.to_datetime(df["time_utc"], errors="coerce")
    df = df.dropna(subset=["time_utc"])
    df["week_start"] = df["time_utc"].dt.to_period(C.WEEK_RULE).dt.start_time
    df["pn_total_cells_l"] = df[C.PN_COLUMNS].sum(axis=1, min_count=1)
    return df


def build_calhabmap_weekly(path=None, sites=None, save=True) -> pd.DataFrame:
    df = load_calhabmap_raw(path)
    if sites is not None:
        df = df[df["location_code"].isin(sites)]
    agg = {**C.GROUND_TRUTH_HAZARD, **C.INSITU_PROXY}
    present = {src: dst for src, dst in agg.items() if src in df.columns}
    g = df.groupby(["location_name", "location_code", "week_start"])
    weekly = g[list(present)].mean().rename(columns=present).reset_index()
    weekly["n_samples"] = g.size().values
    weekly = weekly.sort_values(["location_code", "week_start"])

    def status(cols):
        have = weekly[[c for c in cols if c in weekly]].notna().sum(axis=1)
        return have.map(lambda k: "none" if k == 0 else ("has" if k == len(cols) else "partial"))
    weekly["hazard_truth_status"] = status(["pda_water_ng_ml", "pn_total_cells_l"])
    weekly["insitu_physical_status"] = status(["insitu_sst_c", "insitu_chl_mg_m3"])
    if save:
        weekly.to_csv(interim_file("calhabmap_weekly.csv"), index=False)
    return weekly

def load_cdph_raw(path=None) -> pd.DataFrame:
    """CDPH tissue DA (mussel/oyster, co-located sites) with censoring + LOD."""
    path = path or input_file(C.CDPH_FILE)
    df = pd.read_excel(path)
    df.columns = [c.strip() for c in df.columns]
    df["date"] = pd.to_datetime(df["Date -Sampled-"], errors="coerce")
    df = df.dropna(subset=["date"])
    df["week_start"] = df["date"].dt.to_period(C.WEEK_RULE).dt.start_time
    df["da"] = pd.to_numeric(df["ASP -ug/g-"], errors="coerce")
    df["cens"] = (df["Mod-ASP"].astype(str).str.strip() == "<").astype(int)
    df["lod"] = np.where(df["cens"] == 1, df["da"], np.nan)
    typ = df["Sample Type"].astype(str).str.lower()
    grp = np.where(typ.str.contains("mussel"), "mussel",
                   np.where(typ.str.contains("oyster"), "oyster", "other"))
    df = df[np.isin(grp, C.SHELLFISH_GROUPS)]
    inv = {v: k for k, v in C.CO_LOCATED_SITES.items()}
    df = df[df["Sample Site"].isin(C.CO_LOCATED_SITES.values())].copy()
    df["location_code"] = df["Sample Site"].map(inv)
    df["da_cat"] = np.where(df["cens"] == 1, 1, np.where(df["da"] >= C.CLOSURE_PPM, 3, 2))
    return df.dropna(subset=["da"])


def build_cdph_weekly(path=None, save=True) -> pd.DataFrame:
    df = load_cdph_raw(path)
    wk = df.groupby(["location_code", "week_start"]).agg(
        da=("da", "max"), da_cat=("da_cat", "max"),
        lod=("lod", "max"), n=("da", "size")).reset_index()
    wk["cens"] = (wk["da_cat"] == 1).astype(int)
    if save:
        wk.to_csv(interim_file("cdph_weekly.csv"), index=False)
    return wk

def _add_time_index(df, key="location_code", tcol="week_start"):
    df = df.sort_values([key, tcol]).copy()
    df["t_days"] = df.groupby(key)[tcol].transform(lambda x: (x - x.min()).dt.days)
    df["t_years"] = df["t_days"] / 365.25
    return df


def build_cell_panel(save=True) -> pd.DataFrame:
    weekly = build_calhabmap_weekly(sites=list(C.CO_LOCATED_SITES), save=False)
    cell = weekly[weekly[C.BIOMASS_COL].notna()][
        ["location_code", "location_name", "week_start", C.BIOMASS_COL]].copy()
    cell = cell.rename(columns={C.BIOMASS_COL: "cells"})
    cell["cells_log10"] = np.log10(cell["cells"] + 1.0)
    cell = _add_time_index(cell)
    if save:
        cell.to_csv(interim_file("cell_panel.csv"), index=False)
    return cell


def build_toxicity_panel(save=True) -> pd.DataFrame:
    wk = _add_time_index(build_cdph_weekly(save=False))
    codes = {c: i + 1 for i, c in enumerate(sorted(wk["location_code"].unique()))}
    wk["subj"] = wk["location_code"].map(codes)
    out = wk[["subj", "location_code", "week_start", "t_years",
              "da_cat", "da", "cens", "lod", "n"]]
    if save:
        out.to_csv(interim_file("toxicity_panel.csv"), index=False)
    return out


def build_joint_panel(save=True) -> pd.DataFrame:
    cell = build_cell_panel(save=False)[
        ["location_code", "week_start", "cells", "cells_log10"]]
    tox = build_cdph_weekly(save=False)[
        ["location_code", "week_start", "da", "da_cat", "cens", "lod"]]
    panel = pd.merge(cell, tox, on=["location_code", "week_start"], how="outer")
    panel = _add_time_index(panel)
    if save:
        panel.to_csv(interim_file("joint_panel.csv"), index=False)
    return panel
