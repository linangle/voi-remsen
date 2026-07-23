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

EARTH_R_KM = 6371.0088


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (vectorised over lat1/lon1)."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * EARTH_R_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def station_coords(path=None) -> pd.DataFrame:
    """calHABMAP station positions (one row per location_code)."""
    df = load_calhabmap_raw(path)
    st = (df.groupby(["location_code", "location_name"])
            [["latitude_degrees_north", "longitude_degrees_east"]]
            .mean().reset_index()
            .rename(columns={"latitude_degrees_north": "lat",
                             "longitude_degrees_east": "lon"}))
    return st.dropna(subset=["lat", "lon"]).reset_index(drop=True)


def load_cdph_raw(path=None, radius_km=None, cal_path=None) -> pd.DataFrame:
    """CDPH tissue DA (mussel/oyster) with censoring + LOD, matched SPATIALLY.

    Each tissue sample is assigned to its nearest calHABMAP station within
    `radius_km`. Matching on exact `Sample Site` strings (the previous behaviour)
    silently dropped real tissue records whose site label never equals a station
    name, so coverage was understated; coordinates are the reliable key.
    """
    radius_km = C.MATCH_RADIUS_KM if radius_km is None else radius_km
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
    df = df[np.isin(grp, C.SHELLFISH_GROUPS)].copy()

    # nearest station within the radius (coordinates, not site names)
    st = station_coords(cal_path)
    lat = pd.to_numeric(df["Latitude"], errors="coerce").to_numpy(float)
    lon = pd.to_numeric(df["Longitude"], errors="coerce").to_numpy(float)
    D = np.stack([haversine_km(lat, lon, r.lat, r.lon)
                  for r in st.itertuples()], axis=1)          # (n_obs, n_station)
    D = np.where(np.isfinite(D), D, np.inf)
    j = D.argmin(axis=1)
    dmin = D[np.arange(len(D)), j]
    hit = dmin <= radius_km
    df["location_code"] = np.where(hit, st["location_code"].to_numpy()[j], None)
    df["match_km"] = np.where(hit, dmin, np.nan)
    df = df[df["location_code"].notna()]

    df["da_cat"] = np.where(df["cens"] == 1, 1,
                            np.where(df["da"] >= C.CLOSURE_PPM, 3, 2))
    return df.dropna(subset=["da"])


def build_cdph_weekly(path=None, save=True, radius_km=None) -> pd.DataFrame:
    radius_km = C.MATCH_RADIUS_KM if radius_km is None else radius_km
    df = load_cdph_raw(path, radius_km=radius_km)
    wk = df.groupby(["location_code", "week_start"]).agg(
        da=("da", "max"), da_cat=("da_cat", "max"), lod=("lod", "max"),
        n=("da", "size"), match_km=("match_km", "min")).reset_index()
    wk["cens"] = (wk["da_cat"] == 1).astype(int)
    if save:
        wk.to_csv(interim_file(f"cdph_weekly_{radius_km:g}km.csv"), index=False)
        wk.to_csv(interim_file("cdph_weekly.csv"), index=False)
    return wk

def _add_time_index(df, key="location_code", tcol="week_start"):
    df = df.sort_values([key, tcol]).copy()
    df["t_days"] = df.groupby(key)[tcol].transform(lambda x: (x - x.min()).dt.days)
    df["t_years"] = df["t_days"] / 365.25
    return df


def build_cell_panel(save=True, sites=None) -> pd.DataFrame:
    """Water-column panel: PN cells and particulate DA (pDA).

    sites=None keeps ALL calHABMAP stations, since spatial matching can now pair
    tissue with stations beyond the four hand-named co-located ones. pDA is carried
    alongside cells because pDA ~ biomass x toxin-quota: in logs it is an additive
    observation of (log B + log Q), which is what identifies the quota axis.
    """
    weekly = build_calhabmap_weekly(sites=sites, save=False)
    keep = ["location_code", "location_name", "week_start", C.BIOMASS_COL]
    if C.PDA_COL in weekly.columns:
        keep.append(C.PDA_COL)
    cell = weekly[weekly[C.BIOMASS_COL].notna() |
                  weekly.get(C.PDA_COL, pd.Series(index=weekly.index)).notna()][keep].copy()
    cell = cell.rename(columns={C.BIOMASS_COL: "cells", C.PDA_COL: "pda"})
    cell["cells_log10"] = np.log10(cell["cells"] + 1.0)
    if "pda" in cell:
        # log10 pDA on the positive part; exact zeros stay 0 for the hurdle
        cell["pda_log10"] = np.where(cell["pda"] > 0, np.log10(cell["pda"]), 0.0)
        cell.loc[cell["pda"].isna(), "pda_log10"] = np.nan
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


def build_joint_panel(save=True, radius_km=None, sites=None) -> pd.DataFrame:
    radius_km = C.MATCH_RADIUS_KM if radius_km is None else radius_km
    cols = ["location_code", "week_start", "cells", "cells_log10"]
    cp = build_cell_panel(save=False, sites=sites)
    cols += [c for c in ("pda", "pda_log10") if c in cp.columns]
    cell = cp[cols]
    tox = build_cdph_weekly(save=False, radius_km=radius_km)[
        ["location_code", "week_start", "da", "da_cat", "cens", "lod"]]
    panel = pd.merge(cell, tox, on=["location_code", "week_start"], how="outer")
    panel = _add_time_index(panel)
    if save:
        panel.to_csv(interim_file(f"joint_panel_{radius_km:g}km.csv"), index=False)
        panel.to_csv(interim_file("joint_panel.csv"), index=False)
    return panel
