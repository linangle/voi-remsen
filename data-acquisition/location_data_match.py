from pathlib import Path

import pandas as pd
import numpy as np

# folder containing python file
BASE_DIR = Path(__file__).resolve().parent

# input / output files
INPUT_CSV = BASE_DIR / "calhabmap.csv"
OUTPUT_CSV = BASE_DIR / "matchup_table_sitexweek.csv"

# column to role mapping
GROUND_TRUTH_HAZARD = {
    'pda_ng_ml':'pda_water_ng_ml',
    'tda_ng_ml':'tda_water_ng_ml',
    'pn_total_cells_l':'pn_total_cells_l',
    'pseudo_nitzschia_seriata_group_cells_l':'pn_seriata_cells_l',
    'alexandrium_spp_cells_l':'alexandrium_cells_l',
    'dinophysis_spp_cells_l':'dinophysis_cells_l',
    'total_phytoplankton_cells_l':'total_phyto_cells_l',
}

# define as in-situ to delineate with satellite-derived measurements
INSITU_PROXY = {
    'temp_degree_c':'insitu_sst_c',
    'avg_chloro_mg_m3':'insitu_chl_mg_m3',
    'nitrate_um':'nitrate_um',
    'silicate_um':'silicate_um',
}
agg = {**GROUND_TRUTH_HAZARD, **INSITU_PROXY}


# parse time stamps
def load(path):
    df = pd.read_csv(path)
    df['time_utc'] = pd.to_datetime(df['time_utc'], errors='coerce') # log unparseable timestamps as NaT
    df = df.dropna(subset=['time_utc']) # drop the rows with NaT
    df["week_start"] = df["time_utc"].dt.to_period("W-SUN").dt.start_time
    return df

# combined pseudo-nitzschia count
def combine_pn(df):
    df["pn_total_cells_l"] = df[['pseudo_nitzschia_delicatissima_group_cells_l', 'pseudo_nitzschia_seriata_group_cells_l']].sum(axis=1, min_count = 1)
    return df

# collapse to one row
def reduce_week(g):
    d={}
    for src,dst in agg.items():
        v = g[src].dropna()
        # in a week, there are repeat measurements of same underlying state --> mean to reduce noise
        d[dst] = v.mean() if len(v) else np.nan
    d['n_samples'] = len(g)
    return pd.Series(d)

# build the matchup table
def classify_pair(df, cols):
    # none if neither column is present
    # partial if only one is present
    # has if both are present

    count_present = df[cols].notna().sum(axis=1)

    return np.select(
        [
            count_present == 0,
            count_present == 1,
            count_present == 2,
        ],
        [
            "none",
            "partial",
            "has",
        ],
        default="none"
    )

def build_matchup(df):
    matchup = (
        df.groupby(["location_name", "location_code", "week_start"])
        .apply(reduce_week, include_groups=False)
        .reset_index()
    )

    matchup = matchup.sort_values(["location_name", "week_start"])

    hazard_cols = ["pda_water_ng_ml", "pn_total_cells_l"]

    insitu_cols = ["insitu_sst_c", "insitu_chl_mg_m3"]

    matchup["hazard_truth_status"] = classify_pair(matchup, hazard_cols)
    matchup["insitu_physics_status"] = classify_pair(matchup, insitu_cols)

    return matchup

# output
def main():
    df = load(INPUT_CSV)
    df = combine_pn(df)
    matchup = build_matchup(df)
    matchup.to_csv(OUTPUT_CSV, index=False)

    print(f"Matchup table shape: {matchup.shape}")
    print(f"Saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
