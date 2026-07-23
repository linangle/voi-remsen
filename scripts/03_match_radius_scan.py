"""Compare name-based vs coordinate-based tissue/cell matching at 3 km and 5 km.

    python scripts/03_match_radius_scan.py

Writes, for each radius: a matched weekly tissue table, a joint panel, a coverage
heatmap, and a per-station matchup table. Prints the totals side by side.
"""

import numpy as np
import pandas as pd

from voi_remsen import config as C, paths
from voi_remsen.data import (build_calhabmap_weekly, build_cdph_weekly,
                             build_joint_panel, load_cdph_raw, station_coords)
from voi_remsen.viz import plot_coverage

RADII = [3.0, 5.0]


def summarise(panel):
    has_c = panel["cells"].notna()
    has_d = panel["da"].notna()
    return dict(
        site_weeks=len(panel), cell_weeks=int(has_c.sum()),
        tissue_weeks=int(has_d.sum()), matchups=int((has_c & has_d).sum()),
        closure_weeks=int((panel["da_cat"] == 3).sum()),
        stations=int(panel.loc[has_d, "location_code"].nunique()),
    )


def main():
    build_calhabmap_weekly(sites=None, save=True)      # all 17 stations
    st = station_coords()
    print(f"calHABMAP stations with coordinates: {len(st)}\n")

    rows, per_station = [], {}
    for r in RADII:
        raw = load_cdph_raw(radius_km=r)
        panel = build_joint_panel(save=True, radius_km=r)
        build_cdph_weekly(save=True, radius_km=r)
        s = summarise(panel)
        s["radius_km"] = r
        s["tissue_samples_matched"] = len(raw)
        s["median_match_km"] = round(float(raw["match_km"].median()), 2)
        rows.append(s)

        ps = (panel.assign(both=panel["cells"].notna() & panel["da"].notna())
                   .groupby("location_code")
                   .agg(cell_weeks=("cells", "count"),
                        tissue_weeks=("da", "count"),
                        matchups=("both", "sum"),
                        closures=("da_cat", lambda x: int((x == 3).sum())))
                   .sort_values("matchups", ascending=False))
        ps = ps[ps.tissue_weeks > 0]
        per_station[r] = ps
        ps.to_csv(paths.output_file(f"matchups_by_station_{r:g}km.csv"))

        plot_coverage(cdph_csv=paths.interim_file(f"cdph_weekly_{r:g}km.csv"),
                      out_pdf=paths.output_file(f"coverage_heatmap_{r:g}km.pdf"))

    summ = pd.DataFrame(rows)[["radius_km", "stations", "tissue_samples_matched",
                               "median_match_km", "cell_weeks", "tissue_weeks",
                               "matchups", "closure_weeks", "site_weeks"]]
    summ.to_csv(paths.output_file("match_radius_summary.csv"), index=False)
    print("=== coordinate matching: totals by radius ===")
    print(summ.to_string(index=False))
    for r in RADII:
        print(f"\n=== stations with tissue, {r:g} km ===")
        print(per_station[r].to_string())


if __name__ == "__main__":
    main()
