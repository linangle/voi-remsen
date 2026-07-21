"""Configuration for the voi-remsen world model.

Central place for the modeling constants documented in docs/methodology.md.
Only world-building lives here (biomass + toxicity dynamics); the decision model
(POMDP / accuracy sweep) is downstream and separate.
"""

# --- raw input filenames (in data/input/) --------------------------------
CALHABMAP_FILE = "calhabmap.csv"
CDPH_FILE = "CDPH_tissueDA.xlsx"

# --- co-located anchor sites: calHABMAP station code -> exact CDPH site name
# (mussel/oyster, biomass and toxicity within a few hundred metres)
# Morro Bay dropped: its tissue (1991-2005) and cells (2023+) don't overlap in
# time, so it contributes no paired biomass+toxicity signal (docs/methodology.md).
CO_LOCATED_SITES = {
    "SCW": "Santa Cruz Wharf",                 # anchor: most paired obs + events
    "CPP": "San Luis Obispo, Cal Poly Pier",
    "MW":  "Monterey Bay, Commercial Wharf",
    "SIO": "La Jolla, Scripps Pier",
}

# --- biomass marker (calHABMAP Pseudo-nitzschia cell counts, cells/L) -----
PN_COLUMNS = [
    "pseudo_nitzschia_delicatissima_group_cells_l",
    "pseudo_nitzschia_seriata_group_cells_l",
]

# --- toxicity marker (CDPH domoic acid) ----------------------------------
# Restrict to fast-depurating shellfish (mussels + oysters); razor clams are
# slow (year+ retention) and modelled separately, if ever.
SHELLFISH_GROUPS = ("mussel", "oyster")

CLOSURE_PPM = 20.0        # CDPH action level: tissue DA >= 20 ppm -> closure

# ordered observed DA categories (censoring absorbed as the non-detect category)
#   1 = non-detect (< LOD)   2 = detected & < 20 ppm   3 = closure (>= 20 ppm)
DA_CATEGORIES = {"non_detect": 1, "detected": 2, "closure": 3}

# --- weekly aggregation of calHABMAP (site x week) ------------------------
# calHABMAP has repeat within-week measurements of the same underlying state;
# we collapse to weekly means to reduce noise and put both markers on a common
# weekly grid. Weeks start Sunday (matches the existing data-acquisition code).
WEEK_RULE = "W-SUN"

# raw calHABMAP columns -> canonical names, averaged within a week.
# GROUND_TRUTH_HAZARD: in-situ truth (cells, water toxins).
GROUND_TRUTH_HAZARD = {
    "pda_ng_ml": "pda_water_ng_ml",
    "tda_ng_ml": "tda_water_ng_ml",
    "pn_total_cells_l": "pn_total_cells_l",                      # computed (deli+seriata)
    "pseudo_nitzschia_seriata_group_cells_l": "pn_seriata_cells_l",
    "alexandrium_spp_cells_l": "alexandrium_cells_l",
    "dinophysis_spp_cells_l": "dinophysis_cells_l",
    "total_phytoplankton_cells_l": "total_phyto_cells_l",
}
# INSITU_PROXY: in-situ physical/biogeochemistry — NOT satellite-derived, so
# these are valid *transition* covariates for the deferred covariate-Q step
# (keeping remote sensing solely in the observation channel).
INSITU_PROXY = {
    "temp_degree_c": "insitu_sst_c",
    "avg_chloro_mg_m3": "insitu_chl_mg_m3",
    "nitrate_um": "nitrate_um",
    "silicate_um": "silicate_um",
}

# the biomass marker used by the world model (calHABMAP PN spp. cells/L)
BIOMASS_COL = "pn_total_cells_l"
