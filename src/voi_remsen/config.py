CALHABMAP_FILE = "calhabmap.csv"
CDPH_FILE = "CDPH_tissueDA.xlsx"

CO_LOCATED_SITES = {
    "SCW": "Santa Cruz Wharf",
    "CPP": "San Luis Obispo, Cal Poly Pier",
    "MW":  "Monterey Bay, Commercial Wharf",
    "SIO": "La Jolla, Scripps Pier",
}

PN_COLUMNS = [
    "pseudo_nitzschia_delicatissima_group_cells_l",
    "pseudo_nitzschia_seriata_group_cells_l",
]

SHELLFISH_GROUPS = ("mussel", "oyster")

CLOSURE_PPM = 20.0

DA_CATEGORIES = {"non_detect": 1, "detected": 2, "closure": 3}

WEEK_RULE = "W-SUN"

GROUND_TRUTH_HAZARD = {
    "pda_ng_ml": "pda_water_ng_ml",
    "tda_ng_ml": "tda_water_ng_ml",
    "pn_total_cells_l": "pn_total_cells_l",
    "pseudo_nitzschia_seriata_group_cells_l": "pn_seriata_cells_l",
    "alexandrium_spp_cells_l": "alexandrium_cells_l",
    "dinophysis_spp_cells_l": "dinophysis_cells_l",
    "total_phytoplankton_cells_l": "total_phyto_cells_l",
}
INSITU_PROXY = {
    "temp_degree_c": "insitu_sst_c",
    "avg_chloro_mg_m3": "insitu_chl_mg_m3",
    "nitrate_um": "nitrate_um",
    "silicate_um": "silicate_um",
}

BIOMASS_COL = "pn_total_cells_l"
