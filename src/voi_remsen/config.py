CALHABMAP_FILE = "calhabmap.csv"
CDPH_FILE = "CDPH_tissueDA.xlsx"

CO_LOCATED_SITES = {
    "SCW": "Santa Cruz Wharf",
    "CPP": "San Luis Obispo, Cal Poly Pier",
    "MW":  "Monterey Bay, Commercial Wharf",
    "SIO": "La Jolla, Scripps Pier",
}

# The world model is anchored on a SINGLE site. Pooling the four sites shared one
# generator/emission across very different baselines and sampling regimes, so the
# latent states partly encoded site identity (at K=2 the state was nearly a proxy
# for the site) and toxicity events were too sparse elsewhere to identify a hazard
# state. Santa Cruz Wharf carries 22 of the 34 closure-level tissue observations,
# 663 tissue weeks, and near-balanced marker coverage (667 cell vs 663 DA weeks),
# so the joint state is informed by toxicity rather than dominated by biomass.
# NOTE: SCW cell counts are the *seriata* group only (delicatissima is not
# reported there) -- internally consistent, but not comparable to other sites.
ANCHOR_SITE = "SCW"
MODEL_SITES = [ANCHOR_SITE]

# Tissue samples are matched to calHABMAP stations by COORDINATES (nearest station
# within this radius), not by exact `Sample Site` string. Name matching silently
# dropped tissue records whose label never equals a station name (e.g. sites
# recorded under a county), understating coverage.
MATCH_RADIUS_KM = 3.0

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
PDA_COL = "pda_water_ng_ml"   # particulate DA: observes biomass x toxin quota
