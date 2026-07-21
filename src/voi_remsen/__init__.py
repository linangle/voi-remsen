"""voi-remsen: value of remote sensing for HAB early-closure decisions.

World-model layer (this package): infer the bloom/toxicity dynamics from in-situ
markers — calHABMAP cell counts (biomass) + CDPH tissue domoic acid (toxicity) —
as a continuous-time, multi-state hidden Markov model.

Pipeline:
    data.py            raw markers  -> data/interim/{joint,toxicity}_panel.csv
    R/fit_toxicity_msm.R  toxicity panel -> continuous-time Q reference (K-scan)
    worldmodel/        joint two-marker Bayesian HMM (PyMC)  [next]
    export             posterior Q -> data/output/ (consumed by the decision model)

See docs/methodology.md for the math, decisions, and citations.
"""

from . import config, data, paths  # noqa: F401

__all__ = ["config", "data", "paths"]
