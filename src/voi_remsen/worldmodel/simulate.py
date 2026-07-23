"""Simulate joint two-marker data from a known Q + emissions.

Used to check that `build_model` recovers a known generator before it is trusted
on the real panel. Produces a joint-panel-shaped DataFrame (so it flows through
`build_sequences` unchanged), with whole weeks randomly dropped to create the
irregular gaps the continuous-time model is meant to handle.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .linalg import expm_np

_DA_MAG = {1: 0.5, 2: 5.0, 3: 30.0}   # nominal ppm per category (magnitude only)


def simulate_panel(Q_week, mu, sigma, pi_zero, E_da, n_weeks=400, n_sites=3,
                   p_keep=0.6, p_cell=0.7, p_da=0.5, seed=0):
    """sigma is a shared scalar; pi_zero is (K,) per-state cell non-detect rate."""
    mu, pi_zero, E_da = map(np.asarray, (mu, pi_zero, E_da))
    sigma = float(sigma)
    P1 = expm_np(Q_week)
    K = len(mu)
    rng = np.random.default_rng(seed)
    base = pd.Timestamp("2010-01-03")   # a Sunday
    rows = []
    for site in range(n_sites):
        s = rng.integers(K)
        for w in range(n_weeks):
            if w > 0:
                s = rng.choice(K, p=P1[s])
            if rng.random() > p_keep:
                continue                       # whole week unobserved -> gap
            cells = clog = da = dacat = np.nan
            if rng.random() < p_cell:
                if rng.random() < pi_zero[s]:
                    clog = 0.0; cells = 0.0    # hurdle: exact non-detect
                else:
                    ystar = rng.normal(mu[s], sigma)   # truncated at 0
                    while ystar <= 0:
                        ystar = rng.normal(mu[s], sigma)
                    clog = ystar; cells = 10.0 ** clog - 1.0
            if rng.random() < p_da:
                dacat = int(rng.choice(3, p=E_da[s]) + 1)
                da = _DA_MAG[dacat]
            if np.isnan(cells) and np.isnan(da):
                continue
            rows.append((f"S{site}", base + pd.Timedelta(weeks=w),
                         cells, clog, da, dacat))
    return pd.DataFrame(
        rows, columns=["location_code", "week_start",
                       "cells", "cells_log10", "da", "da_cat"])
