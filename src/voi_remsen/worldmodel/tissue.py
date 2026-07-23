"""One-compartment shellfish tissue-burden model H_t (Mafra et al. 2010).

Tissue domoic acid is a STOCK: produced by ingesting particulate toxin and
removed by approximately first-order elimination,

    dH/dt = alpha * E(t) - k * H(t),

with E(t) = particulate toxin exposure (pDA, or latent B_t*Q_t) and k the
elimination rate. Holding E constant across an interval of length D (days) gives
the exact update actually used here:

    H_{t+D} = H_t exp(-k D) + (alpha E_t / k) (1 - exp(-k D)).

Mafra et al. (2010) support this one-compartment form: their two-compartment
(visceral/non-visceral transfer) extension improved r^2 by only 1-9%, and whole
tissue fits were good (r^2 0.66-0.94 for mussels). They measured M. edulis at
k = 1.4-1.6 /day and cite M. californianus at 0.3-0.5 /day (Whyte et al. 1995).
Santa Cruz tissue is 702/706 M. californianus, so this is a California-mussel
model -- but k is ESTIMATED with a broad prior rather than fixed at a lab value,
since field depuration may be slower than in clean-diet laboratory conditions.

Deliberately excluded (per the paper and the available data): visceral/non-visceral
compartments (not resolved by CDPH whole-tissue assays), and body size (Mafra found
no clear mussel size effect and CDPH records no sizes here).

Observations are log-normal with LEFT CENSORING at each sample's own reported
detection limit, so non-detects contribute P(Y < LOD) rather than a fabricated value.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

H_FLOOR = 1e-3          # keeps log H finite when exposure has been 0 for a long time


def build_exact_series(site="SCW", radius_km=3.0):
    """Exact-date water (exposure) and tissue series for one site.

    Weekly summaries are NOT used: the kinetic model needs real elapsed days.
    """
    from ..data import load_calhabmap_raw, load_cdph_raw

    w = load_calhabmap_raw()
    w = w[w["location_code"] == site][["time_utc", "pn_total_cells_l", "pda_ng_ml"]]
    w = w.dropna(subset=["time_utc"]).copy()
    w["date"] = pd.to_datetime(w["time_utc"]).dt.normalize()
    w = (w.groupby("date")[["pn_total_cells_l", "pda_ng_ml"]].mean()
           .reset_index().sort_values("date"))
    w = w.rename(columns={"pn_total_cells_l": "cells", "pda_ng_ml": "pda"})

    t = load_cdph_raw(radius_km=radius_km)
    t = t[t["location_code"] == site][["date", "da", "cens", "lod"]].copy()
    t["date"] = pd.to_datetime(t["date"]).dt.normalize()
    t = t.groupby("date").agg(da=("da", "max"), cens=("cens", "min"),
                              lod=("lod", "max")).reset_index().sort_values("date")
    # the kinetic model needs an exposure history, so keep the water-covered era
    t = t[t["date"] >= w["date"].min()].reset_index(drop=True)
    return w.reset_index(drop=True), t


def build_timeline(water, tissue):
    """Merged event timeline with elapsed days and the exposure active on each step.

    Returns dt (days between events), E_prev (exposure held over that interval,
    forward-filled from the last water sample) and the event index of each tissue
    observation.
    """
    dates = np.union1d(water["date"].values, tissue["date"].values)
    dates = np.sort(dates)
    dt = np.diff(dates).astype("timedelta64[D]").astype(float)
    dt = np.concatenate([[0.0], dt])

    # exposure carried forward from the most recent water sample at/or before t
    wd = water["date"].values
    E_at = pd.Series(np.nan, index=range(len(dates)), dtype=float)
    pos = np.searchsorted(dates, wd)
    E_at.iloc[pos] = water["pda"].to_numpy()
    E = E_at.ffill().fillna(0.0).to_numpy()
    E_prev = np.concatenate([[0.0], E[:-1]])          # exposure over interval j-1 -> j

    tidx = np.searchsorted(dates, tissue["date"].values)
    return dict(dates=dates, dt=dt, E_prev=E_prev, tissue_idx=tidx,
                n_events=len(dates))


def kinetic_H_np(dt, E_prev, alpha, k, H0):
    """Numpy one-compartment recursion (simulation / checking)."""
    H = np.empty(len(dt))
    h = H0
    for j in range(len(dt)):
        if dt[j] > 0:
            e = np.exp(-k * dt[j])
            h = h * e + (alpha * E_prev[j] / k) * (1.0 - e)
        H[j] = h
    return H


def _censored_lognormal_logp(y, cens, lod, H, sigma):
    """log p(tissue | H): lognormal, left-censored at each sample's own LOD."""
    import pytensor.tensor as pt
    logH = pt.log(pt.maximum(H, H_FLOOR))
    z_obs = (pt.log(pt.maximum(y, 1e-6)) - logH) / sigma
    lp_obs = -0.5 * z_obs ** 2 - pt.log(sigma) - 0.5 * np.log(2 * np.pi)
    z_cen = (np.log(np.maximum(lod, 1e-6)) - logH) / sigma
    lp_cen = pt.log(pt.maximum(0.5 * pt.erfc(-z_cen / np.sqrt(2.0)), 1e-12))
    return pt.switch(cens > 0, lp_cen, lp_obs)


def build_tissue_model(tl, tissue, kinetic=True, k_med=0.4, k_sd=0.8, n_fit=None):
    """PyMC tissue model. kinetic=True -> stock with memory; False -> no-memory.

    The no-memory comparator sets H proportional to CURRENT exposure (the k->inf
    limit of the same kinetics), which is what the previous ordered-logit did.

    n_fit restricts the LIKELIHOOD to the first n_fit tissue observations (for
    leave-future-out refits) while `obs_logp` is still reported for all of them.
    The exposure path is exogenous, so conditioning on y_{1:n} only changes which
    tissue observations inform the parameters -- exactly what LFO-CV requires.
    """
    import pymc as pm
    import pytensor
    import pytensor.tensor as pt

    y = tissue["da"].to_numpy(float)
    cens = tissue["cens"].to_numpy(float)
    lod = np.nan_to_num(tissue["lod"].to_numpy(float), nan=2.5)
    dt = tl["dt"]; E_prev = tl["E_prev"]; tidx = tl["tissue_idx"]

    with pm.Model() as m:
        log_alpha = pm.Normal("log_alpha", np.log(4.0), 2.0)
        alpha = pm.Deterministic("alpha", pt.exp(log_alpha))
        sigma = pm.HalfNormal("sigma_tissue", 1.5)
        H_bg = pm.HalfNormal("H_bg", 1.0)          # analytical/background baseline

        if kinetic:
            # broad prior: lab M. californianus ~0.3-0.5/day, but field may differ
            k = pm.Lognormal("k", np.log(k_med), k_sd)
            pm.Deterministic("half_life_days", np.log(2.0) / k)
            pm.Deterministic("weekly_retention", pt.exp(-7.0 * k))
            H0 = pm.HalfNormal("H0", 1.0)

            def step(dt_j, E_j, h_prev, k_, a_):
                e = pt.exp(-k_ * dt_j)
                return pt.switch(dt_j > 0, h_prev * e + (a_ * E_j / k_) * (1.0 - e),
                                 h_prev)

            H_path, _ = pytensor.scan(
                step, sequences=[pt.as_tensor_variable(dt),
                                 pt.as_tensor_variable(E_prev)],
                outputs_info=[H0], non_sequences=[k, alpha])
            H_obs = H_bg + H_path[tidx]
        else:
            E_now = pt.as_tensor_variable(E_prev)[tidx]
            H_obs = H_bg + alpha * E_now

        lp = _censored_lognormal_logp(pt.as_tensor_variable(y),
                                      pt.as_tensor_variable(cens),
                                      lod, H_obs, sigma)
        pm.Deterministic("obs_logp", lp)           # per-observation, for LFO
        n = len(y) if n_fit is None else int(n_fit)
        pm.Potential("tissue_lik", lp[:n].sum())   # condition on y_{1:n} only
    return m
