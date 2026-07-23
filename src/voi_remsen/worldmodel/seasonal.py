"""Phase D: non-stationary seasonal generator Q(day-of-year).

The birth-death rates become smooth functions of week-of-year via a log link on
low-order Fourier (harmonic) terms, so the generator and its one-week transition
P_t = expm(Q(woy)) breathe with the season. Everything else -- ordered states,
the two censored emissions -- is unchanged from the stationary model; only the
dynamics are seasonal (time-of-year only, no physical covariates).

Because Q now varies week to week, the transition over an unobserved gap is the
PRODUCT of the intervening weeks' matrices, not a power of one matrix. So the
forward filter steps over the full weekly grid (observed + unobserved weeks),
applying P_t at each week and multiplying in an emission only on observed weeks.

There are only 53 distinct week-of-year bins, so we exponentiate 53 generators
once (batched) and index them by each grid week's bin.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .emissions import emission_matrix_pt
from .linalg import birthdeath_Q, expm_np, expm_pt

WOY = 53                       # week-of-year bins (0..52)
PERIOD = 365.25 / 7.0          # weeks per year


def season_design(woy, H, period=PERIOD):
    """Harmonic design matrix: [1, cos, sin, cos2, sin2, ...] -> (n, 1+2H)."""
    woy = np.asarray(woy, float)
    cols = [np.ones_like(woy)]
    for h in range(1, H + 1):
        ang = 2.0 * np.pi * h * woy / period
        cols += [np.cos(ang), np.sin(ang)]
    return np.stack(cols, axis=1)


def _woy_of(dates):
    """Week-of-year bin (0..52) for a datetime64 array."""
    doy = pd.to_datetime(dates).dayofyear.to_numpy()
    return np.minimum((doy - 1) // 7, WOY - 1).astype(int)


def build_grid_sequences(panel, sites=None):
    """Per-site sequences expanded onto the full weekly grid.

    Each grid week carries its week-of-year bin and emission arrays (masks are 0
    on unobserved weeks, so their emission factor is 1). Observed values are
    placed at their grid positions.
    """
    if sites is None:
        sites = sorted(panel["location_code"].dropna().unique())
    gseqs = []
    for code in sites:
        d = panel[panel["location_code"] == code].copy()
        d = d.dropna(subset=["cells", "da"], how="all").sort_values("week_start")
        if len(d) < 2:
            continue
        wk = pd.to_datetime(d["week_start"])
        widx = ((wk - wk.min()).dt.days // 7).to_numpy()
        n_grid = int(widx.max()) + 1
        grid_dates = wk.min() + pd.to_timedelta(np.arange(n_grid) * 7, unit="D")

        y_cell = np.zeros(n_grid)
        cell_mask = np.zeros(n_grid)
        da_onehot = np.zeros((n_grid, 3))
        da_mask = np.zeros(n_grid)

        obs_cellmask = d["cells"].notna().to_numpy().astype(float)
        y_cell[widx] = np.where(obs_cellmask > 0,
                                np.nan_to_num(d["cells_log10"].to_numpy(), nan=0.0), 0.0)
        cell_mask[widx] = obs_cellmask
        dm = d["da"].notna().to_numpy().astype(float)
        cat = np.nan_to_num(d["da_cat"].to_numpy(), nan=1).astype(int)
        da_onehot[widx, cat - 1] = dm
        da_mask[widx] = dm

        gseqs.append(dict(code=code, n=len(d), n_grid=n_grid,
                          woy=_woy_of(grid_dates), y_cell=y_cell,
                          cell_mask=cell_mask, da_onehot=da_onehot,
                          da_mask=da_mask, obs_widx=widx,
                          week_start=grid_dates.to_numpy()))
    return gseqs


def _build_Q_batched(up, down, K):
    """Batched birth-death generator: up/down (B,K-1) -> Q (B,K,K)."""
    import pytensor.tensor as pt
    B = up.shape[0]
    Q = pt.zeros((B, K, K))
    for i in range(K - 1):
        Q = pt.set_subtensor(Q[:, i, i + 1], up[:, i])
        Q = pt.set_subtensor(Q[:, i + 1, i], down[:, i])
    rowsum = Q.sum(axis=2)
    for k in range(K):
        Q = pt.set_subtensor(Q[:, k, k], -rowsum[:, k])
    return Q


def _forward_seasonal(E, woy, pi0, P_woy):
    """Forward filter over the full grid with a per-week transition matrix."""
    import pytensor
    import pytensor.tensor as pt
    tiny = 1e-30                     # guard against all-emissions underflow -> log(0)
    a0 = pi0 * E[0]
    c0 = pt.maximum(a0.sum(), tiny)
    a0 = a0 / c0

    def step(E_t, woy_t, a_prev, P_):
        pred = a_prev.dot(P_[woy_t])
        a = pred * E_t
        c = pt.maximum(a.sum(), tiny)
        return a / c, pt.log(c)

    (_, logc), _ = pytensor.scan(
        step, sequences=[E[1:], woy[1:]],
        outputs_info=[a0, None], non_sequences=[P_woy])
    return pt.log(c0) + logc.sum()


def build_seasonal_model(gseqs, K, H=1, mu_prior=None, sigma_floor=0.15,
                         mu_sep=0.5, amp_sd=0.7):
    """PyMC model with seasonal (harmonic) birth-death rates for K states.

    H is the number of annual harmonics on each log-rate. amp_sd is the prior sd
    on the harmonic coefficients (the seasonal amplitude). Log-rate intercepts use
    a neutral weakly-informative prior N(log 0.1, 1) -- a generic weekly-transition
    scale, NOT the stationary-fit posterior means (which would double-use the same
    data as an empirical-Bayes anchor and understate uncertainty).
    """
    import pymc as pm
    import pytensor.tensor as pt

    if mu_prior is None:
        mu_prior = np.linspace(0.5, 4.0, K)
    Xd = season_design(np.arange(WOY), H)          # (53, 1+2H)
    ncoef = Xd.shape[1]
    a_int, sd_int = np.log(0.1), 1.0                # neutral intercept prior

    with pm.Model() as m:
        # harmonic coefficients: column 0 = intercept, rest = seasonal terms
        mu_c = np.zeros((K - 1, ncoef))
        sd_c = np.full((K - 1, ncoef), amp_sd)
        mu_c[:, 0], sd_c[:, 0] = a_int, sd_int
        coef_up = pm.Normal("coef_up", mu=mu_c, sigma=sd_c, shape=(K - 1, ncoef))
        coef_down = pm.Normal("coef_down", mu=mu_c, sigma=sd_c, shape=(K - 1, ncoef))

        Xt = pt.as_tensor_variable(Xd)
        up_woy = pt.exp(Xt @ coef_up.T)            # (53, K-1)
        down_woy = pt.exp(Xt @ coef_down.T)
        Q_woy = _build_Q_batched(up_woy, down_woy, K)
        P_woy = expm_pt(Q_woy)                      # (53, K, K)

        mu0 = pm.TruncatedNormal("mu0", mu=float(mu_prior[0]), sigma=1.0,
                                 lower=0.0, initval=float(max(mu_prior[0], 0.2)))
        dmu = pm.HalfNormal("dmu", sigma=1.5, shape=K - 1,
                            initval=np.maximum(np.diff(mu_prior) - mu_sep, 0.3))
        mu = pm.Deterministic("mu", mu0 + pt.concatenate(
            [pt.zeros(1), pt.cumsum(mu_sep + dmu)]))
        sigma = pm.Deterministic(
            "sigma", sigma_floor + pm.HalfNormal("sigma_raw", 1.0))
        pi_zero = pm.Beta("pi_zero", 1.0, 1.0, shape=K)
        E_da = pm.Dirichlet("E_da", a=np.ones((K, 3)), shape=(K, 3))
        pi0 = pm.Dirichlet("pi0", a=np.ones(K))

        total = pt.as_tensor_variable(0.0)
        for s in gseqs:
            E = emission_matrix_pt(
                pt.as_tensor_variable(s["y_cell"]),
                pt.as_tensor_variable(s["cell_mask"]),
                pt.as_tensor_variable(s["da_onehot"]),
                pt.as_tensor_variable(s["da_mask"]),
                mu, sigma, pi_zero, E_da)
            total = total + _forward_seasonal(
                E, pt.as_tensor_variable(s["woy"].astype("int64")), pi0, P_woy)
        pm.Potential("loglik", total)
    return m


def seasonal_P_woy(coef_up, coef_down, K, H):
    """Numpy: 53 weekly transition matrices from harmonic coefficients."""
    Xd = season_design(np.arange(WOY), H)
    up = np.exp(Xd @ np.atleast_2d(coef_up).T)      # (53, K-1)
    down = np.exp(Xd @ np.atleast_2d(coef_down).T)
    return np.stack([expm_np(birthdeath_Q(up[w], down[w])) for w in range(WOY)])


def periodic_occupancy(P_woy, iters=80):
    """Propagated periodic prevalence of a time-varying weekly chain.

    Returns pi (53, K): pi[w,k] = P(state k in week-of-year w) in the recurring
    annual cycle, i.e. the fixed point of pi_w = pi_{w-1} P_woy[w] over the year.
    This is the TRUE seasonal occupancy -- unlike the frozen-rate ratio
    u(w)/(u(w)+d(w)), which is only the equilibrium week w would reach if its
    rates never changed, and which leads the propagated prevalence by weeks.

    P_woy may be (53,K,K) or batched (S,53,K,K); the batch axis is preserved.
    """
    P = np.asarray(P_woy)
    batched = P.ndim == 4
    if not batched:
        P = P[None]
    S, W, K, _ = P.shape
    pi = np.full((S, K), 1.0 / K)
    for _ in range(iters):                      # converge the week-0 distribution
        p = pi
        for w in range(W):
            p = np.einsum("sk,skj->sj", p, P[:, w])
        pi = p
    out = np.empty((S, W, K))
    p = pi
    for w in range(W):
        p = np.einsum("sk,skj->sj", p, P[:, w])
        out[:, w] = p
    return out if batched else out[0]


_DA_MAG = {1: 0.5, 2: 5.0, 3: 30.0}


def simulate_seasonal_panel(coef_up, coef_down, mu, sigma, pi_zero, E_da, H,
                            n_weeks=800, n_sites=4, p_keep=0.65, p_cell=0.7,
                            p_da=0.5, seed=0, base="2008-01-06"):
    """Ground-truth generator for the seasonal model (validation only).

    Steps a weekly Markov chain with the seasonal P_woy, emits the two markers
    (hurdle cells + categorical DA), and drops weeks to create gaps. sigma is a
    shared scalar; pi_zero is (K,). Returns a panel shaped like joint_panel.
    """
    mu, pi_zero, E_da = map(np.asarray, (mu, pi_zero, E_da))
    sigma = float(sigma)
    K = len(mu)
    P_woy = seasonal_P_woy(coef_up, coef_down, K, H)
    rng = np.random.default_rng(seed)
    dates = pd.Timestamp(base) + pd.to_timedelta(np.arange(n_weeks) * 7, unit="D")
    woy = _woy_of(dates)
    rows = []
    for site in range(n_sites):
        s = rng.integers(K)
        for w in range(n_weeks):
            if w > 0:
                s = rng.choice(K, p=P_woy[woy[w], s])
            if rng.random() >= p_keep:
                continue
            cells = clog = da = dacat = np.nan
            if rng.random() < p_cell:
                if rng.random() < pi_zero[s]:
                    clog = 0.0; cells = 0.0
                else:
                    ystar = rng.normal(mu[s], sigma)
                    while ystar <= 0:
                        ystar = rng.normal(mu[s], sigma)
                    clog = ystar; cells = 10.0 ** clog - 1.0
            if rng.random() < p_da:
                dacat = int(rng.choice(3, p=E_da[s]) + 1)
                da = _DA_MAG[dacat]
            if np.isnan(cells) and np.isnan(da):
                continue
            rows.append((f"S{site}", dates[w], cells, clog, da, dacat))
    return pd.DataFrame(rows, columns=["location_code", "week_start", "cells",
                                       "cells_log10", "da", "da_cat"])
