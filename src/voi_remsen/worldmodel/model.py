"""Joint two-marker continuous-time HMM (PyMC).

One ordered latent severity S(t) in {0..K-1} drives both markers. Transitions are
a continuous-time birth-death process with generator Q (per-week units); the
discrete states are marginalized by the forward algorithm using the per-gap
transition P(gap) = expm(Q)^gap (gaps are whole weeks on the weekly grid).

Identifiability anchor: ordered cell-emission means mu[0] < ... < mu[K-1].
Cross-checked against the DA-only msm reference (K=2).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .emissions import emission_matrix_pt
from .linalg import expm_pt


# data prep
def build_sequences(panel: pd.DataFrame, sites=None):
    """Split the joint panel into per-site sequences of observed weeks.

    Returns a list of dicts with numpy arrays: y_cell (log10, NaN where absent),
    cell_mask, da_onehot (T,3), da_mask, and gaps (whole weeks to previous obs).
    """
    if sites is None:
        sites = sorted(panel["location_code"].dropna().unique())
    seqs = []
    for code in sites:
        d = panel[panel["location_code"] == code].copy()
        d = d.dropna(subset=["cells", "da"], how="all")
        d = d.sort_values("week_start")
        if len(d) < 2:
            continue
        wk = pd.to_datetime(d["week_start"])
        widx = ((wk - wk.min()).dt.days // 7).to_numpy()
        gaps = np.diff(widx, prepend=widx[0]).astype(int)   # gaps[0]=0

        y_cell = np.where(d["cells"].notna().to_numpy(),
                          d["cells_log10"].to_numpy(), np.nan)
        cell_mask = d["cells"].notna().to_numpy().astype(float)

        da_mask = d["da"].notna().to_numpy().astype(float)
        cat = np.nan_to_num(d["da_cat"].to_numpy(), nan=1).astype(int)  # 1..3
        onehot = np.zeros((len(d), 3))
        onehot[np.arange(len(d)), cat - 1] = da_mask     # zero row where missing

        seqs.append(dict(code=code, y_cell=y_cell, cell_mask=cell_mask,
                         da_onehot=onehot, da_mask=da_mask, gaps=gaps,
                         n=len(d), week_start=wk.to_numpy()))
    return seqs


# helpers
def _build_Q(up, down, K):
    import pytensor.tensor as pt
    Q = pt.zeros((K, K))
    for i in range(K - 1):
        Q = pt.set_subtensor(Q[i, i + 1], up[i])
    for i in range(1, K):
        Q = pt.set_subtensor(Q[i, i - 1], down[i - 1])
    return Q - pt.diag(Q.sum(axis=1))          # diagonal = -offdiag rowsum


def _matrix_powers(P1, gmax, K):
    """Stack [I, P1, P1^2, ..., P1^gmax] -> (gmax+1, K, K).

    gmax is a static Python int (max whole-week gap), so this unrolls to plain
    matmuls -- no inner scan (a nested scan confuses PyMC's logp compilation).
    """
    import pytensor.tensor as pt
    mats = [pt.eye(K)]
    for _ in range(gmax):
        mats.append(mats[-1] @ P1)
    return pt.stack(mats)


def _forward_loglik(E, gaps, pi0, pows):
    """Scaled forward algorithm for one sequence; returns scalar log-likelihood.

    E    (T,K) emission likelihoods; gaps (T,) whole-week gaps as an int tensor
    (gaps[0] unused); pi0 (K,) initial distribution; pows (gmax+1,K,K) powers.

    `pows` is passed as an explicit scan non-sequence and E/gaps as sequences, so
    the scan's inner graph is pure arithmetic -- no RV/RNG lands inside the Scan
    Op (which PyMC's logp compilation cannot rewrite).
    """
    import pytensor, pytensor.tensor as pt
    a0 = pi0 * E[0]
    c0 = a0.sum()
    a0 = a0 / c0

    def step(E_t, gap_t, a_prev, pows_):
        pred = a_prev.dot(pows_[gap_t])
        a = pred * E_t
        c = a.sum()
        return a / c, pt.log(c)

    (_, logc), _ = pytensor.scan(
        step,
        sequences=[E[1:], gaps[1:]],
        outputs_info=[a0, None],
        non_sequences=[pows])
    return pt.log(c0) + logc.sum()


# model
def build_model(seqs, K, mu_prior=None, sigma_floor=0.15, rate_prior=1.0,
                gap_cap=52, mu_sep=0.5):
    """Assemble the PyMC model over the given sequences for K latent states.

    gap_cap bounds how many whole-week transition powers P1^g are built. Blooms
    mix on a scale of weeks, so P1^g is effectively stationary well before a year;
    capping gaps at ~52 weeks is numerically negligible but keeps the graph small
    (a deep matmul chain is what makes compilation slow).

    mu_sep is a minimum gap (log10 cells) between consecutive state means. Without
    it, at K>=3 a state can collapse onto its neighbour (same cell mean, weakly
    told apart by DA) -- a degenerate mode that leaves chains stuck at different
    optima (r-hat >> 1). Two states at the same biomass level aren't identifiable
    from the cell marker, so enforcing a floor separation is the right regulariser;
    it is slack for the well-separated K=2 fit.
    """
    import pymc as pm
    import pytensor.tensor as pt

    gmax = min(int(max(s["gaps"].max() for s in seqs)), gap_cap)
    if mu_prior is None:
        mu_prior = np.linspace(0.5, 4.0, K)          # log10 cells anchors

    with pm.Model() as m:
        up = pm.Exponential("up", rate_prior, shape=K - 1)
        down = pm.Exponential("down", rate_prior, shape=K - 1)
        Q = _build_Q(up, down, K)
        P1 = expm_pt(Q)
        pows = _matrix_powers(P1, gmax, K)

        # Ordered, NON-NEGATIVE cell means built from positive increments:
        # y = log10(cells+1) >= 0, and a near-all-censored low state otherwise
        # lets mu -> -inf (the tobit degeneracy that wrecks convergence).
        mu0 = pm.TruncatedNormal("mu0", mu=float(mu_prior[0]), sigma=1.0,
                                 lower=0.0, initval=float(max(mu_prior[0], 0.2)))
        dmu = pm.HalfNormal("dmu", sigma=1.5, shape=K - 1,
                            initval=np.maximum(np.diff(mu_prior) - mu_sep, 0.3))
        gaps_mu = mu_sep + dmu                     # each gap >= mu_sep
        mu = pm.Deterministic("mu", mu0 + pt.concatenate([pt.zeros(1),
                                                          pt.cumsum(gaps_mu)]))
        sigma = pm.Deterministic(
            "sigma", sigma_floor + pm.HalfNormal("sigma_raw", 1.0, shape=K))
        E_da = pm.Dirichlet("E_da", a=np.ones((K, 3)), shape=(K, 3))
        pi0 = pm.Dirichlet("pi0", a=np.ones(K))

        total = pt.as_tensor_variable(0.0)
        for s in seqs:
            E = emission_matrix_pt(
                pt.as_tensor_variable(np.nan_to_num(s["y_cell"], nan=0.0)),
                pt.as_tensor_variable(s["cell_mask"]),
                pt.as_tensor_variable(s["da_onehot"]),
                pt.as_tensor_variable(s["da_mask"]),
                mu, sigma, E_da)
            gaps = np.minimum(s["gaps"], gmax).astype("int64")
            total = total + _forward_loglik(
                E, pt.as_tensor_variable(gaps), pi0, pows)
        pm.Potential("loglik", total)
    return m
