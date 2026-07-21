"""Per-observation log-likelihoods for PSIS-LOO / WAIC model selection.

The fitted likelihood is a single lumped ``pm.Potential`` (the forward algorithm),
so PyMC stores no per-observation log-likelihood and ``az.loo`` / ``az.waic``
cannot run directly. This module reconstructs them.

For an HMM the joint log-likelihood of a sequence factors, via the scaled forward
filter, into one-step-ahead conditional log predictive densities

    log p(y_t | y_{<t}, theta) = log c_t ,

whose sum over t is the exact sequence log-likelihood the sampler used. These are
the natural pointwise terms for a time-series model, and are what LOO/WAIC compare
here (leave-one-observation-out on the one-step predictive densities; see
Buerkner, Gabry & Vehtari 2020 on LFO-CV for the temporal caveat). We replay the
*same* recursion as ``model._forward_loglik``, vectorised over posterior draws,
so the per-draw sum equals the fitted log-likelihood by construction.
"""

from __future__ import annotations

import numpy as np
import xarray as xr
from scipy.special import erfc

from .emissions import SQRT2, SQRT2PI
from .linalg import birthdeath_Q, expm_np


def _draws(idata, name, *shape):
    """Posterior draws of `name` as (S, *shape), chain-major (chain slowest)."""
    x = idata.posterior[name].values           # (chain, draw, *)
    S = x.shape[0] * x.shape[1]
    return x.reshape(S, *shape)


def _emission_S(seq, mu, sigma, E_da):
    """(S,T,K) emission likelihoods, matching emissions.emission_matrix_pt."""
    y = np.nan_to_num(seq["y_cell"], nan=0.0)               # (T,)
    cm, oh, dm = seq["cell_mask"], seq["da_onehot"], seq["da_mask"]
    S, K = mu.shape
    z = (y[None, :, None] - mu[:, None, :]) / sigma[:, None, :]     # (S,T,K)
    dens = np.exp(-0.5 * z * z) / (sigma[:, None, :] * SQRT2PI)
    cens = 0.5 * erfc(-((0.0 - mu) / sigma) / SQRT2)               # (S,K)
    cell = np.where((y > 0)[None, :, None], dens,
                    np.broadcast_to(cens[:, None, :], dens.shape))
    cell = np.where((cm > 0)[None, :, None], cell, 1.0)
    da_obs = np.einsum("tc,skc->stk", oh, E_da)                    # (S,T,K)
    da = np.where((dm > 0)[None, :, None], da_obs, 1.0)
    return cell * da


def _pows_S(P1, gmax):
    """[I, P1, ..., P1^gmax] per draw -> (S, gmax+1, K, K)."""
    S, K, _ = P1.shape
    pows = np.empty((S, gmax + 1, K, K))
    pows[:, 0] = np.eye(K)[None]
    for g in range(1, gmax + 1):
        pows[:, g] = pows[:, g - 1] @ P1
    return pows


def _forward_S(E, gaps, pi0, pows):
    """Vectorised scaled forward filter -> (S,T) one-step log predictive dens."""
    S, T, K = E.shape
    ll = np.empty((S, T))
    a = pi0 * E[:, 0, :]
    c = a.sum(1)
    ll[:, 0] = np.log(c)
    a = a / c[:, None]
    ar = np.arange(S)
    for t in range(1, T):
        Pg = pows[ar, gaps[t]]                      # (S,K,K)
        pred = np.einsum("sk,skj->sj", a, Pg)
        a = pred * E[:, t, :]
        c = a.sum(1)
        ll[:, t] = np.log(c)
        a = a / c[:, None]
    return ll


def _forward_backward(E, gaps, pi0, pows):
    """Scaled forward-backward -> smoothed posterior gamma (S,T,K).

    gamma[s,t,k] = P(S_t = k | y_1..T, theta_s): the state marginal given the
    WHOLE sequence, so a low reading between two blooms is still pulled toward
    'bloom' by its neighbours (what a hard threshold on y_t cannot do).
    """
    S, T, K = E.shape
    ar = np.arange(S)
    alpha = np.empty((S, T, K))
    cs = np.empty((S, T))
    a = pi0 * E[:, 0, :]
    c = a.sum(1)
    alpha[:, 0] = a / c[:, None]
    cs[:, 0] = c
    for t in range(1, T):
        pred = np.einsum("sk,skj->sj", alpha[:, t - 1], pows[ar, gaps[t]])
        a = pred * E[:, t, :]
        c = a.sum(1)
        alpha[:, t] = a / c[:, None]
        cs[:, t] = c

    beta = np.empty((S, T, K))
    beta[:, T - 1, :] = 1.0
    for t in range(T - 2, -1, -1):
        tmp = E[:, t + 1, :] * beta[:, t + 1, :]
        b = np.einsum("skj,sj->sk", pows[ar, gaps[t + 1]], tmp)
        beta[:, t] = b / cs[:, t + 1][:, None]

    g = alpha * beta
    return g / g.sum(2, keepdims=True)


def smoothed_bloom_prob(idata, seqs, K, gap_cap=52):
    """Posterior P(S_t = bloom | data) per observed week, per draw.

    bloom = the top (highest-cell-mean) state, index K-1. Returns
    (bloom, weeks): bloom is (S, n_obs) with a probability per draw; weeks is
    the (n_obs,) datetime64 array aligned to the obs axis (sequence order).
    """
    up = _draws(idata, "up", K - 1)
    down = _draws(idata, "down", K - 1)
    mu = _draws(idata, "mu", K)
    sigma = _draws(idata, "sigma", K)
    E_da = _draws(idata, "E_da", K, 3)
    pi0 = _draws(idata, "pi0", K)
    S = up.shape[0]

    gmax = min(int(max(s["gaps"].max() for s in seqs)), gap_cap)
    P1 = np.stack([expm_np(birthdeath_Q(up[i], down[i])) for i in range(S)])
    pows = _pows_S(P1, gmax)

    blocks, weeks = [], []
    for s in seqs:
        E = _emission_S(s, mu, sigma, E_da)
        gaps = np.minimum(s["gaps"], gmax).astype(int)
        gamma = _forward_backward(E, gaps, pi0, pows)      # (S,T,K)
        blocks.append(gamma[:, :, K - 1])                  # bloom = top state
        weeks.append(s["week_start"])
    return np.concatenate(blocks, axis=1), np.concatenate(weeks)


def pointwise_loglik(idata, seqs, K, gap_cap=52):
    """Attach idata.log_likelihood['obs'] of shape (chain, draw, n_obs).

    n_obs is the total number of observed weeks across sequences (concatenated
    in `seqs` order). Returns the same idata.
    """
    nchain = idata.posterior.sizes["chain"]
    ndraw = idata.posterior.sizes["draw"]
    up = _draws(idata, "up", K - 1)
    down = _draws(idata, "down", K - 1)
    mu = _draws(idata, "mu", K)
    sigma = _draws(idata, "sigma", K)
    E_da = _draws(idata, "E_da", K, 3)
    pi0 = _draws(idata, "pi0", K)
    S = up.shape[0]

    gmax = min(int(max(s["gaps"].max() for s in seqs)), gap_cap)
    P1 = np.stack([expm_np(birthdeath_Q(up[i], down[i])) for i in range(S)])
    pows = _pows_S(P1, gmax)

    blocks = []
    for s in seqs:
        E = _emission_S(s, mu, sigma, E_da)
        gaps = np.minimum(s["gaps"], gmax).astype(int)
        blocks.append(_forward_S(E, gaps, pi0, pows))       # (S, T_s)
    ll = np.concatenate(blocks, axis=1)                     # (S, n_obs)

    da = xr.DataArray(
        ll.reshape(nchain, ndraw, ll.shape[1]),
        dims=("chain", "draw", "obs"),
        coords={"chain": idata.posterior.chain, "draw": idata.posterior.draw},
    )
    ds = xr.Dataset({"obs": da})
    if hasattr(idata, "add_groups"):        # classic arviz.InferenceData
        idata.add_groups({"log_likelihood": ds})
    else:                                   # arviz 1.x DataTree
        idata["log_likelihood"] = ds
    return idata
