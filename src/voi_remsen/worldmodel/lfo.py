"""Approximate leave-future-out cross-validation (PSIS-LFO-CV).

Buerkner, Gabry & Vehtari (2020), "Approximate leave-future-out cross-validation
for Bayesian time series models". Ordinary PSIS-LOO is POSITIVELY BIASED for time
series -- it lets the posterior see the future -- so it systematically overstates
forecast skill. LFO-CV instead targets

    ELPD_LFO = sum_i log p(y_{i+1} | y_{1:i}),

conditioning both the predictive and the parameter posterior on the past only.

Forward algorithm (their Section 2.1.2): refit at i* = L, then walk forward
reusing that posterior with importance weights

    r_i^(s)  =  prod_{j in (i*+1):i}  p(y_j | y_{1:(j-1)}, theta^(s)),

stabilised by Pareto-smoothed importance sampling. When the PSIS shape parameter
k exceeds tau (default 0.7) the weights have degenerated, so the model is refit at
that point and the process restarts. In their simulations refits were needed for
only ~1-3% of observations.
"""

from __future__ import annotations

import numpy as np


def psis_smooth(logw):
    """Pareto-smooth log importance weights; returns (log weights, Pareto k)."""
    logw = np.asarray(logw, float) - np.max(logw)
    w = np.exp(logw)
    S = len(w)
    M = int(min(0.2 * S, 3 * np.sqrt(S)))          # tail size (Vehtari et al.)
    if M < 5:
        return logw, 0.0
    order = np.argsort(w)
    tail_idx = order[-M:]                          # ascending in w
    mu = w[order[-(M + 1)]]                        # cutoff just below the tail
    xs = np.maximum(w[tail_idx] - mu, 1e-300)      # already ascending
    if xs[-1] <= 1e-300:
        return logw, 0.0
    # Zhang & Stephens (2009) empirical-Bayes GPD estimator
    n = len(xs)
    m_grid = 30 + int(np.sqrt(n))
    bs = (1.0 - np.sqrt(m_grid / (np.arange(1, m_grid + 1) - 0.5)))
    bs = bs / (3.0 * xs[max(int(n / 4 + 0.5) - 1, 0)]) + 1.0 / xs[-1]
    kk = np.array([np.mean(np.log1p(-b * xs)) for b in bs])
    Ls = n * (np.log(np.maximum(-bs / kk, 1e-300)) - kk - 1.0)
    wts = 1.0 / np.array([np.sum(np.exp(Ls - Lj)) for Lj in Ls])
    b = float(np.sum(bs * wts))
    k = float(np.mean(np.log1p(-b * xs)))
    sigma = -k / b if b != 0 else 0.0
    if sigma > 0 and abs(k) > 1e-10:
        p = (np.arange(1, M + 1) - 0.5) / M        # ascending quantiles
        qs = mu + sigma * np.expm1(-k * np.log1p(-p)) / k
        w[tail_idx] = np.minimum(qs, np.max(w))    # both ascending -> aligned
    out = np.log(np.maximum(w, 1e-300))
    return out - np.max(out), k


def psis_lfo(logp, L=None, tau=0.7, refit=None, M=1):
    """PSIS-LFO-CV for a factorised model.

    logp   : (S, N) per-draw, per-observation log predictive densities
             log p(y_i | y_{1:i-1}, theta^(s)) under the CURRENT posterior.
    refit  : callable(i) -> (S, N) logp recomputed from a posterior fit to y_{1:i}.
             If None, the initial posterior is reused throughout (no refits;
             report the returned `n_refit`/k diagnostics with that caveat).
    Returns dict with elpd, pointwise contributions, Pareto k's and refit points.
    """
    logp = np.asarray(logp, float)
    S, N = logp.shape
    L = int(0.3 * N) if L is None else L
    i_star = L
    cur = logp
    elpd_i, ks, refits = [], [], []

    for i in range(L, N - M + 1):
        if i > i_star:
            lw = cur[:, i_star:i].sum(axis=1)      # eq (10): product over the gap
            lw, k = psis_smooth(lw)
        else:
            lw, k = np.zeros(S), 0.0
        ks.append(k)

        if k > tau and refit is not None:          # weights degenerated -> refit
            cur = np.asarray(refit(i), float)
            i_star = i
            refits.append(i)
            lw = np.zeros(S)

        # eq (9): weighted average of the future predictive density
        tgt = cur[:, i:i + M].sum(axis=1)
        a = lw + tgt
        elpd_i.append(float(np.logaddexp.reduce(a) - np.logaddexp.reduce(lw)))

    return dict(elpd=float(np.sum(elpd_i)), pointwise=np.array(elpd_i),
                pareto_k=np.array(ks), refits=refits, n_refit=len(refits), L=L)
