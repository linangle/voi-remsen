"""State-conditional emission likelihoods for the two markers.

Both markers are conditionally independent given the latent severity, so the
per-state emission likelihood is the product of the two channels (a factor of 1
when a marker is missing that week).

Biomass (cells), tobit-Normal on y = log10(cells + 1), left-censored at 0:
    y > 0  -> Normal density  N(y | mu_s, sigma_s)
    y == 0 -> censored mass    Phi((0 - mu_s) / sigma_s)   (the non-detect spike)

Toxicity (DA), 3-category categorical, censoring absorbed as category 1:
    P(cat | s) = E_da[s, cat-1],  cat in {1 non-detect, 2 detected<20, 3 closure}

Each function returns an emission-likelihood matrix of shape (T, K): entry (t, s)
is P(observations at week t | state s).
"""

from __future__ import annotations

import numpy as np

SQRT2 = np.sqrt(2.0)
SQRT2PI = np.sqrt(2.0 * np.pi)

def _norm_pdf(z):
    return np.exp(-0.5 * z * z) / SQRT2PI


def _norm_cdf(z):
    from scipy.special import erfc
    return 0.5 * erfc(-z / SQRT2)


def emission_matrix_np(y_cell, da_cat, mu, sigma, E_da):
    """(T,K) emission likelihoods. y_cell / da_cat use NaN for missing weeks."""
    y_cell = np.asarray(y_cell, float)
    da_cat = np.asarray(da_cat, float)
    mu, sigma = np.asarray(mu, float), np.asarray(sigma, float)
    E_da = np.asarray(E_da, float)
    T, K = y_cell.shape[0], mu.shape[0]

    z = (y_cell[:, None] - mu[None, :]) / sigma[None, :]            # (T,K)
    dens = _norm_pdf(z) / sigma[None, :]
    cens = _norm_cdf((0.0 - mu[None, :]) / sigma[None, :])
    cell = np.where(y_cell[:, None] > 0, dens, cens)
    cell = np.where(np.isnan(y_cell)[:, None], 1.0, cell)

    da = np.ones((T, K))
    obs = ~np.isnan(da_cat)
    idx = (np.nan_to_num(da_cat, nan=1).astype(int) - 1)
    da_obs = E_da[:, idx].T                                          # (T,K)
    da = np.where(obs[:, None], da_obs, 1.0)
    return cell * da

def emission_matrix_pt(y_cell, cell_mask, da_onehot, da_mask, mu, sigma, E_da):
    """(T,K) emission likelihoods in pytensor.

    y_cell    (T,) log10 cells, NaN-safe filled (masked by cell_mask)
    cell_mask (T,) 1.0 where a cell count was observed
    da_onehot (T,3) one-hot of the DA category (zeros where missing)
    da_mask   (T,) 1.0 where DA was observed
    """
    import pytensor.tensor as pt

    z = (y_cell[:, None] - mu[None, :]) / sigma[None, :]            # (T,K)
    dens = pt.exp(-0.5 * z * z) / (sigma[None, :] * SQRT2PI)
    is_pos = (y_cell > 0)[:, None]
    zc = (0.0 - mu[None, :]) / sigma[None, :]
    cens = 0.5 * pt.erfc(-zc / SQRT2)
    cell = pt.switch(is_pos, dens, cens)
    cell = pt.switch(cell_mask[:, None] > 0, cell, 1.0)

    da_obs = da_onehot @ E_da.T                                     # (T,K)
    da = pt.switch(da_mask[:, None] > 0, da_obs, 1.0)
    return cell * da
