"""State-conditional emission likelihoods for the two markers.

Both markers are conditionally independent given the latent severity, so the
per-state emission likelihood is the product of the two channels (a factor of 1
when a marker is missing that week).

Biomass (cells), HURDLE on y = log10(cells + 1):
    y == 0 -> a per-state zero probability   pi0_s          (the non-detect spike)
    y  > 0 -> (1 - pi0_s) * TruncNormal(y | mu_s, sigma)    truncated at 0
The positive part uses ONE SHARED sigma across states, so the cell likelihood
ratio between an ordered pair of states is monotone in y: a larger count can
never make a lower state more likely. (An earlier per-state-sigma tobit had an
inflated low-state sigma that reversed the ordering above ~10^6 cells/L.)

Toxicity (DA), 3-category categorical, censoring absorbed as category 1:
    P(cat | s) = E_da[s, cat-1],  cat in {1 non-detect, 2 detected<20, 3 closure}

Each function returns an emission-likelihood matrix of shape (T, K): entry (t, s)
is P(observations at week t | state s).
"""

from __future__ import annotations

import numpy as np

SQRT2 = np.sqrt(2.0)
SQRT2PI = np.sqrt(2.0 * np.pi)


def _Phi_np(z):
    from scipy.special import erfc
    return 0.5 * erfc(-z / SQRT2)


def emission_matrix_np(y_cell, da_cat, mu, sigma, pi_zero, E_da):
    """(T,K) emission likelihoods. sigma is a shared scalar; pi_zero is (K,).

    y_cell / da_cat use NaN for missing weeks; an exact 0 in y_cell is a cell
    non-detect.
    """
    y_cell = np.asarray(y_cell, float)
    da_cat = np.asarray(da_cat, float)
    mu, pi_zero = np.asarray(mu, float), np.asarray(pi_zero, float)
    E_da = np.asarray(E_da, float)
    sigma = float(sigma)
    T, K = y_cell.shape[0], mu.shape[0]

    z = (y_cell[:, None] - mu[None, :]) / sigma                    # (T,K)
    dens = np.exp(-0.5 * z * z) / (sigma * SQRT2PI)
    Zk = _Phi_np(mu / sigma)                                       # (K,) trunc const
    pos = (1.0 - pi_zero)[None, :] * dens / Zk[None, :]            # y>0 density
    cell = np.where(y_cell[:, None] <= 0, pi_zero[None, :], pos)   # zero vs positive
    cell = np.where(np.isnan(y_cell)[:, None], 1.0, cell)          # missing -> 1

    da = np.ones((T, K))
    obs = ~np.isnan(da_cat)
    idx = (np.nan_to_num(da_cat, nan=1).astype(int) - 1)
    da_obs = E_da[:, idx].T                                         # (T,K)
    da = np.where(obs[:, None], da_obs, 1.0)
    return cell * da


def emission_matrix_pt(y_cell, cell_mask, da_onehot, da_mask, mu, sigma,
                       pi_zero, E_da):
    """(T,K) emission likelihoods in pytensor (hurdle cells + categorical DA).

    y_cell    (T,) log10 cells, NaN-safe filled (masked by cell_mask); 0 = non-detect
    cell_mask (T,) 1.0 where a cell count was observed
    da_onehot (T,3) one-hot of the DA category (zeros where missing)
    da_mask   (T,) 1.0 where DA was observed
    mu        (K,) ordered nonneg cell means; sigma scalar (shared); pi_zero (K,)
    """
    import pytensor.tensor as pt

    z = (y_cell[:, None] - mu[None, :]) / sigma                    # (T,K)
    dens = pt.exp(-0.5 * z * z) / (sigma * SQRT2PI)
    Zk = 0.5 * pt.erfc(-(mu / sigma) / SQRT2)                      # Phi(mu/sigma), (K,)
    pos = (1.0 - pi_zero)[None, :] * dens / Zk[None, :]
    cell = pt.switch((y_cell <= 0)[:, None], pi_zero[None, :], pos)
    cell = pt.switch(cell_mask[:, None] > 0, cell, 1.0)

    da_obs = da_onehot @ E_da.T                                    # (T,K)
    da = pt.switch(da_mask[:, None] > 0, da_obs, 1.0)
    return cell * da
