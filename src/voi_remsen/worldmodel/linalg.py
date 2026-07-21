"""Matrix exponential for the continuous-time transition model.

expm(dt * Q) is the transition probability over an elapsed time dt for a CTMC
with generator Q (Jackson eq. 2). We use scaling-and-squaring with a Taylor
core: expm(A) = (expm(A / 2^s))^(2^s). Both the divide and the s repeated
squarings are exact and differentiable, so this backprops cleanly under NUTS.

`expm_pt` is the pytensor version used inside the PyMC model; `expm_np` is a
numpy mirror used for simulation and testing.
"""

from __future__ import annotations

import numpy as np


def expm_pt(A, order: int = 12, s: int = 10):
    """expm(A) for a small square pytensor matrix (scaling-and-squaring Taylor)."""
    import pytensor.tensor as pt

    n = A.shape[0]
    As = A / (2.0 ** s)
    I = pt.eye(n)
    term = I
    E = I
    for k in range(1, order + 1):
        term = term @ As / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def expm_np(A: np.ndarray, order: int = 12, s: int = 10) -> np.ndarray:
    """Numpy mirror of `expm_pt` (same algorithm), for simulation/tests."""
    A = np.asarray(A, float)
    n = A.shape[0]
    As = A / (2.0 ** s)
    I = np.eye(n)
    term = I.copy()
    E = I.copy()
    for k in range(1, order + 1):
        term = term @ As / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def birthdeath_Q(up: np.ndarray, down: np.ndarray) -> np.ndarray:
    """Assemble a K x K birth-death generator from adjacent-level rates.

    up[i]   = rate i -> i+1   (i = 0..K-2)
    down[i] = rate i -> i-1   (i = 1..K-1, down[0] unused)
    Diagonal is set so rows sum to zero.
    """
    up = np.asarray(up, float)
    down = np.asarray(down, float)
    K = up.shape[0] + 1
    Q = np.zeros((K, K))
    for i in range(K - 1):
        Q[i, i + 1] = up[i]
    for i in range(1, K):
        Q[i, i - 1] = down[i - 1]
    for i in range(K):
        Q[i, i] = -Q[i].sum()
    return Q
