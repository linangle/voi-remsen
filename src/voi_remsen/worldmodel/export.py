"""Export the fitted world model for the decision layer.

Writes posterior draws of the generator Q and the weekly transition
P = expm(Q) (mean + samples, for Arcieri-style robust planning), plus the
emission parameters, to data/output/worldmodel.npz. This is the world-building
boundary: the decision model consumes only this file.
"""

from __future__ import annotations

import numpy as np

from ..paths import output_file
from .linalg import birthdeath_Q, expm_np


def _flat(idata, name):
    x = idata.posterior[name].values
    return x.reshape(-1, *x.shape[2:])


def export_worldmodel(idata, K, out="worldmodel.npz"):
    up = _flat(idata, "up").reshape(-1, K - 1)
    down = _flat(idata, "down").reshape(-1, K - 1)
    Qs = np.stack([birthdeath_Q(up[i], down[i]) for i in range(len(up))])
    P1s = np.stack([expm_np(Qs[i]) for i in range(len(Qs))])

    payload = dict(
        K=K,
        Q_mean=Qs.mean(0), Q_samples=Qs,
        P_week_mean=P1s.mean(0), P_week_samples=P1s,
        mu=_flat(idata, "mu"), sigma=_flat(idata, "sigma"),
        E_da=_flat(idata, "E_da"),
        rate_units="per_week",
    )
    path = output_file(out)
    np.savez(path, **payload)
    return path
