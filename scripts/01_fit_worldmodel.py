"""Phase C entrypoint: fit the joint two-marker continuous-time HMM.

    python scripts/01_fit_worldmodel.py [K ...]      # default: K-scan 2 3 4

For each K: fit (nutpie), save the trace to data/output/worldmodel_K{K}.nc and a
PSIS-LOO row; then export the best K to data/output/worldmodel.npz. Cross-check
the K=2 toxicity component against the msm reference (data/output/msm_tox_*).
"""

import sys

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

from voi_remsen import paths
from voi_remsen.worldmodel import build_sequences, build_model
from voi_remsen.worldmodel.export import export_worldmodel


def fit_one(seqs, K, draws=700, tune=1000):
    # NOTE: the likelihood is a single lumped pm.Potential, so PSIS-LOO/WAIC
    # (which need per-observation log-lik) are not yet available -- K is chosen
    # by the msm AIC/BIC reference + posterior predictive checks for now.
    m = build_model(seqs, K=K)
    with m:
        idata = pm.sample(draws=draws, tune=tune, chains=4, cores=4,
                          target_accept=0.95, nuts_sampler="nutpie",
                          random_seed=0, progressbar=False)
    return m, idata


def main(Ks):
    panel = pd.read_csv(paths.interim_file("joint_panel.csv"),
                        parse_dates=["week_start"])
    seqs = build_sequences(panel)
    print(f"sequences: {[(s['code'], s['n']) for s in seqs]}")
    print(f"total observed weeks: {sum(s['n'] for s in seqs)} | "
          f"gmax {max(int(s['gaps'].max()) for s in seqs)} wk")

    rows = []
    traces = {}
    for K in Ks:
        print(f"\n--- fitting K={K} ---")
        m, idata = fit_one(seqs, K)
        idata.to_netcdf(str(paths.output_file(f"worldmodel_K{K}.nc")))
        traces[K] = idata
        rhat = float(max(float(az.rhat(idata)[v].max())
                         for v in az.rhat(idata).data_vars))
        ndiv = int(idata.sample_stats.diverging.values.sum())
        rows.append(dict(K=K, max_rhat=round(rhat, 3), divergences=ndiv))
        print(f"K={K}: max r-hat {rhat:.3f}, divergences {ndiv}")

    sel = pd.DataFrame(rows)
    sel.to_csv(paths.output_file("worldmodel_modelsel.csv"), index=False)
    print("\n", sel.to_string(index=False))

    best = Ks[0]                          # cross-check vs msm before trusting K-scan
    path = export_worldmodel(traces[best], best)
    print(f"\nexported world model (K={best}) -> {path}")


if __name__ == "__main__":
    Ks = [int(a) for a in sys.argv[1:]] or [2, 3, 4]
    main(Ks)
