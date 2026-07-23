"""Phase D entrypoint: fit the non-stationary seasonal Q(day-of-year) model.

    python scripts/02_fit_seasonal.py [H]        # H = annual harmonics, default 1

Fits the seasonal model at the production K (read from worldmodel.npz), cross-
checks its per-observation log-likelihood, scores it against the stationary fit at
the same K with PSIS-LOO (does seasonality help?), and exports the 53 weekly
transition matrices P_woy to seasonal_worldmodel.npz.
"""

import sys

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

from voi_remsen import config, paths
from voi_remsen.worldmodel import build_sequences
from voi_remsen.worldmodel.seasonal import (
    build_grid_sequences, build_seasonal_model, season_design,
    periodic_occupancy, WOY)
from voi_remsen.worldmodel.linalg import expm_np_batched
from voi_remsen.worldmodel.selection import (
    pointwise_loglik, seasonal_pointwise_loglik)



def _max_rhat(idata):
    rh = az.rhat(idata)
    return float(max(float(rh[v].max()) for v in rh.data_vars))


def main(H, K=None):
    if K is None:                      # follow the production world model
        K = int(np.load(paths.output_file("worldmodel.npz"))["K"])
    print(f"seasonal fit at production K={K}")
    panel = pd.read_csv(paths.interim_file("joint_panel.csv"),
                        parse_dates=["week_start"])
    gseqs = build_grid_sequences(panel, sites=config.MODEL_SITES)
    print(f"sites: {[(g['code'], g['n'], g['n_grid']) for g in gseqs]}")
    print(f"observed weeks {sum(g['n'] for g in gseqs)} | "
          f"grid weeks {sum(g['n_grid'] for g in gseqs)} | harmonics H={H}")

    m = build_seasonal_model(gseqs, K=K, H=H)
    with m:
        idata = pm.sample(draws=500, tune=1000, chains=4, cores=4,
                          target_accept=0.95, nuts_sampler="nutpie",
                          random_seed=0, progressbar=False)
    seasonal_pointwise_loglik(idata, gseqs, K, H)          # attach BEFORE saving
    idata.to_netcdf(str(paths.output_file(f"seasonal_K{K}_H{H}.nc")))
    rhat, ndiv = _max_rhat(idata), int(idata.sample_stats.diverging.values.sum())
    print(f"\nmax r-hat {rhat:.3f} | divergences {ndiv}")
    if rhat >= 1.05 or ndiv > 0:
        raise SystemExit(f"seasonal fit not converged (r-hat {rhat:.3f}, "
                         f"div {ndiv}); trace saved, but not exporting.")

    lo_s = az.loo(idata, var_name="obs")
    print(f"seasonal LOO: elpd {lo_s.elpd:.1f} +/- {lo_s.se:.1f}, "
          f"p_loo {lo_s.p:.1f}, max k {float(lo_s.pareto_k.max()):.2f}")

    # stationary K=2 for the comparison
    stat = az.from_netcdf(paths.output_file(f"worldmodel_K{K}.nc"))
    pointwise_loglik(stat, build_sequences(panel, sites=config.MODEL_SITES), K)
    lo0 = az.loo(stat, var_name="obs")
    print(f"stationary LOO: elpd {lo0.elpd:.1f} +/- {lo0.se:.1f}, "
          f"p_loo {lo0.p:.1f}")
    cmp = az.compare({"seasonal": idata, "stationary": stat}, var_name="obs")
    cmp.to_csv(paths.output_file("seasonal_compare.csv"))
    print("\n=== az.compare (PSIS-LOO) ===")
    print(cmp.to_string())

    # export the seasonal transition cycle + occupancy curves (per draw)
    ncoef = season_design(np.arange(WOY), H).shape[1]
    cu = idata.posterior["coef_up"].values.reshape(-1, K - 1, ncoef)   # (S,K-1,ncoef)
    cd = idata.posterior["coef_down"].values.reshape(-1, K - 1, ncoef)
    Xd = season_design(np.arange(WOY), H)
    up_w = np.exp(np.einsum("wc,skc->swk", Xd, cu))                    # (S,53,K-1)
    down_w = np.exp(np.einsum("wc,skc->swk", Xd, cd))
    S = up_w.shape[0]
    # per-draw weekly generators, then the TRUE propagated occupancy of the top state
    Q = np.zeros((S, WOY, K, K))
    for i in range(K - 1):
        Q[:, :, i, i + 1] = up_w[:, :, i]
        Q[:, :, i + 1, i] = down_w[:, :, i]
    for k in range(K):
        Q[:, :, k, k] = -Q[:, :, k, :].sum(-1)
    P_woy = expm_np_batched(Q)                                          # (S,53,K,K)
    # frozen-rate equilibrium of a birth-death chain: pi_k ~ prod_{i<k} up_i/down_i
    ratio = np.cumprod(up_w / down_w, axis=2)                           # (S,53,K-1)
    w = np.concatenate([np.ones((S, WOY, 1)), ratio], axis=2)
    occ_frozen = (w / w.sum(2, keepdims=True))[:, :, K - 1]             # (S,53) diagnostic
    occ_prop = periodic_occupancy(P_woy)[:, :, K - 1]                   # (S,53) TRUE prevalence
    np.savez(paths.output_file("seasonal_worldmodel.npz"),
             K=K, H=H, coef_up=cu, coef_down=cd,
             P_woy_mean=P_woy.mean(0),                                 # posterior-mean matrices
             P_woy_samples=P_woy,
             up_woy_all=up_w.mean(0), down_woy_all=down_w.mean(0),
             mu=idata.posterior["mu"].values.reshape(-1, K),
             sigma=idata.posterior["sigma"].values.reshape(-1),        # shared scalar
             pi_zero=idata.posterior["pi_zero"].values.reshape(-1, K),
             pi0=idata.posterior["pi0"].values.reshape(-1, K),
             E_da=idata.posterior["E_da"].values.reshape(-1, K, 3),
             up_woy=up_w[:, :, -1].mean(0), down_woy=down_w[:, :, -1].mean(0),
             occ_frozen=occ_frozen.mean(0), occ_propagated=occ_prop.mean(0),
             rate_units="per_week", cell_transform="log10(cells+1)")
    upm, opm = up_w[:, :, -1].mean(0), occ_prop.mean(0)
    print(f"\ntop-state onset rate/wk range {upm.min():.3f}-{upm.max():.3f} "
          f"(peak woy {upm.argmax()}); propagated top-state occupancy "
          f"{opm.min():.2f}-{opm.max():.2f} (peak woy {opm.argmax()})")
    print(f"exported -> {paths.output_file('seasonal_worldmodel.npz')}")


if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    main(H)
