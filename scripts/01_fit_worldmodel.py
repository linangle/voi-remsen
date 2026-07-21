"""Phase C entrypoint: fit the joint two-marker continuous-time HMM + K-scan.

    python scripts/01_fit_worldmodel.py [K ...]      # default: K-scan 2 3 4

For each K: fit (nutpie), reconstruct per-observation log-likelihoods from the
forward filter (selection.pointwise_loglik), and score with PSIS-LOO. Saves the
trace to worldmodel_K{K}.nc, a per-K diagnostics row (worldmodel_kscan.csv), the
az.compare ranking (worldmodel_compare.csv), and exports the LOO-best K to
worldmodel.npz. (WAIC is not in arviz 1.2; PSIS-LOO supersedes it.)
"""

import sys

import arviz as az
import pandas as pd
import pymc as pm

from voi_remsen import paths
from voi_remsen.worldmodel import build_sequences, build_model
from voi_remsen.worldmodel.export import export_worldmodel
from voi_remsen.worldmodel.selection import pointwise_loglik


def fit_one(seqs, K, draws=500, tune=800):
    m = build_model(seqs, K=K)
    with m:
        idata = pm.sample(draws=draws, tune=tune, chains=4, cores=4,
                          target_accept=0.95, nuts_sampler="nutpie",
                          random_seed=0, progressbar=False)
    pointwise_loglik(idata, seqs, K)          # attach log_likelihood['obs']
    return idata


def _max_rhat(idata):
    rh = az.rhat(idata)
    return float(max(float(rh[v].max()) for v in rh.data_vars))


def main(Ks):
    panel = pd.read_csv(paths.interim_file("joint_panel.csv"),
                        parse_dates=["week_start"])
    seqs = build_sequences(panel)
    print(f"sequences: {[(s['code'], s['n']) for s in seqs]}")
    print(f"total observed weeks: {sum(s['n'] for s in seqs)} | "
          f"gmax {max(int(s['gaps'].max()) for s in seqs)} wk")

    rows, traces = [], {}
    for K in Ks:
        print(f"\n--- fitting K={K} ---")
        idata = fit_one(seqs, K)
        idata.to_netcdf(str(paths.output_file(f"worldmodel_K{K}.nc")))
        traces[f"K{K}"] = idata
        lo = az.loo(idata, var_name="obs")
        ndiv = int(idata.sample_stats.diverging.values.sum())
        rows.append(dict(K=K, elpd_loo=round(float(lo.elpd), 1),
                         se=round(float(lo.se), 1), p_loo=round(float(lo.p), 1),
                         max_pareto_k=round(float(lo.pareto_k.max()), 2),
                         max_rhat=round(_max_rhat(idata), 3), divergences=ndiv))
        print(f"K={K}: elpd_loo {lo.elpd:.1f} +/- {lo.se:.1f}, p_loo {lo.p:.1f}, "
              f"max k {float(lo.pareto_k.max()):.2f}, "
              f"r-hat {rows[-1]['max_rhat']}, div {ndiv}")

    sel = pd.DataFrame(rows)
    sel.to_csv(paths.output_file("worldmodel_kscan.csv"), index=False)
    print("\n=== K-scan diagnostics ===")
    print(sel.to_string(index=False))

    cmp = az.compare(traces, var_name="obs")     # LOO ranking (best first)
    cmp.to_csv(paths.output_file("worldmodel_compare.csv"))
    print("\n=== az.compare (PSIS-LOO) ===")
    print(cmp.to_string())

    # Only select among CONVERGED fits: LOO from unmixed chains is meaningless,
    # and a p_loo far above the parameter count flags an unreliable estimate.
    RHAT_MAX, PLOO_MAX = 1.05, 3 * max(Ks) ** 2
    ok = sel[(sel.max_rhat < RHAT_MAX) & (sel.p_loo < PLOO_MAX)]
    dropped = sel[~sel.K.isin(ok.K)]
    if len(dropped):
        print("\nEXCLUDED from selection (not converged / unreliable LOO):")
        print(dropped[["K", "max_rhat", "p_loo"]].to_string(index=False))
    if ok.empty:
        raise SystemExit("No K converged; not exporting a world model.")
    best = int(ok.sort_values("elpd_loo").iloc[-1].K)   # best LOO among converged
    path = export_worldmodel(traces[f"K{best}"], best)
    print(f"\nselected K={best} (LOO-best among converged) -> {path}")


if __name__ == "__main__":
    Ks = [int(a) for a in sys.argv[1:]] or [2, 3, 4]
    main(Ks)
