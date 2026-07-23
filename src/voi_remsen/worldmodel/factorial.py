"""Factorial (biomass x toxin-quota) world model.

Particulate DA is NOT a pure toxicity readout: pDA ~ cell abundance x toxin per
cell. Treating it as an observation of a toxicity axis would smuggle biomass
variance into that axis. So the latent state factorises into two chains

    B_t : Pseudo-nitzschia biomass level      (log10 cells)
    Q_t : toxin quota / toxigenic propensity  (log10 toxin per cell)

and, because the relationship is multiplicative, in logs it is additive:

    cells  observes   mu_B[b]
    pDA    observes   mu_B[b] + nu_Q[q]        <- this offset identifies Q
    tissue observes   a censored, monotone function of exposure e = mu_B + nu_Q

Empirically at SCW the quota is emphatically not constant (log10 quota spans
~4 orders, var 0.45 vs 0.88 for log cells; the regression of log pDA on log
cells has slope 0.73, not 1), which is what makes Q identifiable.

B and Q evolve as independent birth-death chains, so the joint generator is a
Kronecker sum and the joint transition factorises: P_joint(g) = P_B^g (x) P_Q^g.
That lets the existing scaled forward filter run unchanged on the K_B*K_Q joint
state.

DEFERRED: shellfish tissue burden H_t has uptake/depuration memory and is not
Markov in (B,Q). Here tissue responds to CURRENT exposure through an ordered
logit; adding H_t as a third, accumulating chain is the next step.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .linalg import expm_pt

SQRT2, SQRT2PI = np.sqrt(2.0), np.sqrt(2.0 * np.pi)


def build_bq_sequences(panel, sites=None):
    """Per-site sequences carrying cells, pDA and tissue with masks."""
    if sites is None:
        sites = sorted(panel["location_code"].dropna().unique())
    seqs = []
    for code in sites:
        d = panel[panel["location_code"] == code].copy()
        d = d.dropna(subset=["cells", "pda", "da"], how="all").sort_values("week_start")
        if len(d) < 2:
            continue
        wk = pd.to_datetime(d["week_start"])
        widx = ((wk - wk.min()).dt.days // 7).to_numpy()
        gaps = np.diff(widx, prepend=widx[0]).astype(int)

        cell_mask = d["cells"].notna().to_numpy().astype(float)
        y_cell = np.nan_to_num(d["cells_log10"].to_numpy(), nan=0.0) * cell_mask
        pda_mask = d["pda"].notna().to_numpy().astype(float)
        y_pda = np.nan_to_num(d["pda_log10"].to_numpy(), nan=0.0) * pda_mask
        pda_zero = ((d["pda"].to_numpy() == 0) & (pda_mask > 0)).astype(float)

        da_mask = d["da"].notna().to_numpy().astype(float)
        cat = np.nan_to_num(d["da_cat"].to_numpy(), nan=1).astype(int)
        onehot = np.zeros((len(d), 3))
        onehot[np.arange(len(d)), cat - 1] = da_mask
        # continuous tissue channel: actual DA value, censoring flag, per-sample LOD
        da_val = np.nan_to_num(d["da"].to_numpy(float), nan=1.0)
        da_cens = np.nan_to_num(d["cens"].to_numpy(float), nan=0.0) * da_mask
        da_lod = np.nan_to_num(d["lod"].to_numpy(float), nan=2.5)
        da_lod = np.where(da_lod > 0, da_lod, 2.5)

        seqs.append(dict(code=code, n=len(d), gaps=gaps, week_start=wk.to_numpy(),
                         y_cell=y_cell, cell_mask=cell_mask,
                         y_pda=y_pda, pda_mask=pda_mask, pda_zero=pda_zero,
                         da_onehot=onehot, da_mask=da_mask,
                         da_val=da_val, da_cens=da_cens, da_lod=da_lod))
    return seqs


def _kron_pt(A, B):
    """Kronecker product of two square pytensor matrices."""
    import pytensor.tensor as pt
    m, n = A.shape[0], B.shape[0]
    return (A[:, None, :, None] * B[None, :, None, :]).reshape((m * n, m * n))


def _bd_Q(up, down, K):
    import pytensor.tensor as pt
    Q = pt.zeros((K, K))
    for i in range(K - 1):
        Q = pt.set_subtensor(Q[i, i + 1], up[i])
        Q = pt.set_subtensor(Q[i + 1, i], down[i])
    return Q - pt.diag(Q.sum(axis=1))


def build_bq_model(seqs, K_B=3, K_Q=2, sigma_floor=0.15, mu_sep=0.5, nu_sep=0.3,
                   quota_centre=-4.85, gap_cap=52, tissue_mode="lognormal"):
    """Factorial (B,Q) PyMC model. Joint state index j = b*K_Q + q.

    tissue_mode:
      "lognormal" -- tissue burden H = H_bg + c * 10^(mu_B+nu_Q), observed as a
        LEFT-CENSORED lognormal at each sample's own detection limit. Keeps the DA
        magnitude and the real LODs instead of collapsing to 3 categories, and
        P(closure) is then DERIVED as P(H >= 20 ppm) rather than fitted directly.
        c absorbs alpha/k: with weekly retention ~4% (fitted k ~ 0.45/day, see
        tissue.py) the kinetic stock is at quasi-equilibrium between weekly
        samples, so H is a function of current exposure at this resolution.
      "ordered" -- the earlier 3-category ordered logit (kept for comparison).
    """
    import pymc as pm
    import pytensor.tensor as pt
    from .model import _forward_loglik

    KJ = K_B * K_Q
    gmax = min(int(max(s["gaps"].max() for s in seqs)), gap_cap)
    mu_prior = np.linspace(0.5, 4.5, K_B)

    with pm.Model() as m:
        # --- dynamics: two independent birth-death chains -------------------
        upB = pm.Exponential("up_B", 1.0, shape=K_B - 1)
        downB = pm.Exponential("down_B", 1.0, shape=K_B - 1)
        upQ = pm.Exponential("up_Q", 1.0, shape=K_Q - 1)
        downQ = pm.Exponential("down_Q", 1.0, shape=K_Q - 1)
        PB = expm_pt(_bd_Q(upB, downB, K_B))
        PQ = expm_pt(_bd_Q(upQ, downQ, K_Q))
        # P_joint^g = P_B^g (x) P_Q^g  -- build the power stack directly
        powsB, powsQ = [pt.eye(K_B)], [pt.eye(K_Q)]
        for _ in range(gmax):
            powsB.append(powsB[-1] @ PB)
            powsQ.append(powsQ[-1] @ PQ)
        pows = pt.stack([_kron_pt(powsB[g], powsQ[g]) for g in range(gmax + 1)])

        # --- biomass levels (ordered, non-negative) -------------------------
        mu0 = pm.TruncatedNormal("mu0", mu=float(mu_prior[0]), sigma=1.0, lower=0.0,
                                 initval=float(max(mu_prior[0], 0.2)))
        dmu = pm.HalfNormal("dmu", 1.5, shape=K_B - 1,
                            initval=np.maximum(np.diff(mu_prior) - mu_sep, 0.3))
        mu_B = pm.Deterministic("mu_B", mu0 + pt.concatenate(
            [pt.zeros(1), pt.cumsum(mu_sep + dmu)]))

        # --- toxin quota levels (ordered; sets the pDA offset) --------------
        nu0 = pm.Normal("nu0", mu=quota_centre, sigma=1.0, initval=quota_centre)
        dnu = pm.HalfNormal("dnu", 1.0, shape=K_Q - 1, initval=np.full(K_Q - 1, 0.6))
        nu_Q = pm.Deterministic("nu_Q", nu0 + pt.concatenate(
            [pt.zeros(1), pt.cumsum(nu_sep + dnu)]))

        sig_c = pm.Deterministic("sigma_cell", sigma_floor + pm.HalfNormal("sc_raw", 1.0))
        sig_p = pm.Deterministic("sigma_pda", sigma_floor + pm.HalfNormal("sp_raw", 1.0))
        pz_B = pm.Beta("pi_zero_cell", 1.0, 1.0, shape=K_B)   # cell non-detect | B

        # exposure e = log10 predicted pDA, per joint state
        e_j = (mu_B[:, None] + nu_Q[None, :]).reshape((KJ,))
        mu_c_j = pt.repeat(mu_B, K_Q)                          # cells depend on B only
        pz_c_j = pt.repeat(pz_B, K_Q)

        # pDA non-detect probability falls with exposure (logistic link)
        a0 = pm.Normal("pda_zero_a0", 0.0, 2.0)
        a1 = pm.HalfNormal("pda_zero_a1", 2.0)
        pz_p_j = pm.math.sigmoid(a0 - a1 * (e_j - quota_centre))

        if tissue_mode == "lognormal":
            # tissue burden as a stock at quasi-equilibrium with current exposure
            log_c = pm.Normal("log_c", 0.0, 3.0)        # absorbs alpha/k and units
            H_bg = pm.HalfNormal("H_bg", 1.0)
            sig_t = pm.HalfNormal("sigma_tissue", 2.0)
            H_j = H_bg + pt.exp(log_c) * pt.power(10.0, e_j)           # (KJ,)
            logH = pt.log(pt.maximum(H_j, 1e-4))
            pm.Deterministic("H_state", H_j)
            pm.Deterministic("p_closure",                              # derived
                             0.5 * pt.erfc((pt.log(20.0) - logH) / (sig_t * SQRT2)))
        else:
            beta = pm.HalfNormal("tissue_beta", 2.0)
            c1 = pm.Normal("tissue_c1", 0.0, 2.0)
            dc = pm.HalfNormal("tissue_dc", 2.0)
            s_det = pm.math.sigmoid(beta * (e_j - quota_centre - c1))
            s_clo = pm.math.sigmoid(beta * (e_j - quota_centre - c1 - dc))
            E_da = pt.stack([1.0 - s_det, s_det - s_clo, s_clo], axis=1)   # (KJ,3)

        pi0 = pm.Dirichlet("pi0", a=np.ones(KJ))

        total = pt.as_tensor_variable(0.0)
        for s in seqs:
            yc = pt.as_tensor_variable(s["y_cell"])
            cm = pt.as_tensor_variable(s["cell_mask"])
            yp = pt.as_tensor_variable(s["y_pda"])
            pm_ = pt.as_tensor_variable(s["pda_mask"])
            pzero = pt.as_tensor_variable(s["pda_zero"])
            oh = pt.as_tensor_variable(s["da_onehot"])
            dm = pt.as_tensor_variable(s["da_mask"])

            # cells: hurdle + shared-sigma truncated normal (monotone in count)
            zc = (yc[:, None] - mu_c_j[None, :]) / sig_c
            dc_ = pt.exp(-0.5 * zc * zc) / (sig_c * SQRT2PI)
            Zc = 0.5 * pt.erfc(-(mu_c_j / sig_c) / SQRT2)
            cell = pt.switch((yc <= 0)[:, None], pz_c_j[None, :],
                             (1.0 - pz_c_j)[None, :] * dc_ / Zc[None, :])
            cell = pt.switch(cm[:, None] > 0, cell, 1.0)

            # pDA: hurdle + normal centred on mu_B + nu_Q (the B*Q product)
            zp = (yp[:, None] - e_j[None, :]) / sig_p
            dp_ = pt.exp(-0.5 * zp * zp) / (sig_p * SQRT2PI)
            pda = pt.switch(pzero[:, None] > 0, pz_p_j[None, :],
                            (1.0 - pz_p_j)[None, :] * dp_)
            pda = pt.switch(pm_[:, None] > 0, pda, 1.0)

            if tissue_mode == "lognormal":
                yv = pt.as_tensor_variable(s["da_val"])
                cflag = pt.as_tensor_variable(s["da_cens"])
                lodv = s["da_lod"]
                zo = (pt.log(pt.maximum(yv, 1e-6))[:, None] - logH[None, :]) / sig_t
                dens = pt.exp(-0.5 * zo ** 2) / (sig_t * SQRT2PI * pt.maximum(yv, 1e-6)[:, None])
                zc = (np.log(np.maximum(lodv, 1e-6))[:, None] - logH[None, :]) / sig_t
                cenp = 0.5 * pt.erfc(-zc / SQRT2)
                da = pt.switch(cflag[:, None] > 0, cenp, dens)
                da = pt.switch(dm[:, None] > 0, da, 1.0)
            else:
                da = pt.switch(dm[:, None] > 0, oh @ E_da.T, 1.0)

            gaps = np.minimum(s["gaps"], gmax).astype("int64")
            total = total + _forward_loglik(
                cell * pda * da, pt.as_tensor_variable(gaps), pi0, pows)
        pm.Potential("loglik", total)
    return m
