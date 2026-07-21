# voi-remsen

**Value of remote sensing for the HAB early-closure decision.**

This project quantifies what a *theoretical* remote-sensing (RS) tool would be
worth for the decision to close a shellfish harvest early under harmful-algal-
bloom (HAB) risk on the California coast. The guiding principle: **remote sensing
appears only in the observation channel `O`**; the world model `T` (bloom /
toxicity dynamics) is built entirely from *in-situ* truth — calHABMAP
*Pseudo-nitzschia* cell counts (biomass) and CDPH tissue domoic acid (toxicity).

The value question is then: **at what RS accuracy `p` does the tool add enough
value to be worth adopting**, on top of the current CDPH tissue-testing policy?

See [`docs/methodology.md`](docs/methodology.md) for the full math, modeling
decisions, and citations, and [`docs/HANDOFF.md`](docs/HANDOFF.md) for the phase
roadmap and next steps.

## Where the phases stand

- **Phase A — data assembly: done.** Weekly panels (cell / toxicity / joint) at
  the co-located anchor sites; coverage EDA with a CDPH-tissue overlay.
- **Phase B — continuous-time toxicity reference (R `msm`): done.** DA-only fit
  selects **K = 2** (clean / toxic-episode) by AIC & BIC. Validation oracle for
  the joint model.
- **Phase C — joint two-marker continuous-time Bayesian HMM (PyMC): next.**
  Shared latent severity driving both markers, `P(Δt)=expm(Δt·Q)`, censored
  emissions. Expected to resolve `K > 2` (cells add resolution DA lacks).
- **Then — decision model:** optimal-stopping POMDP on the joint state, RS as a
  swept-accuracy observation channel; `VoI(p)` + break-even. The reusable engine
  (`hab_pomdp/`) is already built and PBVI-based.

## Layout

```
src/voi_remsen/      world-model package (config, data, paths, viz)
hab_pomdp/           decision engine: PBVI solver, parametric sensor, HMM inference
scripts/             00_build_panels.py  (Phase A entrypoint)
R/                   fit_toxicity_msm.R  (Phase B continuous-time msm reference)
docs/                methodology.md, HANDOFF.md, accuracy_sweep.md
notebooks/           HAB_Sensor_Accuracy_Sweep.ipynb
notebooks/reference/ Arcieri's original HMM tutorials (reference only; see below)
data/input|interim|output/   raw markers -> weekly panels -> fitted artifacts
derived/             decision-engine artifacts (hmm_K*.npz, sweep_accuracy.npz)
```

## Reproduce

```bash
pip install -e .                         # puts voi_remsen on the path
python scripts/00_build_panels.py        # -> data/interim/*.csv  (expects ~558 joint weeks)
Rscript R/fit_toxicity_msm.R             # -> data/output/msm_tox_*  (K=2 for DA alone)
```

Phase-C inference additionally needs `pip install -e ".[worldmodel]"`
(pymc, nutpie, arviz).

## Data attribution

- CDPH tissue domoic acid: © California Department of Public Health,
  Environmental Management Branch.
- calHABMAP shore-station monitoring data.
- The `notebooks/reference/` tutorials are from Giacomo Arcieri's
  [Hidden-Markov-Models](https://github.com/giarcieri/Hidden-Markov-Models)
  repository, kept here as learning reference (they are not part of the pipeline).
  This project's POMDP inference follows the method of Arcieri et al. (2023, 2024).

The large C-HARM NetCDF cubes (`data/charmv3-*.nc`) are deprecated iteration-1
satellite data, kept locally but git-ignored (they are not used by the in-situ
world model).
