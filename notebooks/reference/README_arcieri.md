# Hidden-Markov-Models

This repository contains tutorials for the inference of transition matrices and observation generating process of Markov chains, Markov Decision Processes (MDPs), and Partially Observable Markov Decision Processes (POMDPs), via Markov Chain Monte Carlo (MCMC) inference of Hidden Markov models (HMM). While other examples for the inference of transition and observation models for Markov chains are available, the tutorials here included extend this inference to the case of Markov chains conditioned on actions, namely MDPs and POMDPs.

The code used in the notebook `HMMs for deterioration process Truncated Normal Process.ipynb` is close to the code used for the POMDP inference in the paper   https://arxiv.org/abs/2212.07933. The differences are represented by the use of real world data, instead of the simulated data of the notebooks, and the slightly different HMM used, fully described in the paper (e.g., Student's t steps instead of the more general Gaussian steps of the notebook).

It is suggested to go over the notebooks in the order: `HMMs for Markov Chains.ipynb`, `HMMs for MDP and POMDP.ipynb`, and `HMMs for deterioration process Truncated Normal Process.ipynb`.

---

## Application: Value of Information of C-HARM remote sensing for California HABs

This fork extends the POMDP machinery above to a real problem: the **early-harvest
decision under Harmful Algal Bloom (HAB) risk on the California coast**, using it
to quantify the **Value of Information (VoI)** of the **C-HARM** satellite product.

**Iteration 2 (current) — "at what accuracy is a remote sensor worth it?"**
An optimal-stopping POMDP (do-nothing / pay-to-look / close-terminal) on a bloom
world model **inferred from calHABMAP by MCMC** (Arcieri's method, PyMC+nutpie),
solved with **point-based value iteration**, sweeping a *theoretical* sensor
across all accuracies to find the break-even. Result: a sensor is worthless
below ~0.6 accuracy; break-even ~0.63–0.73.
- **`README_ACCURACY_SWEEP.md`** — full write-up.
- **`HAB_Sensor_Accuracy_Sweep.ipynb`** — narrative + plots.
- **`hab_pomdp/`** modules `inference`, `stopping`, `sensor`, `run_stopping`.
  Run: `python -m hab_pomdp.inference 3 400 600` then `python -m hab_pomdp.run_stopping`.

**Iteration 1 — VoI of the actual C-HARM product (its measured skill as the observation).**
- **`README_HAB_POMDP.md`** — mapping (states = binned *Pseudo-nitzschia*; C-HARM
  as observation), how the Arcieri / Bashar / Song papers are used, results, caveats.
- **`HAB_POMDP_CHARM_VoI.ipynb`** — walkthrough on real calHABMAP + C-HARM data.
- **`hab_pomdp/`** modules `config`, `data`, `pomdp`, `solve`, `run`
  (`python -m hab_pomdp.run`).

- **`data/`** — calHABMAP shore-station monitoring + C-HARM v3.1 NetCDF fields.
- **`derived/`** — cached match-ups, the MCMC world model, and sweep outputs.
