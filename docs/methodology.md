# voi-remsen: Methodology and Handoff

**Value of remote-sensing information for harmful-algal-bloom harvest decisions**

Status at handoff: the *world model* (the ecological process and its observation
channels) is built, validated, and corrected. The *decision layer* (actions,
rewards, policy, value-of-information) is **not yet implemented**. This document
records what exists, the mathematical and scientific basis for every modelling
choice, what was found to be wrong and fixed, what the numbers currently are, and
what remains.

---

## Contents

1. [The question and the design principle](#1-the-question-and-the-design-principle)
2. [Primer: the concepts, from the ground up](#2-primer-the-concepts-from-the-ground-up)
3. [Data: sources, matching, and what is derived](#3-data-sources-matching-and-what-is-derived)
4. [The world model and why it has this form](#4-the-world-model-and-why-it-has-this-form)
5. [Validation methodology](#5-validation-methodology)
6. [Current results](#6-current-results)
7. [Errors found and corrected](#7-errors-found-and-corrected)
8. [Known limitations](#8-known-limitations)
9. [What to do next](#9-what-to-do-next)
10. [Code map and reproduction](#10-code-map-and-reproduction)

---

## 1. The question and the design principle

### 1.1 The question

Shellfish harvests are closed when domoic acid (DA) in shellfish tissue reaches
**20 µg/g** ("20 ppm"), the regulatory action level. Confirming that requires
laboratory tissue assays, which are slow, sparse, and expensive. Satellites cannot
measure tissue toxicity, but they can (imperfectly) inform estimates of
phytoplankton biomass.

> **How much is a remote-sensing signal worth for an early-closure / harvest decision?**

Value of information is *not* accuracy, correlation, or mutual information. It is
the improvement in an **optimised decision objective**:

$$\mathrm{VoI} = \big[\text{optimal value with the information option}\big] - \big[\text{optimal value without it}\big].$$

If the information action is optional and the with-information action set contains
every baseline action, gross VoI cannot be negative under exact optimisation — the
decision maker can always ignore the new source. Any fixed adoption cost is
subtracted separately.

### 1.2 The design principle

> **Remote sensing appears only in the observation channel $O$. The world model
> $T$ (the bloom / toxicity dynamics) is built entirely from in-situ truth.**

This keeps the causal object (what the ocean does) separate from the informational
object (what a sensor tells you about it), so that sweeping sensor quality does not
silently change the ecology being modelled. Every choice below respects it:
in-situ cell counts, in-situ particulate DA, and laboratory tissue assays build $T$;
nothing derived from satellites enters $T$.

---

## 2. Primer: the concepts, from the ground up

### 2.1 The ecological variables are related but not interchangeable

Four distinct quantities appear in this problem:

| quantity | what it is | who measures it |
|---|---|---|
| **PN cell abundance** | *Pseudo-nitzschia* cells per litre in the water | calHABMAP, microscopy |
| **particulate DA (pDA)** | toxin bound in particulate matter, ng/mL of seawater | calHABMAP |
| **tissue DA** | toxin accumulated in shellfish flesh, µg/g | CDPH laboratory assay |
| **remote sensing** | ocean colour → chlorophyll, SST → *inferred* biomass | satellite |

Critically, **not every bloom is toxic**. Toxin production per cell varies with
species, strain, physiology, and nutrient stress. Treating all four as noisy readings
of one scalar "severity" is a substantive scientific hypothesis, not a neutral
simplification — and §4.5 shows the data rejects it.

Two relationships matter mechanistically:

$$\text{pDA} \;\approx\; \underbrace{B}_{\text{biomass}} \times \underbrace{Q}_{\text{toxin per cell}}, \qquad
\text{tissue DA} \;=\; \text{a stock accumulated from pDA exposure}.$$

### 2.2 Hidden Markov models and belief

A Markov chain says the next state depends only on the current state. For discrete
states $S_t$, the transition matrix is $T_{ij} = P(S_{t+1}=j \mid S_t = i)$.

A **hidden** Markov model adds that $S_t$ is not observed; instead we see $Y_t$ with
emission probability $P(Y_t \mid S_t)$. The analyst carries a **belief** $b_t$, a
probability distribution over states, updated by predict-then-correct:

$$\tilde b_{t+1}(j) = \sum_i b_t(i)\,T_{ij},
\qquad
b_{t+1}(j) = \frac{P(Y_{t+1}\mid S_{t+1}=j)\,\tilde b_{t+1}(j)}{\sum_k P(Y_{t+1}\mid S_{t+1}=k)\,\tilde b_{t+1}(k)}.$$

A belief is strictly more useful than a decoded label: "state 2" hides whether the
model is 51% or 99.9% sure, and those should often produce different decisions.

### 2.3 Continuous time and irregular observations

Field sampling is irregular — weekly at best, with multi-year gaps. Treating a
4-week gap as one transition would be wrong. Instead we model the latent process in
continuous time with a **generator** $Q$, whose off-diagonal entries $q_{ij}\ge 0$
are instantaneous rates and whose rows sum to zero. The transition over elapsed time
$\Delta$ is the matrix exponential

$$P(\Delta) = \exp(\Delta Q).$$

This is the standard construction for panel data (Jackson, *msm*). It correctly
allows that an observed $1 \to 3$ move over a long interval may have passed through
state 2 unobserved.

We restrict $Q$ to **birth–death** form (adjacent levels only):

$$Q_{i,i+1} = u_i, \qquad Q_{i+1,i} = d_i, \qquad Q_{ii} = -\!\!\sum_{j\neq i} Q_{ij},$$

which encodes that severity escalates and relaxes through intermediate levels rather
than jumping. A useful consequence: a birth–death chain satisfies detailed balance,
so its stationary distribution is available in closed form,

$$\pi_k \;\propto\; \prod_{i<k} \frac{u_i}{d_i}.$$

### 2.4 The forward algorithm

We never sample the discrete states. They are marginalised analytically by the
**scaled forward filter**. With $E_t(k) = P(Y_t \mid S_t = k)$ and gap $g_t$:

$$\alpha_1 \propto \pi_0 \odot E_1, \qquad
\alpha_t \propto \big(\alpha_{t-1} P(g_t)\big) \odot E_t,$$

rescaling at each step by $c_t = \sum_k \big(\alpha_{t-1}P(g_t)\big)_k E_t(k)$. Then

$$\log P(Y_{1:T}) = \sum_t \log c_t,
\qquad\text{and importantly}\qquad
\log c_t = \log P(Y_t \mid Y_{1:t-1}).$$

**Each scaling constant is a one-step-ahead predictive density.** This fact is used
twice below: to reconstruct per-observation log-likelihoods (§5.3), and to run
leave-future-out cross-validation almost for free (§2.7).

### 2.5 Missing markers and gaps

When a marker is not measured in a week, its emission factor is set to **1 for every
state**. Multiplying every state by 1 leaves the belief's shape untouched: a
measurement that was never taken says nothing. This is the correct likelihood, not a
patch — $P(\text{observed data}\mid\text{state})$ genuinely excludes it. The state
still propagates through $P(g)$, so a tissue-only week still informs the latent
process and still constrains the transition rates.

**Assumption this rests on:** data are missing at random given the model. If
sampling were triggered by visible blooms, absence would itself be informative and
this would bias the fit. The observed gaps here are multi-year, programme-level ones
(e.g. Morro Bay's median closure-to-cell gap is ~21 years), which look like funding
and monitoring eras rather than event-triggered sampling — but this has not been
formally audited (§8).

### 2.6 Censoring

Both markers are censored, and in different ways.

**Cell counts** have a large spike of exact zeros (66% of SCW cell weeks). Two
candidate treatments:
- *Tobit*: a latent normal censored at zero, $Y=\max(Y^\*,0)$, so
  $P(Y=0) = \Phi(-\mu/\sigma)$.
- *Hurdle*: a separate probability of a structural zero, and a distribution over
  positives.

We use the **hurdle** — §4.3 shows the tobit form caused a serious pathology.

**Tissue DA** is left-censored at a laboratory detection limit that *changed over
time* (1.0 → 2.5 µg/g around 2010). A non-detect is not a value; it is the statement
$Y < \mathrm{LOD}$, contributing

$$P(Y < \mathrm{LOD}) = \Phi\!\left(\frac{\log \mathrm{LOD} - \log H}{\sigma}\right)$$

with **each sample's own LOD**, not a pooled constant.

### 2.7 Model comparison: why LOO is not enough

Leave-one-out cross-validation (PSIS-LOO) estimates expected log predictive density
by removing one observation at a time. For **time series this is positively biased**:
the posterior still sees the future, so LOO systematically overstates forecast skill
(Bürkner, Gabry & Vehtari 2020).

The right target for a forecasting claim is **leave-future-out**:

$$\mathrm{ELPD}_{\mathrm{LFO}} = \sum_{i=L}^{N-M} \log p\big(y_{i+1:i+M} \mid y_{1:i}\big),$$

where *both* the predictive density and the parameter posterior condition only on the
past. Exact LFO needs a refit at every $i$. The PSIS approximation refits at $i^\*=L$
and then reuses that posterior with importance ratios

$$r_i^{(s)} \;\propto\; \prod_{j=i^\*+1}^{i} p\big(y_j \mid y_{1:j-1}, \theta^{(s)}\big),$$

Pareto-smoothed for stability, refitting only when the fitted Pareto shape $\hat k$
exceeds $\tau = 0.7$ (weights have degenerated). Because §2.4 gives us
$p(y_j\mid y_{<j},\theta) = \exp(\log c_j)$ directly from the forward filter, these
ratios are available at no extra modelling cost.

---

## 3. Data: sources, matching, and what is derived

### 3.1 Sources

| source | content | scale |
|---|---|---|
| **calHABMAP** (`calhabmap.csv`) | PN cell counts, particulate DA, SST, chlorophyll, nutrients | 7,203 samples, 17 shore stations, with lat/lon |
| **CDPH** (`CDPH_tissueDA.xlsx`) | shellfish tissue DA, censoring flag, detection limit, species, lat/lon, county | 12,682 records (12,407 geolocated), 447 distinct sample sites |

### 3.2 Tissue-to-station matching is spatial, not by name

**Original behaviour (defective).** Tissue was joined to stations by exact
`Sample Site` string equality against four hand-named sites. Any tissue record whose
label never equals a station name was silently dropped — for example sites recorded
under a county rather than a place name.

**Correction.** Each tissue sample is assigned to its **nearest calHABMAP station
within a radius**, by great-circle (haversine) distance:

$$d = 2R \arcsin\sqrt{\sin^2\!\tfrac{\Delta\phi}{2} + \cos\phi_1\cos\phi_2 \sin^2\!\tfrac{\Delta\lambda}{2}}, \quad R = 6371.0088\ \mathrm{km}.$$

**Effect** (all stations pooled):

| | exact name | **3 km** | 5 km |
|---|---|---|---|
| stations with tissue | 4 | **15** | 15 |
| tissue samples matched | 1,159 | **4,513** | 4,695 |
| tissue weeks | 1,110 | **3,412** | 3,522 |
| matchups (both markers, same week) | 558 | **758** | 773 |
| closure weeks (≥20 ppm) | 34 | **66** | 68 |

**Why 3 km is the default.** Going 3→5 km adds only 15 matchups and 2 closures, while
the median match distance is already 0.73 km (90th percentile 2.36 km) — nearly every
tissue site is effectively *at* a station. Meanwhile 5 km creates genuine ambiguity
where stations are closer together than the radius: **MBB↔MBF are 4.73 km apart** and
BBB↔BML 1.12 km, so a sample can fall within reach of two stations and be assigned by
nearest-wins. Set by `config.MATCH_RADIUS_KM`.

### 3.3 A hard limit discovered, and a negative result worth recording

Of the 66 closure weeks, **only 20 co-occur with a cell count**. We tested whether a
lag window (motivated by bioaccumulation) recovers more:

| window | tissue paired | **closures paired** |
|---|---|---|
| same week | 758 | **20** |
| ±1 week | 790 | **20** |
| cells lead 0–2 wk | 786 | **20** |
| ±2 weeks | 799 | **20** |

**A lag window recovers exactly zero additional closure events.** The reason is that
the gap distribution is bimodal: 20 closures have a same-week cell count, **none at
1–4 weeks**, and 44 are >26 weeks away or at stations with no cell data at all
(Humboldt has 15 closures and *zero* cell counts ever). This is a **monitoring-coverage
problem, not a timing-alignment problem**. Do not re-attempt lag windows to harvest
more data.

### 3.4 Site anchoring: Santa Cruz Wharf

The model is fit to **SCW alone** (`config.ANCHOR_SITE`). Rationale:

- **Pooling confounded state with site.** With four sites sharing one generator and
  one emission system, the latent state partly encoded *site identity*: at pooled
  $K=2$ the smoothed occupancy was 0.89 at Scripps versus 0.26 at Santa Cruz. The
  state was close to a proxy for which site you were looking at.
- **SCW carries the toxicity signal**: 22 of the 34 closure observations in the
  original panel, 663 tissue weeks, and near-**balanced marker coverage** (667 cell
  vs 663 DA weeks, ≈1:1, against ≈3:1 pooled). Balance matters: with cells
  outnumbering tissue 3:1, the shared latent state is dominated by biomass.
- **Species homogeneity**: 702 of 706 SCW tissue records are *Mytilus californianus*,
  so tissue kinetics are a single-species model rather than a mixture.

**Documented caveat:** SCW cell counts are the *seriata* group only (*delicatissima*
is not reported there). Internally consistent, but **not comparable across sites** —
this is why the totals must not be pooled without care.

### 3.5 Exact dates for kinetics

Weekly summaries are adequate for the latent-state models but **not** for tissue
kinetics, where within-week ordering matters. `tissue.build_exact_series` preserves
true sample dates and per-sample detection limits. Facts about the SCW design:

- water sampling is essentially weekly (median interval 7 days);
- **97% of tissue samples are collected the same day as a water sample** (max offset
  5 days), so the feared within-week misalignment is largely absent at SCW;
- the kinetic model is limited to the water era (2011-10 onward): **355 tissue
  observations, 11 closures**, ~80% censored.

---

## 4. The world model and why it has this form

The model evolved through several corrections. Each subsection states the choice, the
mathematics, and the evidence for it.

### 4.1 Skeleton

One ordered latent severity state per week, continuous-time birth–death transitions,
markers conditionally independent given the state:

$$P(Y^{\text{cell}}_t, Y^{\text{DA}}_t \mid S_t) = P(Y^{\text{cell}}_t\mid S_t)\,P(Y^{\text{DA}}_t\mid S_t).$$

Gaps use $P(g) = \exp(Q)^g$ on the weekly grid. Powers are capped at
`gap_cap = 52` weeks: the chain mixes on a scale of weeks, so $P^{52}$ is already
essentially stationary — verified as numerically negligible (max difference between
$P^{52}$ and $P^{255}$ ≈ $1.4\times10^{-4}$) while keeping the computational graph small.

### 4.2 Identifiability constraints

Label switching is the classic HMM pathology. We anchor states by **ordering the cell
means** and enforcing a **minimum separation**:

$$\mu_0 \ge 0, \qquad \mu_k = \mu_0 + \sum_{i<k}(\texttt{mu\_sep} + \delta_i), \quad \delta_i \ge 0 .$$

Both parts were necessary and both were discovered empirically:

- **Non-negativity.** $y = \log_{10}(\text{cells}+1) \ge 0$ by construction. Without the
  constraint the near-all-censored low state drove $\mu \to -\infty$, because under a
  tobit the zero-probability $\Phi\big((0-\mu)/\sigma\big) \to 1$ as $\mu\to-\infty$.
  That is a flat ridge, and chains parked at different points along it
  (max $\hat R = 1.55$ with **zero** divergences — the Stan-documented signature of
  multimodality rather than bad geometry).
- **Minimum separation.** At $K \ge 3$ a state would collapse onto its neighbour
  (identical cell mean, weakly told apart by DA). Diagnosed by per-chain posterior
  means: three chains found $\mu \approx [0, 3.70, 4.77]$ while one found
  $[0, 0.05, 4.07]$.

### 4.3 The cell emission: why hurdle with a *shared* variance

**The defect.** With a per-state-$\sigma$ tobit, the fitted low state had
$\sigma_0 = 3.35$ — inflated to absorb both the zero spike and site heterogeneity.
The consequence is not cosmetic. For two normals the log-likelihood ratio is

$$\log\frac{\mathcal N(y\mid\mu_1,\sigma_1)}{\mathcal N(y\mid\mu_0,\sigma_0)}
= \underbrace{-\tfrac12\left(\tfrac{1}{\sigma_1^2}-\tfrac{1}{\sigma_0^2}\right)y^2}_{\text{quadratic}} + \cdots$$

**If $\sigma_1 \neq \sigma_0$ the ratio is quadratic in $y$, hence non-monotone.**
Empirically the "bloom" state stopped being the more likely state above
$y \approx 5.93$ (~857,000 cells/L), and at the observed maximum
(3.2 M cells/L) the bloom-state cell likelihood was only **15%** of the low state's.
The low state also placed 3.7% of its mass above $10^6$ cells/L and 0.36% above
$10^9$ — biologically impossible. **A larger bloom could read as less severe**, which
is disqualifying for a hazard model.

**The fix.** A hurdle with one shared $\sigma$:

$$P(y\mid S=k) = \begin{cases}
\pi^{0}_k & y = 0\\[2pt]
(1-\pi^{0}_k)\,\dfrac{\mathcal N(y\mid \mu_k,\sigma)}{\Phi(\mu_k/\sigma)} & y > 0
\end{cases}$$

Separating the zero process stops $\sigma$ from having to absorb the spike, and the
**shared** $\sigma$ makes the quadratic term vanish, so the likelihood ratio is
$\exp\!\big((\mu_1-\mu_0)(y - \bar\mu)/\sigma^2\big)$ — **monotone increasing in $y$ by
construction**. Verified after refitting: $P(\text{top state}\mid \text{cells})$ rises
monotonically to 0.998 at the observed maximum.

This change had a large downstream consequence: it **removed the $K\ge3$
non-convergence entirely**. What had looked like fundamental multimodality was an
artefact of the broken emission.

### 4.4 Seasonality: $Q(\text{day-of-year})$

Rates become smooth functions of week-of-year through a log link on annual Fourier
terms:

$$\log u_i(w) = \beta^u_{i0} + \sum_{h=1}^{H}\Big[\beta^u_{ih,c}\cos\tfrac{2\pi h w}{52.18} + \beta^u_{ih,s}\sin\tfrac{2\pi h w}{52.18}\Big].$$

Because $Q$ now changes weekly, **a gap's transition is an ordered product, not a
power**:

$$P(t \to t+g) = \prod_{s=t}^{t+g-1} \exp\!\big(Q(w_s)\big).$$

The forward filter therefore runs over the **full weekly grid** (observed and
unobserved weeks), applying the appropriate weekly matrix and multiplying an emission
only on observed weeks. There are only 53 distinct week-of-year bins, so the 53
generators are exponentiated once (batched) and indexed.

**Frozen rates are not prevalence.** A subtle but material point: the ratio
$u(w)/(u(w)+d(w))$ is only the equilibrium week $w$ *would* reach if its rates were
frozen forever. Actual prevalence has memory and must be **propagated**: the periodic
occupancy is the fixed point of

$$\pi_w = \pi_{w-1}\,P_{\text{woy}}(w) \quad \text{cyclically over the year}.$$

The two differ substantially here — the frozen ratio peaks ~3 weeks earlier and
overstates the seasonal swing. For $K>2$ the frozen equilibrium must also use the
birth–death form $\pi_k \propto \prod_{i<k} u_i/d_i$, not the two-state formula.

### 4.5 The factorial model: separating biomass from toxin quota

**Why.** pDA is *not* a toxicity readout. Since $\text{pDA} \approx B\cdot Q$, using it
as an observation of a toxicity axis would smuggle biomass variance into that axis.

**The key structural fact** is that in logs the product is additive:

$$\log \text{pDA} = \log B + \log Q,$$

so with two latent chains $B_t$ (biomass, log-mean $\mu_B$) and $Q_t$ (quota, log-offset
$\nu_Q$):

- **cells** observe $\mu_B[b]$;
- **pDA** observes $\mu_B[b] + \nu_Q[q]$ — *this offset is what identifies $Q$*;
- **tissue** observes a censored function of exposure $e = \mu_B + \nu_Q$.

**Empirical support for a distinct $Q$** (SCW): the regression of $\log$ pDA on
$\log$ cells has slope **0.733**, not 1.0 (so not pure proportionality); correlation
0.742 (related but not redundant); $\log_{10}$ quota spans **4.3 orders of magnitude**;
and the variance decomposition is
$\mathrm{var}(\log \text{pDA}) = 0.864$, $\mathrm{var}(\log \text{cells}) = 0.884$,
$\mathrm{var}(\log Q) = 0.451$. The quota carries roughly half the variance of biomass:
it is a major independent axis, not noise.

**Tractability.** Two independent chains have a Kronecker-sum generator, so the joint
transition factorises:

$$Q_{\text{joint}} = Q_B \oplus Q_Q \;\Longrightarrow\; P_{\text{joint}}(g) = P_B^{\,g} \otimes P_Q^{\,g},$$

which lets the existing forward filter run unchanged on the $K_B K_Q$ joint state.

### 4.6 Tissue burden: one-compartment kinetics

Following Mafra et al. (2010), tissue DA is a **stock**, produced by ingesting
particulate toxin and removed by approximately first-order elimination:

$$\frac{dH}{dt} = \alpha E(t) - k H(t), \qquad E = B\cdot Q .$$

Holding $E$ constant across an interval of length $D$ days gives the exact update used
in the code:

$$H_{t+D} = H_t e^{-kD} + \frac{\alpha E_t}{k}\big(1 - e^{-kD}\big).$$

Observations are lognormal with **left censoring at each sample's own LOD**.

**Choices justified by the paper:**

| finding | implementation |
|---|---|
| particulate toxic cells are the principal DA source | uptake driven by $E=BQ$, not biomass alone |
| elimination ≈ exponential | first-order persistence $e^{-kD}$ |
| two-compartment model improved $r^2$ by only 1–9% | **one** compartment |
| no clear mussel body-size effect; CDPH records no sizes | **no** size covariate |
| substantial individual variability | explicit observation noise $\sigma_{\text{tissue}}$ |
| species differ enormously | single-species (*M. californianus*) model; do not pool bivalves |

$k$ is **estimated with a broad prior** ($\mathrm{Lognormal}(\log 0.4,\,0.8)$, spanning
≈0.1–1.5 d⁻¹) rather than fixed at a laboratory value: Mafra measured *M. edulis* at
1.4–1.6 d⁻¹ and cites *M. californianus* at 0.3–0.5 d⁻¹ (Whyte et al. 1995), but field
depuration may be slower than under clean-diet laboratory conditions.

---

## 5. Validation methodology

Nothing was trusted on real data before it recovered known truth on simulated data.

### 5.1 Numerical primitives

The matrix exponential is a differentiable scaling-and-squaring Taylor
implementation, $\exp(A) = \big(\exp(A/2^s)\big)^{2^s}$, needed because NUTS requires
gradients through it. Checked against `scipy.linalg.expm` over hundreds of random
birth–death generators: **max absolute error 1.4×10⁻¹³** (numpy) and **1.7×10⁻¹³**
(PyTensor), with all outputs valid stochastic matrices.

### 5.2 Parameter recovery under the real design

| model | truth → estimate | diagnostics |
|---|---|---|
| stationary hurdle HMM | up 0.050→0.061, down 0.300→0.363, $\mu$ [0.3,3.5]→[0.36,3.48], $\sigma$ 0.80→0.735, $\pi^0$ [.55,.05]→[.52,.095] | $\hat R$ 1.007, 0 div |
| seasonal $Q$ | onset cycle recovered; rate range 0.060–0.303 → 0.062–0.351/wk; peak week 4 → 6 | $\hat R$ 1.008, 0 div |
| tissue kinetics (fast) | $k$ 0.40→**0.382** [0.33,0.45], $\alpha$ 4.0→3.61 | $\hat R$ 1.003, 0 div |
| tissue kinetics (slow) | $k$ 0.10→**0.105** [0.093,0.119], $\alpha$ 1.0→1.02 | $\hat R$ 1.004, 0 div |

The tissue recovery deliberately used the **real** exposure path, **real** sampling
dates and **real** detection limits, simulating only the tissue values — so it tests
identifiability of the actual field design. Both fast and slow depuration are
recoverable, which is what licenses interpreting the fitted $k$ as data-driven.

### 5.3 Likelihood reconstruction

The fitted likelihood is a single lumped `pm.Potential`, so PyMC stores no
per-observation log-likelihood. `selection.py` reconstructs the one-step predictive
densities $\log c_t$ (§2.4) by replaying the same recursion vectorised over posterior
draws. Verified against an independently written scalar forward filter:

- stationary model: agreement to **4.5×10⁻¹¹**;
- seasonal model: agreement to **4.4×10⁻¹⁶**.

Their per-draw sum equals the fitted log-likelihood by construction.

### 5.4 The LFO implementation

The Pareto-smoothing shape estimator (Zhang & Stephens) was unit-tested against known
truth: $k_{\text{true}} = 0.1/0.3/0.5/0.7/1.0 \to$ estimated
$0.10/0.30/0.49/0.69/1.01$. The full LFO routine was checked on synthetic
conditionally-i.i.d. data, where it should reduce to the naive mean log predictive:
$-201.2$ vs $-201.5$, max $\hat k = 0.01$, zero refits.

### 5.5 Convergence gates

Selection considers only fits with $\hat R < 1.05$ and zero divergences. Model
comparisons are never computed or saved across unconverged fits. `p_loo` far
*exceeding* the parameter count is treated as a misspecification flag (Vehtari).

---

## 6. Current results

### 6.1 State count (SCW)

$K=4$ is the **only** value that converges, and is also best by LOO — bracketed by
misspecification below and non-identifiability above:

| $K$ | elpd | $p_{\text{loo}}$ (params) | max $\hat R$ | verdict |
|---|---|---|---|---|
| 2 | −907.0 | 64.2 (15) | 1.534 | misspecified, unmixed |
| 3 | −838.4 | 34.3 (23) | 1.554 | misspecified, 4 divergences |
| **4** | **−802.1** | **17.4 (31)** | **1.007** | **converged, healthy** |
| 5 | −815.1 | 25.5 (39) | 2.148 | over-parameterised; a ~11-cell state with 26% closure (impossible) |

$K=4$: min ESS 543, max Pareto $\hat k$ 0.34, 0 divergences.

### 6.2 The factorial model — the central scientific result

Fitted SCW, $K_B=3$, $K_Q=2$, continuous censored tissue ($\hat R$ 1.008, 0 div,
$\sigma_{\text{tissue}} = 1.10$):

biomass levels ≈ 162 / 1,491 / 39,862 cells; quota levels **9.4× apart**;
high-quota regimes persist ~7 weeks and hold ~31% of the time.

| biomass | $Q$ low: burden → P(closure) | $Q$ high: burden → P(closure) |
|---|---|---|
| ~162 cells | 0.04 ppm → 0.0% | 0.09 ppm → 0.0% |
| ~1,491 cells | 0.08 ppm → 0.0% | 0.50 ppm → 0.1% |
| **~39,862 cells** | **1.34 ppm → 0.8%** | **12.40 ppm → 33.0%** [24–43%] |

**At identical biomass, tissue burden differs 9.3× and closure risk ~40×.** Biomass
does not determine toxicity. Since a remote sensor observes $B$ only, this gap is an
**irreducible ceiling on the value of any biomass-only sensor**, and it is the
mechanism behind the expected ordering $\mathrm{EVPI}(\text{biomass}) < \mathrm{EVPI}(\text{toxicity})$.

Two structural features matter for the decision layer:

1. **Asymmetry.** At low and medium biomass, risk is ~0 *regardless of quota*. Biomass
   is decisive for ruling danger **out**, and nearly uninformative for ruling it **in**.
2. **Closure is a tail event.** Even in the worst state the *median* burden (12.4 ppm)
   is **below** the 20 ppm threshold. Closures occur because of the spread
   ($\sigma_{\text{tissue}} = 1.10$), not because the typical bad week crosses the line.
   A decision rule must therefore target tail mass, not the mean.

Robustness: $\mu_B$ and the quota spread barely moved when the tissue link was
completely respecified (ordered logit → continuous censored lognormal), so the $B$/$Q$
separation is not an artefact of the tissue model.

### 6.3 Tissue kinetics

Fitted on SCW: $k = 0.448$ /day [0.230, 1.086] — **half-life 1.55 days**, landing on
the literature range for *M. californianus* despite the broad prior. Weekly retention
$e^{-7k} = 0.043$.

**Model comparison (LFO): kinetic vs no-memory = +0.8 ± 1.3 — indistinguishable.**
The reason is structural: tissue memory (~1.5 days) is *shorter than the sampling
interval* (7 days), so at weekly resolution the stock is at quasi-equilibrium with
current exposure. This is a measurement, not an identification failure — §5.2 shows
slow depuration would have been detected had it been present.

**Consequence:** $H_t$ should **not** be a separate POMDP state on this data. It is
correctly represented as a function of current exposure, which is what the continuous
tissue link encodes.

### 6.4 Seasonality, re-scored honestly

| estimator | seasonal − stationary |
|---|---|
| in-sample PSIS-LOO (biased) | +14.3 ± 5.5 (2.6 SE) |
| **PSIS-LFO-CV (unbiased)** | **+5.0 ± 2.8 (1.8 SE)** |

Seasonality **survives but is ~3× weaker** than the in-sample estimate suggested. It
should be described as *supported but not established*. Further, the ±2.8 comes from
pointwise differences that are autocorrelated in a time series, so true uncertainty is
likely **wider**. The earlier "az.compare weight 0.96" framing was overconfident and
must not be quoted.

---

## 7. Errors found and corrected

Recorded because a defensible paper must show what was wrong and how it was caught.
Identifiers refer to an external code review.

| id | defect | correction |
|---|---|---|
| **C-3** | per-state-$\sigma$ tobit made severity non-monotone: a bigger bloom could read as a lower state | hurdle + shared $\sigma$ (§4.3); monotonicity now holds by construction |
| **C-2** | PSIS-LOO on forward-filter factors reported as predictive validation | implemented PSIS-LFO-CV (§2.7); seasonality re-scored and corrected downward |
| **C-4** | saved comparison table ranked an **unconverged** $K=4$ first with weight 1.0 | comparisons now computed only across converged fits; stale artefact deleted |
| **H-7** | frozen-rate ratio $u/(u+d)$ labelled as "bloom occupancy" | propagated periodic occupancy (§4.4); frozen curve retained only as a labelled diagnostic |
| **H-9** | a same-data HMM smoother described as "independent validation" | relabelled as internal consistency |
| **H-10** | observation-weighted posterior mean labelled "= stationary occupancy" | relabelled; they are different quantities |
| **H-8** | seasonal priors hard-coded the stationary posterior means (empirical-Bayes double use) | neutral weakly-informative prior $\mathcal N(\log 0.1, 1)$ |
| **H-5** | one generator/emission pooled across four sites; state acted as a site proxy | anchored on SCW alone (§3.4) |
| **H-1** | tissue matched by exact site name, silently dropping records | coordinate matching (§3.2) |
| **H-4** | tissue detection limits collapsed into one "non-detect" category | continuous censored lognormal using each sample's own LOD |
| **M-5** | ad-hoc, candidate-dependent $p_{\text{loo}}$ gate | standard criteria: $\hat R$ and divergences |
| **M-7** | seasonal trace saved before likelihood attached; exported regardless of convergence | likelihood attached first; export refuses unconverged fits |
| **M-9** | forward filter could underflow to $\log 0$ | floor guard in both filters |
| **M-10** | referenced `docs/methodology.md` did not exist | this document |

Two claims of my own were also corrected by later evidence and should not be reused:
the pooled "$K\ge3$ are non-identifiable" conclusion (an artefact of the broken
emission), and the toxic-state closure probability, which fell from 51% → 40.6% → **33%**
as the tissue link became continuous and censoring-aware.

---

## 8. Known limitations

1. **No decision layer.** No actions, rewards, belief-based policy, or VoI computation
   exists. No current claim about sensor break-even accuracy is supported by this code.
2. **Closure evidence is thin.** 22 closure observations at SCW; only 20 closure
   *matchups* pooled across all stations. The 33% estimate carries [24–43%].
3. **$Q$ conflates toxin-per-cell with species composition.** Not all *Pseudo-nitzschia*
   are toxigenic, so $Q$ is best read as "toxigenic propensity" rather than strictly
   toxin quota.
4. **Missing-at-random is assumed, not audited** (§2.5).
5. **Single site.** Results are SCW-specific by design; generalisation is untested.
   *M. californianus* only.
6. **In-sample LOO figures elsewhere in the artefacts** (the K-scan elpd column) remain
   in-sample and must not be quoted as forecast skill.
7. **Week convention.** `to_period("W-SUN").start_time` yields *Monday* week starts;
   timestamps are timezone-naive. Harmless internally, but must be made explicit before
   any daily satellite alignment.
8. **SCW cells are *seriata*-only** and not comparable to other sites' totals.

---

## 9. What to do next

**Immediate (completes the world model).**

1. **Re-score the $K$-scan with LFO** before quoting any elpd gap predictively. The
   $K=4$ choice itself rests on convergence, not LOO, so the selection stands — but the
   gaps should not be cited as forecast skill.
2. **Test $K_Q = 3$**, to check whether the quota ladder is finer than low/high.
3. **Audit the sampling process** for informative missingness (§2.5) by site, season,
   prior result, and programme era.

**The decision layer (the actual research question).**

4. **Write the decision contract first**: who decides, what action set (do nothing /
   order a tissue test / close / reopen), the exact within-epoch timing, the baseline
   monitoring policy, the horizon, and the loss function. Song et al.'s
   observe → update → act ordering is the right template; the reward should depend on
   $P(H \ge 20)$, which §6.2 already produces.
5. **Add the remote-sensing observation channel** as a noisy observation of $B$ only,
   with a sweepable accuracy parameter. Keep it strictly in $O$ (§1.2).
6. **Compute VoI as tail-risk reduction**, not mean-risk reduction (§6.2, point 2), and
   decompose the value loss into *sensor error* versus *biomass–toxicity decoupling* —
   the latter is the irreducible ceiling and is the paper's central quantitative claim.

**Scientific extensions, in order of expected value.**

7. **Test whether SST/upwelling predicts $Q$** (not $B$). This is the high-leverage
   covariate question: $Q$ is what a biomass sensor cannot see, so a satellite-observable
   predictor of *toxigenicity* would materially raise the achievable VoI. Adding
   covariates to predict $B$ is expected to add little, since they largely duplicate the
   seasonal term.
8. **Sub-weekly tissue sampling during events** would resolve the kinetics that weekly
   sampling cannot (§6.3) — a monitoring-design recommendation that follows directly
   from the measured 1.55-day half-life.
9. **Extend to a second site** (Morro Bay Back Bay has 166 matchups) with site-specific
   or hierarchical parameters, to test whether the $B$/$Q$ structure generalises.

---

## 10. Code map and reproduction

| module | role |
|---|---|
| `config.py` | anchor site, match radius, thresholds, column maps |
| `data.py` | ingestion, haversine station matching, weekly and exact-date panels |
| `worldmodel/linalg.py` | differentiable (batched) matrix exponential, birth–death $Q$ |
| `worldmodel/emissions.py` | hurdle cell emission, categorical DA emission |
| `worldmodel/model.py` | sequence prep, stationary HMM, scaled forward filter |
| `worldmodel/seasonal.py` | harmonic $Q(\text{woy})$, grid forward filter, propagated occupancy |
| `worldmodel/factorial.py` | factorial $(B,Q)$ model, continuous censored tissue link |
| `worldmodel/tissue.py` | one-compartment kinetics on exact dates |
| `worldmodel/selection.py` | per-observation log-likelihoods, forward–backward smoother |
| `worldmodel/lfo.py` | PSIS-LFO-CV and Pareto smoothing |
| `worldmodel/export.py` | world-model export for the decision layer |
| `viz.py` | coverage heatmap, seasonal figures |

```bash
python scripts/00_build_panels.py          # raw -> interim panels
python scripts/03_match_radius_scan.py     # spatial matching, 3 km vs 5 km
python scripts/01_fit_worldmodel.py 2 3 4 5  # K-scan with convergence gating
python scripts/02_fit_seasonal.py 1        # seasonal Q at the production K
```

**Reproducibility gaps to close:** no test suite, no environment lock, no artefact
manifest distinguishing current from stale outputs. These were flagged in review and
remain open.
