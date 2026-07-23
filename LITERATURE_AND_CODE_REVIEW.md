# From HAB observations to decision-grade Value of Information

## A deep literature review and forensic review of the `voi-remsen` repository

**Review date:** 21 July 2026

**Repository snapshot:** the current working tree at `/Users/angelinale/voi-remsen`, including the uncommitted seasonal-model files and artifacts present at review time

**Scope:** seven supplied publications/manuals; all source modules and scripts; all three Arcieri reference notebooks, including their saved outputs; input/interim data structure; saved fitted-model and model-selection artifacts; and the four locally stored C-HARM v3.1 NetCDF products

**Review type:** read-only methodological, statistical, logical, and reproducibility review. No model was altered, no new implementation was written, and the expensive models were not re-fit.

---

## Abstract

This repository is intended to quantify the value of remote-sensing information for harmful algal bloom (HAB) decisions. At the reviewed snapshot, it does not yet implement that decision problem. It implements a data-preparation pipeline and two related hidden-state ecological models: a stationary continuous-time hidden Markov model (HMM) and an experimental seasonally varying extension. These models combine weekly *Pseudo-nitzschia* cell counts with California Department of Public Health shellfish-tissue domoic-acid measurements at four sites. A separate R model fits tissue toxicity alone. The repository also contains three older Arcieri tutorial notebooks about Bayesian HMM inference for Markov chains and action-conditioned transition models.

The strongest part of the current work is its structural handling of irregular gaps through a continuous-time generator and its explicit propagation of posterior parameter draws. The stationary two-state sampler appears numerically converged, and the custom matrix exponential is accurate in the fitted parameter range. The current model nevertheless falls well short of a decision-grade HAB simulator. Its single latent severity variable is driven primarily by cell biomass, not closure-level toxicity; its high-biomass state still emits a tissue non-detect roughly 77% of the time. It treats same-week water-cell and shellfish-tissue measurements as conditionally independent contemporaneous manifestations of one state, despite biological accumulation and depuration lags. It pools sites, decades, assay regimes, and highly unequal sampling processes into one set of transitions and emissions. It also converts a historically changing tissue detection limit into one undifferentiated “non-detect” category and constructs total *Pseudo-nitzschia* counts in a way that can turn a missing taxonomic component into an implicit zero.

Two evaluation defects are especially consequential. First, the repository applies ordinary pointwise PSIS-LOO to one-step forward-filter likelihood factors. Removing one factor does not remove that observation’s effect on the filter state and on later likelihood factors; the resulting score is neither genuine leave-one-observation-out validation nor leave-future-out forecasting validation. Second, the seasonal visualization calls the frozen-rate equilibrium ratio \(q_{01}(w)/(q_{01}(w)+q_{10}(w))\) the seasonal “bloom occupancy.” In a time-varying Markov chain this is only the equilibrium that would apply if that week’s rates remained frozen. A 53-bin plug-in cycle propagated through the exported weekly matrices lags that curve: the plotted ratio peaks at week 17 while the propagated approximation peaks at week 20. Because the fitted calendar mapping uses actual day-of-year bins, a 52.1786-week harmonic period, leap years, and occasionally repeated or skipped bins, even that propagated cycle is a diagnostic approximation rather than the exact calendar-specific occupancy of every fitted sequence.

The saved higher-state models do not converge: maximum \(\hat R\) is 1.544 for \(K=3\) and 1.536 for \(K=4\). Although the export script excludes them, the saved comparison table still ranks \(K=4\) first with weight 1.0. That artifact must not be treated as scientific evidence. The converged \(K=2\) fit is an exploratory biomass-regime model, not a validated toxic-risk state model.

The supplied literature provides useful components but not a turnkey solution. Arcieri et al. show how an action-conditioned HMM, posterior uncertainty, belief filtering, and control can be connected; their papers and tutorial notebooks also contain important mathematical and implementation problems. Song et al. give the cleanest observe-update-act timing and a useful optional-information comparison, but their one-way aging layers are inappropriate for recurrent seasonal HABs. Kim et al. demonstrate a covariate-driven network HMM for blooms but not a decision model, and their validation weakens substantially for rare high thresholds. Jackson supplies the most disciplined treatment of continuous-time panel models, irregular gaps, misclassification, and informative observation times. Bashar and Torres-Machi demonstrate an empirical satellite observation channel inside a POMDP, but their scalar “accuracy” thresholds are artifacts of a particular synthetic channel and cost structure, not universal sensor break-even values.

The project can become a defensible VoI study, but several missing layers must be built and independently validated: a precisely defined management decision and horizon; a spatial/temporal C-HARM matchup; a calibrated remote-sensing observation model; an ecologically adequate, site-aware and lag-aware latent process; explicit actions and their causal or physical effects; a reward model; a belief update; a policy solver; and a nested monitoring comparison that distinguishes gross information value from fixed deployment cost. Until those layers exist, the repository supports exploratory world-model research only.

## Contents

1. A primer for readers new to the subject
2. What a complete HAB VoI system would contain
3. Deep review of the supplied literature
4. Repository anatomy
5. Data pipeline and scientific implications
6. Stationary, tissue-only, and seasonal world models
7. Deep audit of all Arcieri notebook cells
8. Literature-to-implementation crosswalk
9. Consolidated finding register
10. Research path to a defensible HAB VoI study
11. Presently defensible and indefensible conclusions
12. Appendices: execution flow, numerical audit, reproducibility, glossary, and sources

---

## Executive verdict

### The short answer

The code currently answers:

> “Can two sparse monitoring streams be summarized by a shared, ordered, two-state weekly latent regime with continuous-time transitions, and do those inferred transitions appear seasonal?”

It does **not** currently answer:

> “How valuable is C-HARM information for an early-harvest, closure, sampling, or other management decision?”

No checked-in code reads the C-HARM scientific variables, matches C-HARM pixels to shore sites, estimates a C-HARM observation model, defines management actions or rewards, updates a POMDP belief after a C-HARM observation, solves a policy, or calculates VoI. The four C-HARM NetCDF files are present, but they are disconnected data assets. The current reference README describes missing modules, notebooks, and break-even results that are not present in this snapshot ([reference README, lines 11–36](/Users/angelinale/voi-remsen/notebooks/reference/README_arcieri.md:11)).

### Maturity assessment

| Layer needed for a VoI study | Present? | Assessment |
|---|---:|---|
| Raw field-monitoring ingestion | Yes | Functional, but important missingness, assay-era, taxonomic, timing, and site-comparability issues remain. |
| C-HARM ingestion and site/pixel matchup | No | NetCDF files exist locally; no source code uses their HAB variables. |
| Latent ecological world model | Partial | A pooled one-dimensional HMM exists. The converged two-state fit is mainly a biomass-regime model and is not independently validated. |
| Nonstationarity | Experimental | A one-harmonic seasonal transition model exists in the uncommitted working tree; its evaluation and occupancy interpretation are flawed. |
| Remote-sensor observation model | No | There is no \(P(\text{C-HARM output}\mid\text{latent hazard})\), calibration analysis, or uncertainty model. |
| Actions and action effects | No | The current HMM is uncontrolled; no harvest, closure, testing, or mitigation action changes the system. |
| Reward/cost model | No | No public-health loss, harvest value, testing cost, closure cost, false-alarm cost, or discounting is encoded. |
| Belief-state filter for deployment | No | Filtering machinery exists internally for model likelihood/smoothing, but no decision-layer interface or exported initial belief exists. |
| POMDP policy solver | No | No dynamic programming, point-based solver, QMDP, RL agent, or stopping policy is present. |
| VoI calculation | No | No with-information versus without-information policy comparison is implemented. |
| Decision validation | No | There is no held-out policy evaluation, operational comparator, or uncertainty analysis for decision value. |

### Highest-priority findings

The findings below use four labels:

- **Confirmed defect:** the implementation or stated result is internally inconsistent or mathematically wrong.
- **Unsupported inference:** the calculation may run, but the interpretation exceeds what it establishes.
- **Modeling assumption:** potentially defensible, but it needs explicit scientific justification and sensitivity analysis.
- **Unresolved risk:** the relevant implementation is absent, so the literature raises an audit question rather than proving a code defect.

| Priority | Finding | Classification | Consequence |
|---|---|---|---|
| Critical | The current repository contains no C-HARM observation model, POMDP, reward, policy, or VoI computation. | Confirmed scope gap | Any current claim about C-HARM break-even accuracy or decision value is unsupported by this code snapshot. |
| Critical | Ordinary PSIS-LOO is applied to dependent forward-filter factors and presented as model-selection evidence. | Confirmed evaluation defect | The \(K\)-ranking and stationary-versus-seasonal comparison are not valid estimates of future predictive performance. |
| Critical | Nonconverged \(K=3,4\) fits remain in a table that ranks \(K=4\) best with weight 1.0. | Confirmed artifact defect | A reader can easily draw the opposite conclusion from the one warranted by diagnostics. |
| Critical | The one-state-axis model conflates water biomass and tissue toxicity despite different biology and lags. | Modeling assumption with strong contrary evidence | The inferred “bloom” state is not a reliable closure-risk state, so a policy built on it could act on the wrong hazard. |
| High | Missing one *Pseudo-nitzschia* taxonomic group is implicitly treated as absence when the other group is present. | Confirmed data-construction behavior | Cell totals change meaning across sites and eras, creating artificial state/site effects. |
| High | All tissue non-detects are collapsed into one category although the recorded detection limit shifts from about 1.0 to 2.5 µg/g. | Confirmed information loss | Assay-era changes can be misread as ecological changes and the likelihood cannot use censoring correctly. |
| High | Same-week cell means and maximum tissue DA are treated as contemporaneous conditionally independent emissions. | Modeling assumption | Bioaccumulation/depuration lag and within-week timing are erased. |
| High | One transition process, emission system, and initial distribution are shared across four sites and decades. | Modeling assumption | Geographic, protocol, sampling, and calendar effects can be absorbed into latent states and transition rates. |
| High | The seasonal plot labels a frozen-rate equilibrium target as actual bloom occupancy and derives an “onset peak” annotation from it. | Confirmed interpretation defect | Relative to the 53-bin propagated plug-in cycle, the reported prevalence phase is shifted by about three weeks; the week-17 onset annotation happens to coincide with the maximum fitted onset rate here but is not established by the plotted calculation. |
| High | The exported world models omit the learned initial-state distribution and essential provenance/transform metadata. | Confirmed interface defect | A decision layer cannot reproduce the fitted filtering problem safely from the export alone. |
| High | Arcieri’s reference notebooks contain global-variable leakage, a forward-algorithm off-by-one error, and an ignored deterioration-mean parameter. | Confirmed notebook defects | They are unsuitable as trusted templates without correction and revalidation. |
| Medium | Weekly dates are described/configured as `W-SUN`, but pandas’ period start is Monday. | Confirmed documentation/timing mismatch | A future satellite matchup can be shifted at week boundaries unless the convention is made explicit. |
| Medium | Monitoring times are treated as non-informative. | Unchecked assumption | Event-triggered or risk-triggered sampling can bias both transition and emission estimates. |
| Medium | No tests, complete methodology document, environment lock, or artifact manifest is present. | Confirmed reproducibility gap | Results are difficult to reproduce, audit, or distinguish from stale outputs. |

---

## 1. A primer for readers new to the subject

### 1.1 The ecological variables are related, but they are not interchangeable

This project contains at least four conceptually distinct quantities:

1. **Water-column *Pseudo-nitzschia* cell abundance.** This measures how many cells from two microscopy size groups were counted per liter. High cell abundance can signal a bloom, but not every bloom is highly toxic.
2. **Domoic acid in the water.** Particulate and cellular domoic-acid products concern toxin associated with phytoplankton. Toxicity depends on species/strain, physiology, environmental conditions, and toxin per cell, not cell count alone.
3. **Domoic acid in shellfish tissue.** Mussels and oysters integrate exposure through uptake and depuration. Tissue concentration can lag, persist after a water-column signal falls, and differ by species and location.
4. **C-HARM outputs.** C-HARM v3.1 does not directly observe a hidden “safe/unsafe” truth. The local NetCDF metadata describe daily 3 km model probabilities of: *Pseudo-nitzschia* exceeding 10,000 cells/L, particulate DA exceeding 500 ng/L, and cellular DA exceeding 10 pg/cell. These are three correlated model products built from gap-filled satellite ocean color and ocean-model inputs.

A defensible decision model must state which of these is the **hazard state**, which are **observations**, which may be **predictors of transitions**, and which outcome triggers the management cost. Treating all of them as noisy readings of one scalar severity is convenient, but it is a substantive scientific hypothesis, not a neutral simplification.

The initial C-HARM skill paper explicitly states that the optical signals of *Pseudo-nitzschia* abundance and DA could not be isolated directly from bulk chlorophyll using the available multispectral sensors; C-HARM therefore combines remote sensing, hydrographic modeling, and empirical ecological models (Anderson et al., 2016, Introduction, manuscript pp. 2–5). It also reports lead-lag structure: modeled particulate-DA probability led Santa Cruz mussel-tissue DA by about five days, while SPATT and tissue measures could be separated by much longer intervals (Results/Discussion, manuscript p. 24). That evidence directly argues against a same-week, conditionally independent, common-state emission model. See the [NOAA-hosted accepted manuscript](https://repository.library.noaa.gov/view/noaa/33076/noaa_33076_DS1.pdf) and [DOI record](https://doi.org/10.1016/j.hal.2016.08.006).

### 1.2 Markov chains, hidden states, and beliefs

A Markov model says that the next system state depends on the information in the current state, rather than on the full past. For a discrete state \(S_t\), a transition matrix contains

\[
T_{ij}=P(S_{t+1}=j\mid S_t=i).
\]

An HMM adds the fact that \(S_t\) is not observed directly. Instead the system emits a measurement \(Y_t\) according to an observation or emission model,

\[
P(Y_t\mid S_t).
\]

The analyst carries a **belief**, a probability distribution over possible hidden states. A correct predict-update step is

\[
\tilde b_{t+1}(j)=\sum_i b_t(i)T_{ij},
\qquad
b_{t+1}(j)=
\frac{P(Y_{t+1}\mid S_{t+1}=j)\tilde b_{t+1}(j)}
{\sum_k P(Y_{t+1}\mid S_{t+1}=k)\tilde b_{t+1}(k)}.
\]

The belief is more useful than a hard decoded label. A label such as “state 2” hides whether the model is 51% or 99.9% certain; those situations should often produce different decisions.

### 1.3 Continuous-time models and irregular observations

The current repository models latent movement in continuous time with a generator \(Q\). Off-diagonal \(q_{ij}\) values are instantaneous transition rates and each diagonal makes its row sum to zero. The transition matrix over elapsed time \(\Delta\) is

\[
P(\Delta)=\exp(\Delta Q).
\]

This is useful when observations occur at irregular dates: a four-week gap uses \(P(4)\), not four unrelated transitions or an assumption that the gap never occurred. Jackson’s manual explains this construction, exponential sojourn times, panel-data likelihoods, misclassification models, and the conditions under which observation times may be informative (Jackson, manual pp. 1–13). The repository implements this mathematical core correctly for a stationary generator, apart from its deliberate one-year cap on long gaps.

### 1.4 From an HMM to an MDP and POMDP

An HMM describes evolution and observation. A decision process also needs:

- **Actions**: for example, do nothing, order a tissue test, issue an advisory, close harvest, or reopen.
- **Action-dependent dynamics**: what the action changes, if anything.
- **Rewards or costs**: public-health harm, foregone harvest, false-alarm cost, testing cost, reopening delay, and possibly equity or risk constraints.
- **A horizon and discount convention**: how far ahead decisions matter.
- **A policy**: a rule mapping current information to an action.

In an MDP the state is observed. In a POMDP it is hidden, so the policy normally acts on the belief. A complete POMDP must make the ordering explicit. For a monitoring decision, the clean order is often:

\[
\text{choose whether/how to observe}
\rightarrow
\text{receive observation}
\rightarrow
\text{update belief}
\rightarrow
\text{choose management action}.
\]

Song et al.’s two-substep construction follows this order. Bashar and Torres-Machi instead combine inspection and maintenance into one pre-observation action, so the inspection informs only a later epoch. These are different decision problems and can yield different VoI.

### 1.5 Value of information

Value of information is not simply classification accuracy, mutual information, or correlation. It is the improvement in the **optimized decision objective** when an information source is available:

\[
\mathrm{VoI}
=
\text{optimal value with information option}
-
\text{optimal value without it}.
\]

If the information action is optional and the with-information action set contains every baseline action, gross VoI cannot be negative under exact optimization: the decision maker can ignore the new source. Fixed adoption cost should be subtracted separately. Negative values in a nested optional design indicate Monte Carlo error, solver error, inconsistent objectives, or an implementation problem—not intrinsically harmful optional information.

The same sensor can have high value in one decision setting and almost none in another. Bashar and Torres-Machi find little incremental value when satellite data compete with a 90%-accurate annual survey, but much greater value on otherwise unmonitored roads. Song et al. show non-monotonic value as replacement cost changes. Both results illustrate why a universal “break-even accuracy” is not meaningful without a specific observation matrix, event prevalence, costs, timing, baseline monitoring regime, and policy horizon.

### 1.6 Bayesian parameter uncertainty is different from hidden-state uncertainty

Two uncertainties coexist:

- **State uncertainty:** which ecological state is present now?
- **Parameter/model uncertainty:** what are the transition rates, emission probabilities, lag structure, and site effects?

A belief over state with fixed parameters does not represent uncertainty about the model. Arcieri et al. emphasize posterior parameter draws and later domain randomization. That is valuable, but averaging model-specific optimal value functions is not the same as solving the Bayes-adaptive problem in which new data update both state and parameter beliefs. The distinction matters throughout the supplied literature and for any eventual HAB policy.

---

## 2. What a complete HAB VoI system would contain

The literature, local data, and current code imply the following conceptual chain:

\[
\begin{aligned}
&\text{field, laboratory, environmental, and remote-sensing data}\\
&\quad\downarrow\\
&\text{quality control, spatial/temporal alignment, and censoring model}\\
&\quad\downarrow\\
&\text{latent ecological process with site, season, lag, and uncertainty}\\
&\quad\downarrow\\
&\text{observation models for cells, water DA, tissue DA, and C-HARM}\\
&\quad\downarrow\\
&\text{belief over decision-relevant hazard state}\\
&\quad\downarrow\\
&\text{observe-update-act policy with explicit costs and constraints}\\
&\quad\downarrow\\
&\text{nested policy comparison and VoI uncertainty distribution}.
\end{aligned}
\]

The current repository implements parts of the first four lines for two field/laboratory channels. It stops before the remote-sensor channel, decision state, action, reward, policy, and VoI layers.

---

## 3. Deep review of the supplied literature

### 3.1 Arcieri et al. (2023): Bayesian HMM inference joined to maintenance control

**Source.** Giacomo Arcieri et al., “Bridging POMDPs and Bayesian decision making for robust maintenance planning under model uncertainty: An application to railway systems,” *Reliability Engineering & System Safety* 239 (2023), 109496. [Local PDF](</Users/angelinale/Downloads/Bridging POMDPs and Bayesion Decision Making (Arcieri, 2023).pdf>). Relevant core sections are §§2–6, PDF pp. 2–14; posterior distributions are in Appendix B, pp. 14–16.

#### Purpose and contribution

The paper connects three tasks that are too often treated separately:

1. infer hidden condition, transition behavior, and an observation process from monitoring and logged maintenance;
2. retain a posterior distribution over uncertain transition and observation parameters; and
3. use that posterior in an MDP/POMDP maintenance calculation.

That architecture is directly relevant to HAB management. It insists that the ecological “world model” and the information model be inferred rather than assumed, and that parameter uncertainty reach the decision layer. Those are sound principles.

The data consist of 62 Swiss Federal Railways track sequences from 2008–2018, each described as roughly 20 half-year observations and 20 actions after aggregating track geometry over 150 m. The paper states that the selected tracks were maintained in 2019 (pp. 4 and 7), meaning the sample is not representative of all track. Four latent condition states are used: very good, good, bad, and very bad. Actions are do nothing, tamping/minor repair, and renewal/major repair. A continuous negative “fractal value” is the observation. Costs are elicited in abstract units (Table 1, p. 5), with condition cost increasing sharply in the worst state.

#### Statistical model

The transition system is action-conditioned:

\[
s_t\mid s_{t-1},a_{t-1},T
\sim \operatorname{Categorical}(T_{a_{t-1},s_{t-1},:}),
\]

with Dirichlet priors on the rows of \(T\) (Eq. 10, p. 6). The priors encode strong domain expectations: deterioration under no action, improvement under repair, and low probability for physically implausible moves. That can stabilize sparse rows, but it also means the result is not “purely learned from data.” The paper does not report the numerical concentration hyperparameters, which prevents exact scrutiny of prior influence.

The observation model is autoregressive, bounded above by zero, heavy-tailed, and action-specific (Eqs. 11–14, p. 6). Under no action, the change in observation follows a truncated Student-\(t\) whose state-specific mean and scale describe deterioration. Under repair, the next observation depends on a fraction of the prior observation plus a state-specific repair effect. This is much richer than the current HAB repository’s conditionally independent contemporaneous emissions. It recognizes both memory and a distinct intervention response.

The reported posterior means are plausible in broad ordering: no-action transitions are mostly persistent; tamping and renewal increasingly move bad states upward; emission locations become more negative and variable in worse states (Figs. 7–9 and B.15–B.17). The paper used four chains, 4,000 warm-up iterations and 3,000 retained draws per chain. It reports no divergences and visually satisfactory chains, but provides no numerical \(\hat R\), effective sample size, Monte Carlo error, or quantitative held-out predictive score. Validation is essentially one observed-versus-simulated deterioration example and one repair example (Figs. 10–11).

#### Decision calculation

For a fully observed MDP, the authors solve each posterior transition draw separately and average the resulting model-specific optimal Q-functions. The posterior-averaged policy differs from the policy using posterior-mean parameters in state 2: it chooses minor rather than major repair (Table 2, p. 9). In simulation, the posterior-averaged policy has a slightly better mean return than the posterior-mean policy, but a worse lower endpoint of its 95% interval (Table 3). This is expected-value optimization, not tail-risk or worst-case robustness.

For partial observability, the paper uses QMDP (Eq. 16, p. 10). QMDP values each action as if state uncertainty will disappear after the next step. It is computationally convenient and can be useful when observations are frequent and informative, but it cannot value actions whose purpose is to gain information. The authors acknowledge this limitation. The partial-observability simulations in Table 4 do **not** show the “robust” policy dominating the posterior-mean policy: mean return is \(-14{,}526\) versus \(-14{,}478\), respectively, with Monte Carlo uncertainty large enough that the small difference is not compelling.

#### Confirmed and conceptual problems

1. **The published belief-update pseudocode is wrong as written.** Algorithm 1 on PDF p. 12 computes a sum of transition probabilities over prior states but omits the prior belief weight \(b_t(s_t)\). It therefore uses a transition-matrix column sum rather than a predicted belief. This contradicts the correct Eqs. 5–6. It may be a typesetting/pseudocode error rather than the authors’ executed implementation, but copying the algorithm literally would be mathematically wrong.

2. **The claimed sufficient belief omits observed context.** The likelihood is autoregressive: \(P(z_{t+1}\mid s_{t+1},a_t,z_t)\). Two histories can have the same distribution over \(s_t\) but different current \(z_t\), leading to different future observation distributions. A state-only belief \(b_t(s)\) is therefore not by itself a sufficient Markov information state for the claimed POMDP or for an exact extension. The decision input must include \(z_t\), or the state must be augmented so that the belief summarizes all predictive context. Arcieri’s particular QMDP action calculation ignores future observation dynamics after the next step, so this omission need not change that heuristic’s Q-values; it remains a flaw in the general information-state claim.

3. **Averaging \(Q^*_{\theta}\) is not generally Bayes-adaptive optimal control.** For each parameter draw \(\theta\), \(Q^*_{\theta}\) assumes future actions can follow the policy optimal for that known \(\theta\). A real controller does not know which fixed draw is true. Averaging these clairvoyant continuation values generally exchanges expectation and future optimization in a way that no single deployable policy can realize. A Bayes-adaptive controller would carry and update a posterior over \(\theta\) as part of its state. “Posterior-averaged heuristic” is a more precise description than “robust optimal.”

4. **Logged actions are not automatically causal treatments.** Operators take maintenance action because of observed and unrecorded condition information. The HMM conditions on actions but does not model the behavior policy or unmeasured reasons for action. Repair rows have a causal interpretation only under a strong sequential-ignorability assumption: given modeled state/history, no unmeasured factor affects both action selection and future condition. Confounding by indication and regression to the mean are plausible.

5. **The state count is not selected predictively.** Four states are retained because three- and five-state versions behaved less well computationally and four aligns with domain grades. That can be pragmatic, but it is not evidence that four latent physical states generated the data.

6. **The policy claims remain simulator claims.** Real data inform the fitted simulator. All comparisons between policies then occur inside that simulator and cost model, over a 25-year horizon that extrapolates beyond the roughly ten-year data record. The paper does not demonstrate improved field outcomes.

#### What transfers to this HAB project

The transferable core is the separation of transition inference, observation inference, posterior uncertainty, belief updating, and control. The project should also preserve joint posterior draws rather than independently mixing marginal parameter samples. The parts that should **not** be imported without modification are the published Algorithm 1, the claim that posterior-averaged model-specific optima are Bayes-optimal, the state-only policy input under autoregressive observations, and a causal interpretation of action-conditioned transitions from observational behavior data.

### 3.2 Arcieri et al. (2024): belief-input PPO and domain randomization

**Source.** Giacomo Arcieri et al., “POMDP inference and robust solution via deep reinforcement learning: an application to railway optimal maintenance,” *Machine Learning* 113 (2024), 7967–7995. [Local PDF](</Users/angelinale/Downloads/POMDP Inference & Deep RL (Arcieri, 2024).pdf>). The statistical model is on PDF pp. 8–12; Algorithm 1 and nominal experiments are pp. 13–15; domain randomization is pp. 16–17; hyperparameters are in Appendix B, p. 25.

#### What the paper adds

The 2024 paper largely reuses the 2023 railway HMM and simulator, then compares three PPO policy representations:

- a feed-forward network that receives the explicit belief;
- an LSTM that receives the last three observations/actions; and
- a gated Transformer-XL that receives a longer history.

In the nominal posterior-mean simulator, the exact-belief QMDP benchmark achieves about \(-14{,}374\), belief-input PPO \(-14{,}677\), the coded current policy \(-16{,}295\), GTrXL \(-17{,}196\), and LSTM \(-18{,}167\) over 100,000 simulated trials (Table 2, p. 15). This is useful evidence that handing a reasonably accurate filter to a policy learner can be far more sample-efficient than asking a generic recurrent network to learn filtering and control simultaneously.

The paper then samples a posterior parameter draw at the start of each episode. A single belief-input PPO policy is trained across these environments. Domain randomization improves the simulated mean from roughly \(-14{,}901\) to \(-14{,}648\) (Table 3, p. 17). This supports a modest, in-distribution conclusion: training across the same parametric posterior used for evaluation can improve expected simulator performance.

#### Important qualifications

1. **The architecture comparison is not controlled enough for a general ranking.** The LSTM sees only three steps; GTrXL sees up to 50. Learning rates, PPO clipping values, architecture sizes, and tuning procedures differ. Parameter counts and equalized hyperparameter-search budgets are not reported. The result is about these implementations on this simulator, not a theorem that beliefs beat sequence models.

2. **The selected-best-evaluation procedure is optimistic.** Models are repeatedly evaluated on 500 episodes during training and the best is selected. Repeated selection over noisy scores creates a winner’s-curse effect. Locked validation and test seeds are needed.

3. **The “current real-life policy” comparison remains entirely simulated.** The comparator is a threshold rule coded inside the inferred simulator. No real tracks were controlled by the learned policy, and no causal off-policy evaluation establishes what historical outcomes would have been under it.

4. **The PPO pseudocode is not an adequate on-policy PPO algorithm.** Algorithm 1 on PDF p. 13 initializes a replay buffer, keeps adding \((y_t,a_t,r_t)\), never clears it, and omits old log probabilities, values, advantages, terminal flags, and next observations. PPO ordinarily trains on fresh on-policy rollouts. The actual unpublished training code may use a correct library; the printed algorithm cannot establish that.

5. **There is a serious unresolved oracle-leakage question.** The domain-randomization description says the episode model \(\hat\theta\) is sampled and otherwise follows Algorithm 1, whose belief update uses environment parameters. If the filter receives the exact sampled \(\hat\theta\), it knows the episode’s true model—information unavailable at deployment. A deployable design would use a fixed posterior-mean filter, a posterior-mixture filter, or a joint state-parameter filter. The RL code is not present in this repository, so this is an unresolved risk, not a confirmed implementation defect.

6. **The autoregressive sufficiency problem persists.** A policy that sees only \(b_t(s)\) omits \(z_t\), even though the next-observation distribution depends on \(z_t\). The belief representation is therefore incomplete for the authors’ own observation model.

#### What transfers to HAB

The broad idea—separate ecological filtering from policy learning and train policies across posterior world-model draws—is promising. It should come only after the filter is prospectively calibrated and after the parameter draw used by the simulator is hidden from the filter. For a modest state/action HAB problem, interpretable dynamic programming or point-based methods should be established as strong baselines before deep RL. RL does not repair a misspecified ecological model or missing reward definition.

### 3.3 Arcieri et al. (2026): Deep Belief Markov Models

**Source.** Giacomo Arcieri et al., “Deep belief Markov models for POMDP inference,” *Neural Networks* 196 (2026), 108386. [Local PDF](</Users/angelinale/Downloads/Deep Belief Markov Models for PODMDP (Arcieri, 2026).pdf>). Core definitions and DBMM factorization are PDF pp. 2–5; experiments are pp. 5–10; Appendix A is pp. 11–13.

#### Purpose

The DBMM is intended to replace a hand-specified HMM with neural modules that learn:

- a transition from the previous belief and action to a predicted belief;
- an observation model; and
- an inference/update map from predicted belief and current observation to posterior belief.

It is trained variationally and then used as the belief input to PPO. Experiments cover a discrete bridge-like benchmark, a nonlinear continuous benchmark, and the railway simulator inherited from the 2023 work. The motivation is attractive: a reusable approximate belief operator could make partially observed control more tractable when exact Bayesian filtering is unavailable.

#### What the experiments show

On the discrete benchmark, aggregate cross-entropy becomes close to that of exact beliefs, but per-class results are uneven: one class is never recovered by the DBMM and is barely recoverable even by the exact filter (Table 1, p. 6). On the railway simulator, the DBMM is materially worse than the exact filter—especially for one state—and has KL 1.28 versus 0.89 (Table 2, p. 10). Its PPO curve is nevertheless closer to exact-belief PPO than the LSTM/GTrXL curves in the plotted simulator experiment.

The continuous experiment reports lower mean-squared error than an ensemble Kalman filter (EnKF) and apparently good calibration. The EnKF is not an exact gold standard for a nonlinear one-dimensional problem; numerical quadrature or a well-tuned particle filter would be a more informative benchmark. A Gaussian DBMM belief also cannot represent a multimodal posterior.

The railway DBMM is not trained directly on the 62 field sequences. It is trained on random-policy trajectories generated from the earlier fitted railway simulator, and the DBMM observation family is a truncated Normal approximation to the simulator’s truncated Student-\(t\). Calling this a “real-world” benchmark blurs the distinction between a simulator inferred from real data and direct real-data validation.

#### Fundamental formulation concerns

1. **A continuous-state belief is not generally a finite vector in \(\mathbb R^{|S|}\).** The claim on p. 3 is valid only after choosing a finite-dimensional parametric family, such as Gaussian mean/covariance. A general continuous-state belief is a probability measure or density.

2. **Equation 16 appears unnormalized and double-counts the prior.** Equation 15 already defines \(q_\psi(s_t\mid\tilde b_t,o_t)\) as the posterior belief. Equation 16 then multiplies it by a prior factor \(q_\omega(s_t\mid b_{t-1},a_{t-1})\). The product of two normalized densities over the same variable is not normalized without an explicit normalizing constant. If the second network returns an updated posterior, the factorization should ordinarily contain the updated posterior once, with the predicted belief as an input.

3. **The printed generative process may lose trajectory dependence.** In Eq. 13, \(s_t\) is sampled from a distribution produced by the previous belief, but the realized \(s_t\) does not feed the next transition; the next distribution is again propagated from a belief. This can match marginal beliefs while failing to generate correct serially correlated state trajectories. That is a major concern for claims of high-fidelity simulation.

4. **The generative and inference recursions use different previous beliefs.** The generative transition in Eq. 12 propagates a previous prior, while Eq. 14 propagates the inference model’s previous posterior. The probabilistic joint and ELBO need a more rigorous normalization and dependency derivation.

5. **The inference update omits action as a direct input.** The observation likelihood may depend on action, but \(\mathcal Q_\psi\) takes only predicted belief and observation. If two actions produce the same predicted belief but different observation mechanisms, correct posteriors differ. Action cannot always be compressed through the predicted belief.

6. **Unsupervised latent states are not automatically interpretable physical states.** Label permutation is acknowledged. In continuous models, more general transformations can preserve fit. Without anchors, monotonicity, known semantics, or labeled states, a neural latent coordinate can be useful for prediction while lacking the physical meaning needed by a reward function.

7. **Variational state uncertainty is not epistemic model uncertainty.** A distribution over \(s_t\) with point-estimated neural weights does not retain the 2023 posterior over transition and emission parameters. The 2026 railway experiment omits the prior paper’s model-uncertainty layer.

8. **“Near-optimal” is not established.** The exact-belief comparator is another PPO policy, not the solved optimal POMDP. It is an oracle-input learning baseline.

9. **The PPO pseudocode repeats the accumulating replay-buffer problem.** As in 2024, the printed algorithm is not a complete or clearly on-policy PPO procedure.

Appendix A also contains direct specification problems: Eq. A.3 labels the replacement standard deviation as a deterioration quantity; Eq. A.4 omits how stochastic variance is interpolated across actions; the printed do-nothing dynamics appear inconsistent with Fig. 14; and Eq. A.5 is ambiguous about whether its second Normal parameter is variance or standard deviation.

#### What transfers to HAB

The DBMM suggests a future direction if the ecology truly requires nonlinear high-dimensional filtering. It should not be the first move for this dataset. With only four sites, sparse tissue observations, strong protocol shifts, and a still-undefined decision state, a neural latent model can hide rather than solve identifiability. If explored later, it needs a normalized temporal generative model, explicit action in the update when appropriate, flexible heavy-tailed/multimodal emissions, state-semantic anchors, held-out likelihood and calibration, and a separate treatment of epistemic parameter uncertainty.

### 3.4 Song et al. (2022): nonstationary reliability-assisted POMDP and VoI

**Source.** Chaolin Song et al., “Value of information analysis in non-stationary stochastic decision environments: A reliability-assisted POMDP approach,” *Reliability Engineering & System Safety* 217 (2022), 108034. [Local PDF](</Users/angelinale/Downloads/Non-Stationary POMDP for VoI (Song).pdf>). VoI/POMDP definitions are PDF pp. 3–4; the layered construction and timing are pp. 5–8; the beam example and results are pp. 9–11.

#### Central idea

Song et al. want to use a stationary POMDP solver for a physical deterioration process whose rate accelerates with age. They enlarge the latent state space. Each physical condition is copied into \(k\) deterioration-rate layers. Within a layer, the system remains in the same physical condition or progresses to the next worse condition; it can also move one-way to a faster layer. The enlarged transition system is stationary, while the marginal distribution over physical condition changes nonstationarily as probability moves through layers (Figs. 2–4, pp. 5–6).

Replacement resets every transition row to the initial belief, representing a new component. This is conceptually clean for an aging asset. It is not a natural model for recurring annual HAB risk: season and circulation can move risk both upward and downward, and the regime repeats rather than irreversibly aging through layers.

#### VoI timing

The strongest transferable feature is the two-substep decision epoch (Fig. 5, p. 7): inspection is selected, an observation is received, belief is updated, and maintenance is then selected. This ensures that current information can change the current management action. The paper distinguishes the value of one current inspection from the value of a continuing information flow (Eqs. 8–9).

Its gross VoI is nonnegative because inspection is optional: if inspection costs more than it helps, the optimal policy declines it. That statement would not hold for a mandatory system or after imposing an unavoidable fixed adoption cost.

#### Fitting the transition and observation models

The transition parameters are fitted so the layered Markov model reproduces time-varying physical-state probabilities obtained from a separate structural reliability calculation (Eqs. 10–12, pp. 7–8). This matches marginal prevalence curves, not necessarily conditional path dynamics. Many transition systems can generate the same marginals but different observation value and policies.

The observation matrix is derived from a continuous measurement-error model with Bayes’ theorem and nested Monte Carlo (Eqs. 13–21, p. 8). This is an important principle for the HAB project: the sensor matrix should be connected to an actual measurement or forecast likelihood, rather than set to a symmetric “correct with probability \(p\)” template. The paper is less clear about the exact rule converting a continuous measurement into a discrete signal. An implementable model needs an explicit classifier \(g(y)\) or should retain a continuous observation.

#### Numerical example and results

The example is a corroding beam with three physical conditions, four deterioration layers, inspection every two years, and replace/do-nothing maintenance. Remaining depth and width are measured with independent Normal error (p. 9). The fitted layered process tracks the same Monte Carlo reliability curve used for calibration; this is in-sample reconstruction, not independent validation.

The current inspection has zero value at the initial step, while the continuing information flow has reported value \$3,284: expected negative cost improves from \(-\$9{,}461\) without inspections to \(-\$6{,}177\) with them. VoI falls as measurement error or inspection cost increases. It is non-monotonic in replacement cost: information matters little when replacement is almost free, matters most when inspection can change a borderline replacement decision, and again matters little when replacement becomes uneconomic (Figs. 7–8, pp. 10–11).

#### Errors and limitations

- Equation 12 minimizes an object called a likelihood; maximum likelihood ordinarily maximizes it unless the object is a negative log-likelihood, which is not stated.
- The prose gives failure cost as \(1.16\times10^5\), while Table 2 prints \(1.16\times10^6\). The sensitivity discussion is more consistent with the former.
- Table 1 calls the example a crack-growth model although the equations describe uniform corrosion.
- The definition \(w=\theta\sqrt t\) is called an instantaneous rate, then multiplied by \(t\) as though it were an average rate; the units and integration are inconsistent.
- Relative-error Normal fitting is unstable when target probabilities are zero/small and ignores the simplex dependence between state probabilities.
- The claim that ordinary time augmentation is inapplicable is too strong. Calendar time can be a perfectly observed component of a mixed-observable process, and a time-indexed finite-horizon POMDP is standard.
- Persistent physical parameters create history dependence. A coarse physical state plus a layer may still be non-Markov unless the posterior over those parameters is represented.
- Solver settings, belief-point count, tolerance, reliability Monte Carlo size, and policy-stability analysis are not reported.

#### HAB lesson

Use Song’s observe-update-act ordering and its distinction between one inspection and a continuing information service. Do not copy its one-way aging layers for seasonality. For HAB, calendar season is observed and environmental regimes recur; time/covariate-dependent transitions, a mixed-observable state, or a semi-Markov ecological model are more plausible.

### 3.5 Kim et al. (2022): a covariate-driven network HMM for algal blooms

**Source.** Kue Bum Kim et al., “A multivariate Chain-Bernoulli-based prediction model for cyanobacteria algal blooms at multiple stations in South Korea,” *Environmental Pollution* 313 (2022), 120078. [Local PDF](</Users/angelinale/Downloads/multivariate HAB prediction (Kim, 2022).pdf>). Data and model are PDF pp. 2–3; validation design is pp. 4–6; results are pp. 6–10; metric formulas are p. 11.

#### What it models

Kim et al. analyze bloom/no-bloom occurrence at eight Nakdong River stations over 72 modeled monthly time points from 2013–2018. The source first describes the water-quality inputs as weekly, but its HMM, figures, and six-year sequence use 72 monthly steps without clearly explaining the weekly-to-monthly aggregation or terminology change. A common hidden network state \(H_t\) drives the eight-dimensional binary occurrence vector. Transition probabilities are nonhomogeneous and depend on current water temperature, total phosphorus, total nitrogen, and modeled water velocity through multinomial logistic regression (Eqs. 1–5, p. 3).

This is a potentially useful template for a spatially coordinated ecological regime: sites share a large-scale latent condition while retaining site-specific emission behavior. It is not a POMDP. There are no management actions, rewards, belief-based policies, or VoI.

The transition coefficient for a predictor is indexed by destination state, not separately by origin-destination pair. That constrains a predictor’s effect on entering a destination to be the same regardless of the prior state. Necessary baseline-category constraints are not stated. Equations 1 and 5 and the accompanying prose imply conditionally independent Bernoulli station outcomes, but the branded Chain-Bernoulli emission is not given as an explicit station-wise factorization with its full parameterization and constraints. The main text therefore leaves important implementation detail underspecified.

#### Model selection and validation

BIC selects five hidden states, but the candidate range, BIC values, parameter counts, starting-value sensitivity, and local-optimum checks are not reported. With only 72 network-level time points, five transition states plus covariate and site-emission parameters may be weakly identified.

The authors use leave-one-year-out validation. This is better than randomly splitting months, but five of six folds train partly on years that occur after the held-out year. It is not a pure forward-deployment test. The paper also does not establish whether state count, scaling, medians, and preprocessing were reselected within each fold.

Performance deteriorates sharply for rare, high thresholds. At 10,000 cells/mL, held-out probability of detection ranges 0.333–0.750, false-alarm ratio 0–0.667, and critical success index 0.286–0.500 (Table 1, p. 6). The abstract’s much stronger summary mainly reflects calibration, not the hardest held-out setting. No uncertainty intervals, proper probabilistic scores, calibration curves, or persistence/seasonal baselines are supplied.

The scenario analysis randomly varies predictor pairs while fixing other predictors at medians, simulates labels from the fitted HMM, and uses a linear SVM to describe separation. It explains the fitted model rather than providing independent evidence. Random predictor combinations can also violate the empirical seasonal joint distribution. Associations should not be read as causal effects of manipulating flow, nutrient, or temperature.

Several figure references on p. 5 are swapped, and Fig. 8’s caption describes only the 1,000-cell threshold although all three thresholds appear. The text describes Viterbi as estimating posterior state probabilities; Viterbi finds a most-probable joint path, whereas marginal posterior probabilities require forward-backward inference.

#### HAB lesson

Kim contributes a candidate covariate-driven state estimator, not a decision model. A HAB POMDP should feed the full filtered state distribution to a policy, not a Viterbi label. For genuine forecasting, transitions should use lagged or separately forecast environmental inputs. A network model should also allow residual spatial dependence, transport lags, site-specific effects, and uncertainty in modeled covariates.

### 3.6 Jackson: continuous-time multi-state modeling with `msm`

**Source.** Christopher Jackson, *Multi-state modelling with R: the msm package*, v1.8.2, 7 November 2024. [Local PDF](</Users/angelinale/Downloads/multi-state modelling (Jackson).pdf>). This attachment is a package manual, not itself a peer-reviewed article. Its peer-reviewed foundation is Jackson (2011), “Multi-State Models for Panel Data: The msm Package for R,” *Journal of Statistical Software* 38(8), 1–28 ([official record](https://doi.org/10.18637/jss.v038.i08)).

#### Why it matters here

Jackson supplies the most disciplined statistical foundation among the supplied sources for irregular panel observations. It defines continuous-time transition intensities, \(P(t)=\exp(tQ)\), exact versus panel-observed transition likelihoods, covariate-dependent intensities, hidden/misclassified states, continuous emissions, and forward likelihood evaluation (manual pp. 1–13).

The key practical distinction is that an observed move from category 1 to category 3 over a long interval does not imply a direct instantaneous \(1\to3\) jump; unobserved state 2 can have been traversed. That is why the repository’s continuous-time generator is preferable to raw observed transition counts.

The manual explicitly distinguishes fixed or state-conditionally scheduled visits from self-selected visits triggered by unobserved worsening. `msm` does not handle informative observation times (pp. 4–5). HAB sampling may be triggered by visible scum, complaints, weather, prior results, or alerts. That possibility is central, not peripheral: absence of a sample is then informative about the operational process.

The manual recommends multiple starting values, scrutiny of sparse transition pairs, convergence/Hessian checks, scaling, observed-versus-expected prevalence and transition diagnostics, time-varying intensity or semi-Markov alternatives, and bootstrap uncertainty (pp. 18 and 28–37). It also clearly distinguishes a Viterbi path from posterior state probabilities (pp. 43–44).

#### Boundaries and direct errors

`msm` estimates an uncontrolled stochastic process. It does not define rewards, solve a POMDP, or give observational covariates causal action meaning. Its exponential-sojourn assumption may also fail when duration in a HAB regime matters.

The manual itself contains errors: its matrix-exponential series on p. 6 omits the linear \(X\) term; the uniform support in Table 1 is reversed; optimizer-default descriptions conflict; and several displayed emission estimates on p. 41 fall outside their printed confidence intervals. Those issues counsel verification of software output rather than blind transcription, but they do not undermine the standard CTMC foundation.

#### HAB lesson

Use Jackson’s framework to separate latent dynamics from observation error, exploit irregular elapsed times, test homogeneity, inspect sparse transitions, consider duration dependence, and investigate informative sampling. Treat fitted \(Q\), emissions, and their uncertainty as inputs to—not substitutes for—the eventual POMDP.

### 3.7 Bashar and Torres-Machi (2023): satellite monitoring inside a POMDP

**Source.** Mohammad Bashar and Cristina Torres-Machi, “Quantifying the Value of Satellite-Based Pavement Monitoring in Partially Observable Stochastic Environments,” *Journal of Computing in Civil Engineering* 37(3) (2023), 04023004. [Local PDF](</Users/angelinale/Downloads/VoI in POMDP (Bashar & Torres).pdf>). POMDP and data definitions are PDF pp. 2–6; results are pp. 7–9; limitations are p. 10.

#### Design

The paper values a WorldView-3 brightness-based observation channel for three pavement condition states. It considers:

1. highways already receiving assumed 90%-accurate annual ground surveys, where satellite monitoring is an optional supplement; and
2. otherwise unmonitored local/ancillary roads, where satellite monitoring is compared with blind management.

The transition matrices are empirical three-state deterioration/repair matrices. Actions combine one of four maintenance choices with no inspection, satellite inspection, or annual survey. Rewards include maintenance, inspection, and user-penalty costs; the discount factor is 0.95 (Table 1, p. 5). Policies are solved for an infinite horizon and then simulated over 25 years.

The actual satellite observation matrix uses five brightness bins estimated from only 42 road sections and one image date (Fig. 4, p. 6). It is not summarized well by one diagonal accuracy because the observation alphabet has five bins and their full distributions distinguish states.

#### Results

In the already-monitored setting, the policy generally continues selecting annual surveys; actual satellite data add little. In the unmonitored setting, satellite data are selected annually and reduce simulated 25-year cost by about 6.5%: \(-349.4\) versus \(-373.8\) thousand dollars per lane-km (Table 2, p. 9). In the first setting, the apparent 0.75% difference has overlapping reported intervals.

The paper also constructs an artificial one-parameter observation family. “Accuracy” \(p\) places probability on one chosen bin and spreads error across selected neighboring bins. Satellite information becomes attractive above roughly 70% in the already-monitored setting and around 30% in the unmonitored setting (Fig. 9, pp. 7–8). These are properties of that particular channel topology, costs, transition model, baseline monitoring, and policy—not universal remote-sensing thresholds.

#### Problems and limitations

- Inspection and maintenance are selected together before the observation, so current inspection cannot change current maintenance. If the operational decision is sample first and close/harvest second, Song’s two-substep timing is needed.
- The policy is optimized for an infinite horizon but evaluated over 25 years. It is not necessarily the optimal finite-horizon policy, and the discussion inconsistently attributes action choices to a 25-year horizon.
- In the nested optional setting, exact gross VoI cannot be negative; small negative values in Fig. 10 must be simulation/solver uncertainty or objective inconsistency.
- Equation 4 repeats the state index incorrectly; the outer sum should be over next state \(s'\).
- The actual observation matrix is based on a very small, single-date sample, without held-out validation or uncertainty.
- Pixel brightness is confounded by surface material, lane markings, vehicles, shadow, illumination, repair history, and view geometry.
- Transition matrices contain structural zeros and are presented without sample sizes, smoothing, uncertainty, or adjustment for maintenance selection.
- Initial belief is arbitrarily uniform; parameter and solver uncertainty do not appear in the reported intervals.
- Network budgets, spatial correlation, shared image acquisition cost, crew constraints, and road prioritization are absent.

#### HAB lesson

The transferable pattern is the comparison of nested monitoring regimes. For HAB, the remote product’s full observation distribution should be estimated on held-out time/site data, not replaced by a symmetric accuracy slider. Gross optional information value, per-use cost, and fixed system cost should be kept distinct. The operating threshold should be chosen by the public-health/economic loss, not generic accuracy.

### 3.8 Cross-literature synthesis

The sources occupy different layers:

| Source | Primary contribution | What it does not supply |
|---|---|---|
| Arcieri 2023 | Bayesian action-conditioned HMM; posterior uncertainty; DP/QMDP connection | Causal identification, Bayes-adaptive optimality, held-out field policy evidence |
| Arcieri 2024 | Belief-input PPO; posterior domain randomization | Fair universal architecture comparison; verified deployable parameter filtering; field evaluation |
| Arcieri 2026 | Learned belief-transition/update architecture | Fully convincing normalized temporal factorization; identifiable physical states; epistemic model uncertainty |
| Song 2022 | Nonstationary augmented state; observe-update-act timing; optional information flow | A seasonally recurring ecological model; independently validated transitions |
| Kim 2022 | Covariate-driven multivariate bloom HMM | Actions, rewards, belief policy, VoI, true forward validation |
| Jackson | CTMC/HMM inference, irregular gaps, misclassification, diagnostics | A decision model or causal action effects |
| Bashar & Torres | Empirical satellite channel and nested monitoring comparison | General accuracy threshold; inspect-then-act timing; parameter/solver uncertainty |

The coherent synthesis is not to copy any one paper. It is to use Jackson’s inference discipline, Kim’s covariate/spatial motivation, Song’s timing, Bashar’s monitoring-regime comparison, and Arcieri’s insistence on belief and parameter uncertainty—while correcting the documented mathematical and validation weaknesses.

---

## 4. Repository anatomy: what exists and what does not

### 4.1 Current source layout

| Component | Role in the current snapshot |
|---|---|
| [data configuration](/Users/angelinale/voi-remsen/src/voi_remsen/config.py:1) | Names raw files, four co-located sites, two *Pseudo-nitzschia* columns, tissue threshold, and weekly rule. |
| [data pipeline](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:9) | Parses calHABMAP and CDPH files; aggregates by week; constructs cell, tissue, and joint panels. |
| [panel entry point](/Users/angelinale/voi-remsen/scripts/00_build_panels.py:3) | Runs the data builders and a coverage visualization. |
| [R `msm` model](/Users/angelinale/voi-remsen/R/fit_toxicity_msm.R:14) | Fits two- through four-state continuous-time HMMs to tissue categories alone. |
| [stationary joint HMM](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/model.py:21) | Fits a pooled birth-death CTMC with cell and tissue emissions. |
| [emissions](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/emissions.py:1) | Tobit-Normal log-cell channel and three-category tissue channel. |
| [matrix exponential](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/linalg.py:1) | Differentiable Taylor scaling-and-squaring approximation to \(\exp(Q)\). |
| [selection/evaluation](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/selection.py:1) | Reconstructs forward likelihood factors, applies smoothing, and prepares inputs for ArviZ PSIS-LOO. |
| [stationary fit script](/Users/angelinale/voi-remsen/scripts/01_fit_worldmodel.py:24) | Fits \(K=2,3,4\), writes traces/diagnostics/comparison, and exports one converged model. |
| [seasonal HMM](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/seasonal.py:1) | Makes transition rates Fourier functions of week of year and propagates through missing weeks. |
| [seasonal fit script](/Users/angelinale/voi-remsen/scripts/02_fit_seasonal.py:33) | Fits the seasonal \(K=2\) model, compares it with the stationary model, and exports an annual transition cycle. |
| [visualization](/Users/angelinale/voi-remsen/src/voi_remsen/viz.py:42) | Draws data coverage, stationary-smoothed monthly state probability, and a seasonal-cycle figure. |
| [export](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/export.py:22) | Writes posterior transition and emission arrays for an intended later decision layer. |
| [simulation helpers](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/simulate.py:19) | Generates synthetic data in the assumed model family. |
| Three Arcieri notebooks | Historical tutorials for Markov-chain, MDP/POMDP, and deterioration HMM inference; not integrated with the HAB pipeline. |

The root [README](/Users/angelinale/voi-remsen/README.md) is empty. The world-model package points to `docs/methodology.md`, which does not exist ([package docstring, lines 1–5](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/__init__.py:1)). There are no checked-in tests. The project declaration describes a world model plus POMDP decision model, but only the former exists ([`pyproject.toml`, lines 1–15](/Users/angelinale/voi-remsen/pyproject.toml:1)). `scipy` and `xarray` are directly imported in the selection module but are not declared as direct dependencies; they may arrive transitively through the optional inference stack, which is fragile dependency management.

### 4.2 The C-HARM files are present but unused

Four local NetCDF files contain v3.1 nowcast and one-, two-, and three-day forecast products. The zero-day file spans 1 November 2022 through 10 July 2026, has daily 3 km cells over coastal California and southern Oregon, and contains `pseudo_nitzschia`, `particulate_domoic`, `cellular_domoic`, and gap-filled chlorophyll variables. Across the joint field panel, the local C-HARM era contains 683 site-weeks, including 674 weeks with cells, 155 with tissue DA, and 146 with both.

No source file opens these datasets or mentions their scientific variables. The data pipeline also drops latitude and longitude when it creates weekly tables, and the four-site configuration contains names but no coordinates. The code therefore cannot yet perform even the first observation-model task: selecting a pixel or neighborhood, defining lead time, aligning daily products to decision epochs, and pairing predictions with held-out field/tissue outcomes.

### 4.3 The reference README is stale relative to the snapshot

The reference README says an optimal-stopping POMDP, point-based value iteration, theoretical sensor sweep, actual C-HARM analysis, and `hab_pomdp` package already exist, with a reported break-even near 0.63–0.73 ([lines 17–36](/Users/angelinale/voi-remsen/notebooks/reference/README_arcieri.md:17)). None of the named modules, notebooks, or write-ups is present. This may reflect an incomplete fork or removed earlier work; either way, the current source cannot support those claims. The stale material should be treated as historical intent, not documentation of the reviewed implementation.

---

## 5. Data pipeline: exact behavior and scientific implications

### 5.1 Input scale and coverage

The raw calHABMAP file has 7,203 rows, 43 columns, 17 sites, and dates from 3 June 2005 through 8 July 2026. Weekly aggregation produces 7,004 site-weeks. The joint model then restricts attention to four nominally co-located shore/tissue sites:

| Code | Location | Joint observed weeks | Time span | Cell weeks | Tissue-DA weeks | Both in same week | Recorded closure weeks |
|---|---|---:|---|---:|---:|---:|---:|
| CPP | Cal Poly Pier / San Luis Obispo | 938 | 2004-04-05 to 2026-07-06 | 821 | 228 | 111 | 3 |
| MW | Monterey Commercial Wharf | 635 | 2002-02-25 to 2026-06-22 | 603 | 76 | 44 | 6 |
| SCW | Santa Cruz Wharf | 1,023 | 1995-03-27 to 2026-07-06 | 667 | 663 | 307 | 22 |
| SIO | Scripps Pier / La Jolla | 952 | 1991-11-18 to 2026-06-08 | 905 | 143 | 96 | 3 |

The final joint panel has 3,548 site-weeks: 2,438 with cells only, 552 with tissue DA only, and 558 with both. There are 2,996 cell observations and 1,110 tissue observations. Tissue categories are 986 non-detects, 90 detections below 20 µg/g, and only 34 closure-level observations. The rarity of closure observations is central: a flexible model can fit overall data while learning closure risk very poorly.

Cell counts are extremely skewed: median roughly 2,950 cells/L, mean about 36,316, maximum about 3.25 million, and 914 exact zeros. Zero rates differ enormously by site—approximately 66% at SCW, 43% at CPP, 13% at MW, and 5% at SIO—already indicating that one common zero/cell emission may be absorbing site or protocol differences.

### 5.2 Weekly time construction

Both input streams apply

`to_period("W-SUN").start_time`

([calHABMAP line 14](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:14); [CDPH line 46](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:46)). In pandas, a period ending Sunday starts on Monday. The resulting dates in the interim data are Monday week starts, despite “W-SUN” being easy to read as Sunday-start. Both streams are internally aligned, so the current join is consistent. A future daily C-HARM match can nevertheless move observations across decision weeks unless timezone, week start/end, forecast issue time, and lead time are specified explicitly.

The calHABMAP timestamp is named `time_utc`, but the parser creates a timezone-naive datetime. CDPH sampling dates have no time of day. Near a boundary, California local time versus UTC can change the assigned day/week.

### 5.3 Construction of *Pseudo-nitzschia* total

The pipeline defines total cells as the row sum of the `delicatissima` and `seriata` size groups using `min_count=1` ([data line 15](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:15)). This means:

- if both are missing, total is missing;
- if one is observed and the other is missing, total equals the observed group—as if the missing group contributed zero.

That behavior is especially damaging here. There are 667 rows where `delicatissima` is missing but `seriata` is observed, and those rows are effectively the usable SCW cell series. Other sites generally sum both groups. The meaning of “total *Pseudo-nitzschia*” therefore differs systematically by site and era. A common emission distribution can interpret taxonomic reporting practice as ecological severity.

This is not merely generic missing-data caution. The current high/low states are identified mainly by this constructed cell variable. The very different zero rates and missing-group pattern are plausible direct causes of the broad low-state emission and apparent state occupancy.

### 5.4 Weekly aggregation choices

calHABMAP measurements within a site-week are averaged ([lines 23–28](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:23)). Tissue DA within a site-week is reduced to the maximum concentration and maximum category ([lines 61–66](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:61)). There are 100 site-weeks with multiple calHABMAP samples and 44 weeks with multiple tissue samples.

Mean water cells and maximum tissue toxin answer different within-week questions. The maximum is defensible for a conservative regulatory target, but it shifts the likelihood toward extremes and should be stated as a decision-motivated choice. It also does not ensure that the water and tissue samples occurred in the same part of the week. The outer merge simply labels them as contemporaneous when their Monday week key matches ([lines 101–109](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:101)).

### 5.5 Tissue censoring and assay-era change

The CDPH parser recognizes a `<` modifier, stores the numeric result as a detection limit, and creates three categories: censored/non-detect; detected below 20; and at least 20 µg/g ([lines 47–58](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:47)). This is reasonable preliminary bookkeeping. The HMM then uses only the category and discards the detection-limit magnitude.

In the selected records, censor limits are approximately 1.00 µg/g for 396 samples, 2.50 for 629, and 2.55 for 4. The change is strongly associated with time: about 1 through 2009 and 2.5 from 2010 onward. A reading below 1 and a reading below 2.5 are not equivalent evidence about a latent toxin concentration. Collapsing them creates an assay-era confound that a stationary categorical emission cannot distinguish from ecological change.

Although the configuration allows mussels and oysters, the selected matched records are mussel types. Species-specific accumulation is therefore not currently modeled, but the documentation should not imply an oyster comparison that the fitted data do not contain.

### 5.6 The most relevant water-toxin channel is omitted from the world model

The weekly calHABMAP table retains particulate and total water-column DA variables in its hazard mapping ([configuration lines 24–31](/Users/angelinale/voi-remsen/src/voi_remsen/config.py:24)). The cell panel selects only total cell biomass, and the joint panel joins that biomass with CDPH tissue DA ([data lines 78–109](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:78)). Water-column DA is not used in the HMM.

That omission is consequential. C-HARM has distinct particulate-DA and cellular-DA products. A water-column DA channel is a much closer calibration target than downstream shellfish tissue. The current model instead tries to force cell biomass and lagged tissue burden onto one latent axis while leaving the available intermediate toxin evidence unused.

### 5.7 Missingness and sampling are treated as ignorable

When a channel is missing, its emission factor is set to one ([emissions lines 45–53 and 70–75](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/emissions.py:45)). This is correct under a missing-at-random/ignorable observation process conditional on modeled information. It is not evidence that the assumption holds.

Sampling intensity varies across sites and eras and may rise during suspicious conditions. Pre-2008 portions of the joint sequences are almost entirely tissue-only; cell monitoring begins later. Jackson’s visit-process warning is directly relevant. If sampling depends on unmodeled risk, a week without a sample carries information that the likelihood discards. At minimum, the project needs a sampling-process audit by site, season, prior result, alert status, and protocol era.

### 5.8 Spatial mismatch is currently hidden

The raw calHABMAP file contains latitude and longitude. Weekly aggregation groups by location name/code and drops coordinates. CDPH is linked by exact sample-site names mapped to four codes ([data lines 50–57](/Users/angelinale/voi-remsen/src/voi_remsen/data.py:50)). No uncertainty is retained about the spatial relationship between a pier sample, a shellfish collection location, and a 3 km offshore C-HARM pixel.

Before fitting a remote observation model, the project must define whether it uses the nearest water pixel, a coastal neighborhood, an offshore feature, an advected footprint, or a region. Anderson et al. show that nearshore/offshore and 3 km/pier mismatch materially affects C-HARM skill. Spatial alignment is part of the observation model, not an incidental data join.

---

## 6. The fitted world models: what the code computes

### 6.1 Tissue-only R reference model

The [R script](/Users/angelinale/voi-remsen/R/fit_toxicity_msm.R:14) fits continuous-time hidden multi-state models with \(K=2,3,4\). Only adjacent latent transitions are allowed. Tissue category is a possibly misclassified categorical observation of the latent state. Time is measured in years from each site’s first tissue observation.

The fitted model is

\[
P(S(t+\Delta)=j\mid S(t)=i)=\left[\exp(\Delta Q)\right]_{ij},
\]

and

\[
P(Y_t=c\mid S_t=k)=E_{kc},\qquad c\in\{\text{non-detect},\text{detected},\text{closure}\}.
\]

The script initializes every permitted transition intensity at 0.3/year and gives progressively more toxic starting emission rows to higher states ([lines 14–23](/Users/angelinale/voi-remsen/R/fit_toxicity_msm.R:14)). Those are starting values, not monotonicity constraints.

Saved model-selection results are:

| \(K\) | log likelihood | parameters | AIC | BIC |
|---:|---:|---:|---:|---:|
| 2 | -388.3 | 6 | 788.6 | 818.7 |
| 3 | -387.0 | 10 | 793.9 | 844.1 |
| 4 | -384.7 | 14 | 797.3 | 867.5 |

Both AIC and BIC favor \(K=2\). Its fitted generator is approximately

\[
Q=
\begin{bmatrix}
-1.514&1.514\\
14.256&-14.256
\end{bmatrix}
\quad\text{per year},
\]

so mean sojourn time is about 34.5 weeks in the clean state and 3.65 weeks in the toxic-episode state. The clean state emits 97.76% non-detect, 2.22% detected, and 0.02% closure. The episode state emits 22.17% non-detect, 52.05% detected, and 25.78% closure. This is an interpretable tissue-toxicity episode model.

Important qualifications:

- **Initial state is silently fixed.** Because the call does not estimate initial probabilities, the saved fits use \([1,0,\ldots]\): every site begins in state 1. That material assumption is not documented in the script.
- **Higher-state toxicity is not constrained.** In the saved \(K=4\) fit, state 3 emits about 80% closure while state 4 emits about 98.6% merely detected and only 1.3% closure. A numeric state index is not a guaranteed severity ordering.
- **The model is pooled and uncontrolled.** It has no site effects, season, assay-era effect, sampling model, actions, or tissue kinetics.
- **AIC/BIC are in-sample.** There is no site/time holdout, calibration assessment, or uncertainty propagation into a policy.
- **The Python “cross-check” is not numerical agreement.** The Python joint model’s high-state stationary occupancy is about 57%; the R tissue-episode occupancy is about 9.6%. Their rates and state meanings are very different. The R model is useful evidence that short toxic episodes exist, not validation of the Python state.

### 6.2 Stationary Bayesian joint HMM

#### Sequence construction

The Python model creates one observed sequence for each site ([model lines 22–53](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/model.py:22)). Weeks with neither channel are omitted; the elapsed gap to the next observed week is retained. Cells are represented as \(\log_{10}(\text{cells}+1)\). Tissue categories become a three-element one-hot vector. Missing observations receive a mask.

#### Latent dynamics

For \(K\) ordered states, only adjacent moves are allowed:

\[
Q_{i,i+1}=u_i,\quad Q_{i+1,i}=d_i,\quad
Q_{ii}=-\sum_{j\ne i}Q_{ij}.
\]

The one-week transition is \(P_1=\exp(Q)\), and a \(g\)-week gap uses \(P_1^g\) ([model lines 56–77](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/model.py:56)). Gaps longer than 52 weeks are clipped to 52 ([lines 110–129 and 163](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/model.py:110)). Eleven observed gaps are clipped, with a maximum actual gap of 255 weeks. For the saved \(K=2\) posterior, the chain is already almost stationary by 52 weeks: the maximum difference between \(P^{52}\) and \(P^{255}\) over posterior draws is about \(1.4\times10^{-4}\). The approximation is currently negligible, although it should remain a checked approximation rather than a universal constant.

#### Emission model

For cells, the code assumes a latent Normal on log scale with censoring at zero:

\[
Y_t^*\mid S_t=k\sim N(\mu_k,\sigma_k^2),
\qquad
Y_t=\max(Y_t^*,0).
\]

An observed positive value uses the Normal density; an exact zero uses the left-tail probability \(\Phi(-\mu_k/\sigma_k)\) ([emissions lines 7–15 and 55–75](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/emissions.py:7)). This treats every recorded zero as censoring at exactly zero, not as a structural absence, rounding result, taxonomic non-report, or assay-specific non-detect. No known cell-count detection limit appears in the data model.

Tissue DA is a free categorical row \(E_{k,:}\). The code assumes cell and tissue emissions are independent conditional on the same state, multiplying the two likelihoods. A missing channel multiplies by one.

#### Priors and identification

Onset and recovery rates have Exponential(1) priors per week. Cell means are forced nonnegative and separated by at least 0.5 log10 units; state scales are separately estimated with a 0.15 floor; tissue rows and a common initial distribution have uniform Dirichlet priors ([model lines 110–153](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/model.py:110)). Ordering the means helps label states, but does not guarantee that the full emission distributions are stochastically ordered when their scales differ.

All four sites share the same generator, cell means/scales, tissue rows, and initial distribution. The only site distinction is that observations are separate sequences.

#### Likelihood and sampling

The forward algorithm predicts the next state through the gap transition, multiplies by that week’s emission likelihood, rescales, and sums log scale factors ([model lines 80–106](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/model.py:80)). Hidden states are analytically marginalized. This is a major improvement over the reference notebooks’ explicit sampling of thousands of categorical states.

The fit script runs four chains, 800 tuning iterations and 500 retained draws, target acceptance 0.95, with a fixed sampler seed ([fit script lines 24–31](/Users/angelinale/voi-remsen/scripts/01_fit_worldmodel.py:24)).

### 6.3 What the saved \(K=2\) fit says

The converged two-state posterior has mean weekly rates:

\[
u\approx0.1177,\qquad d\approx0.0896,
\]

with mean one-week transition matrix

\[
P_1\approx
\begin{bmatrix}
0.8938&0.1062\\
0.0809&0.9191
\end{bmatrix}.
\]

If these rates are frozen, the high-state equilibrium fraction is approximately \(u/(u+d)=56.8\%\). Mean sojourns are about 8.5 weeks in state 0 and 11.1 weeks in state 1.

Mean cell-emission parameters are:

| State | \(\mu\), log10(cells+1) | \(\sigma\) | Nominal interpretation |
|---|---:|---:|---|
| 0 | 0.008 | 3.35 | low-mean, extremely diffuse/zero-absorbing regime |
| 1 | 4.08 | 0.749 | high-cell regime centered near 12,000 cells/L |

Mean tissue-emission rows are:

| State | Non-detect | Detected <20 | Closure ≥20 |
|---|---:|---:|---:|
| 0 | 99.34% | 0.48% | 0.18% |
| 1 | 76.67% | 16.78% | 6.55% |

The high state is therefore **not** a toxic or closure state. It is a high-cell regime in which tissue is still non-detect in more than three quarters of observations. Its closure probability is higher than state 0’s, but closure is not its defining feature.

Smoothed mean high-state occupancy also differs dramatically by site: about 43% at CPP, 74% at MW, 29% at SCW, and 90% at SIO. Given the site-specific measurement construction and zero rates, these numbers cannot safely be interpreted as comparable ecological bloom prevalence.

### 6.4 Critical cell-emission pathology

At the posterior-mean emission parameters, the state-0 scale of 3.35 log10 units is not merely “wide.” It makes the supposed low state assign substantial probability to biologically extreme values:

- about 3.68% probability above one million cells/L;
- about 0.36% above one billion cells/L; and
- an implied expected raw count around \(8.3\times10^{12}\) cells/L under the censored-lognormal transformation, vastly beyond the observed maximum of \(3.25\times10^6\).

Because state 1 is much narrower, the cell likelihood ratio is non-monotone. State 1 becomes more likely than state 0 near \(y=2.66\) (roughly 450 cells/L), but state 0 becomes more likely again above \(y=5.93\) (roughly 857,000 cells/L). At the observed maximum \(y\approx6.51\), the direct state-1 cell likelihood is only about 15% of the state-0 likelihood.

This appears in smoothed outputs: mean high-state probability is about 96% for observations between \(10^5\) and \(10^6\) cells/L, then falls to about 73% above \(10^6\), with individual extreme weeks near 25–34%. A “severity” model that can classify a larger count as less severe is unsafe for closure decisions.

The likely mechanism is clear. A nonnegative low-state mean plus a Tobit censoring model must explain many exact zeros, positive low values, site/protocol heterogeneity, and outliers. It does so by inflating \(\sigma_0\). Ordering only \(\mu\) does not impose an ordered likelihood ratio or stochastic severity. Candidate remedies to investigate include a separate zero/hurdle process, explicit detection limits, site/protocol effects, robust or mixture emissions, lagged multi-compartment states, and constraints that make hazard ordering valid over the observed range. This is a model-design requirement, not a request to merely increase sampling iterations.

### 6.5 Sampling diagnostics and state count

Saved diagnostics are:

| \(K\) | reported ELPD | \(p_{\text{LOO}}\) | max Pareto \(k\) | max \(\hat R\) | divergences |
|---:|---:|---:|---:|---:|---:|
| 2 | -4567.1 | 7.3 | 0.39 | 1.013 | 0 |
| 3 | -4420.7 | 121.7 | 0.50 | 1.544 | 0 |
| 4 | -4295.3 | 69.8 | 0.49 | 1.536 | 0 |

The \(K=2\) run is computationally healthy: no divergences or maximum-depth hits and effective sample sizes comfortably above several hundred. The \(K=3\) and \(K=4\) chains are severely unmixed. A longer six-chain \(K=3\) attempt has maximum \(\hat R\) about 1.623, so more runtime did not resolve the structural multimodality.

The fit script correctly refuses to export nonconverged fits ([lines 73–85](/Users/angelinale/voi-remsen/scripts/01_fit_worldmodel.py:73)). The correct scientific conclusion is:

> \(K=2\) is the only current parameterization that was successfully sampled.

It is not:

> \(K=2\) has been validated as the true or predictively optimal number of HAB states.

Pohle et al. provide a further warning: even if \(K=3\) or \(K=4\) eventually converges and scores better, extra states may absorb omitted site, season, lag, assay, tail, or dependence structure rather than represent ecological regimes.

### 6.6 Why the reported PSIS-LOO is not prospective validation

The fitted likelihood is stored as one PyMC potential, so the selection module reconstructs the forward scale factors

\[
\ell_t(\theta)=\log P(y_t\mid y_{<t},\theta)
\]

and passes them to ArviZ as pointwise log likelihoods ([selection lines 1–17 and 218–255](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/selection.py:1)). Their sum is exactly the full-sequence log likelihood. That algebraic decomposition is correct.

The leave-one-out interpretation is not. If observation \(y_t\) is held out, the filter at \(t\) must skip its update, and every later prediction must be recomputed from the altered belief. The stored factor \(\ell_{t+1}\) still conditions on a filtering distribution that incorporated \(y_t\); so do all later factors. ArviZ removes or reweights one factor, not the observation’s downstream effect. Moreover, the posterior draws were estimated using later data.

Bürkner, Gabry, and Vehtari define the relevant forecasting target as

\[
P(y_{i+1:i+M}\mid y_{1:i}),
\]

with both the predictive distribution and parameter posterior conditioned only on past data (2020, §2, pp. 3–4). They show ordinary LOO can be positively biased when future observations inform predictions of the past (§3.1 and §4.1). See the [open manuscript](https://arxiv.org/pdf/1902.06281) and [journal DOI](https://doi.org/10.1080/00949655.2020.1783262).

The current calculation is useful as a forward-factor influence diagnostic. A small Pareto \(k\) means importance weights for that calculation are stable; it does **not** prove that the target is correct. The saved ELPD values and weights must not be cited as future forecast performance.

A valid evaluation program needs at least:

1. rolling-origin/leave-future-out forecasts within site, with horizon matched to the management decision;
2. blocked seasonal/year holdouts that do not train on later years;
3. leave-one-site-out or hierarchical held-out-site validation if geographic transfer matters;
4. calibration and proper scores separately for cells, water DA, tissue DA, and closure events; and
5. downstream policy-value stability, because a small likelihood improvement may not change a decision while a rare-event calibration error can change it greatly.

### 6.7 Misleading and fragile selection artifacts

The saved [`worldmodel_compare.csv`](/Users/angelinale/voi-remsen/data/output/worldmodel_compare.csv:1) ranks \(K=4\) first and assigns it weight 1.0 even though its maximum \(\hat R\) is 1.536. The script calls `az.compare` on every fit before applying its convergence filter ([fit script lines 68–83](/Users/angelinale/voi-remsen/scripts/01_fit_worldmodel.py:68)). A reader opening the comparison file without the console log gets a scientifically invalid ranking.

The convergence gate itself uses \(p_{\text{LOO}}<3\max(K)^2\), an ad hoc threshold that changes with the set of candidates supplied on the command line. It filters rounded diagnostic values and does not include divergence or Pareto-\(k\) criteria. A stale `worldmodel_modelsel.csv` reports still different diagnostics and is not produced by the current script. Artifact names do not distinguish current, stale, converged, or exploratory outputs.

The selection module also creates a log-likelihood variable with the same name as its observation dimension. This is nonstandard xarray structure: specifying `var_name="obs"` works, while automatic `az.loo(idata)` fails. No site/date coordinates are attached, so an influential Pareto point cannot be traced back to a biological observation. These are reproducibility/interface defects even apart from the validation estimand.

### 6.8 Missing posterior predictive checks

No checked-in workflow asks whether replicated data resemble the observations. Essential checks are absent for:

- zero proportions by site and era;
- high cell-count tails and monotonic likelihood ordering;
- tissue-category counts, especially the 34 closure weeks;
- paired and lagged cell–water-DA–tissue relationships;
- state dwell times and transition episode length;
- seasonal residuals;
- sampling/missingness patterns; and
- site-specific predictive calibration.

The state-0 astronomical tail would have been exposed immediately by simulating raw cell counts or plotting state-specific densities over the observed range. Convergence diagnostics tell whether the sampler explored the specified posterior; they do not tell whether the specified model is sensible.

### 6.9 Numerical implementation: what is sound and what is fragile

The custom scaling-and-squaring Taylor approximation to \(\exp(Q)\) is accurate over every saved stationary and seasonal posterior draw: maximum absolute disagreement with SciPy is about \(2.3\times10^{-13}\), with row-sum error about \(2.4\times10^{-13}\). This is a strong implementation result for the fitted range.

The forward recursion works in ordinary probability space between rescaling steps and has no floor or log-sum-exp safeguard. Extreme proposed cell parameters can make all state emissions underflow to zero, producing \(\log 0\) and division by zero ([model lines 90–106](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/model.py:90)). The fixed matrix-exponential approximation is likewise not guarded against invalid probabilities far into the prior/proposal tail. No such failure is evident in the saved \(K=2\) posterior, but a robust fitting engine should defend the whole sampled domain.

### 6.10 Seasonal extension

#### Mathematical construction

The seasonal model replaces each constant adjacent rate with a Fourier function of week of year:

\[
\log u_i(w)=\beta^u_{i0}+
\sum_{h=1}^H
\left[\beta^u_{ih,c}\cos\left(\frac{2\pi hw}{52.18}\right)+
\beta^u_{ih,s}\sin\left(\frac{2\pi hw}{52.18}\right)\right],
\]

with an analogous formula for recovery \(d_i(w)\) ([seasonal lines 26–37 and 122–152](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/seasonal.py:26)). Because \(Q\) changes each week, a gap transition is an ordered product, not a power of one matrix. The code correctly expands each site to its full weekly grid, supplies an emission factor of one on unobserved weeks, and applies the appropriate weekly transition in sequence ([seasonal lines 46–85 and 102–119](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/seasonal.py:46)).

The saved \(K=2,H=1\) fit is computationally converged: maximum \(\hat R\) 1.011, no divergences, and good effective sample sizes. Posterior-mean onset rates range roughly 0.054–0.243/week; recovery rates are much flatter, roughly 0.081–0.099/week. The fitted dynamics therefore attribute most seasonality to entering the high-cell state.

#### Frozen equilibrium is not propagated seasonal prevalence

The export defines

\[
\operatorname{occ}_{\text{frozen}}(w)=\frac{u(w)}{u(w)+d(w)}
\]

([seasonal fit lines 66–80](/Users/angelinale/voi-remsen/scripts/02_fit_seasonal.py:66)). This is the stationary high-state fraction **if week \(w\)’s rates were frozen forever**. In a time-varying chain, prevalence retains memory: it must be propagated through the ordered sequence of transition matrices rather than read from one week’s rate ratio.

For an idealized cycle formed by repeatedly applying the 53 exported week-bin matrices:

| Quantity | Minimum | Maximum | Peak week |
|---|---:|---:|---:|
| Exported frozen-rate ratio, transformed from mean coefficients | 35.6% | 74.8% | 17 |
| Propagated 53-bin plug-in high-state probability | 39.6% | 73.6% | 20 |

The maximum pointwise difference is about 10.4 percentage points and the mean absolute difference about 5.6 points. The propagated plug-in high-state probability peaks around three weeks later than the instantaneous equilibrium target. The plotted curve is a posterior-draw average of frozen-rate ratios rather than the exact exported mean-coefficient transform; the two differ by at most about 0.15 percentage points in this fit and both peak at week 17.

The visualization labels that frozen curve `P(in bloom state)` and derives an “onset peak” annotation from its maximum ([visualization lines 183–245](/Users/angelinale/voi-remsen/src/voi_remsen/viz.py:183)). The probability label is incorrect because the curve is not propagated prevalence. The annotation is also not justified by the quantity being maximized, although it happens to identify week 17 here, which is also the fitted onset-rate maximum.

The 53-bin propagated diagnostic should not itself be called the exact fitted periodic occupancy without qualification. The implementation maps real dates to week-of-year bins while using a 52.1786-week harmonic period; leap years and calendar alignment can repeat or skip bins. Exact deployment beliefs are calendar- and history-conditioned. The propagated 53-bin cycle is nevertheless the appropriate like-for-like diagnostic showing why the frozen ratio is not an occupancy curve.

#### Empirical-Bayes reuse of the same data

The seasonal intercept priors are hard-coded at \(\log(0.12)\) and \(\log(0.09)\), essentially the stationary posterior mean rates from the same observations ([seasonal line 137](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/seasonal.py:137)). The seasonal model then refits those observations as though the anchors were external fixed prior information. This two-stage empirical-Bayes procedure understates uncertainty and makes the stationary-versus-seasonal comparison non-clean. A principled sequential analysis would document the prior construction and propagate first-stage uncertainty; a clean model comparison would use common external priors and held-out data.

#### Circular “cross-check”

The seasonal figure compares the seasonal model with monthly smoothed state probabilities from the stationary HMM and calls the latter “empirical” and “nonparametric” ([visualization lines 183–230](/Users/angelinale/voi-remsen/src/voi_remsen/viz.py:183)). Both curves use the same observations, the same two emission families, and closely related latent states. Agreement is in-sample model agreement, not independent ecological validation.

The stationary monthly smoother yields approximate high-state probabilities of 39% in January, 50% February, 63% March, 68% April, 74% May, 67% June, 61% July, 61% August, 50% September, 48% October, 50% November, and 47% December. These values are also observation-schedule weighted: months/sites with more recorded weeks contribute more. Their narrow 5–95% bars propagate parameter draws through a smoother but not latent-state realization, site/year sampling variability, model-form uncertainty, or missingness uncertainty.

The same stationary visualization contains a separate confirmed labeling error. It calculates `overall` as the observation-weighted mean of smoothed posterior state probabilities over recorded site-weeks, then labels that line “= stationary occupancy” ([visualization lines 139–142](/Users/angelinale/voi-remsen/src/voi_remsen/viz.py:139)). That empirical posterior average is not generally the CTMC equilibrium \(u/(u+d)\): it depends on finite-sequence initialization, unequal site and observation weights, missingness schedules, and the observed data.

#### Additional seasonal-interface issues

- One shared free \(\pi_0\) initializes four sequences that begin in different calendar weeks and years. A periodic model needs week-conditioned starting distributions or sequence-specific initial beliefs.
- The script saves its NetCDF before attaching reconstructed likelihood factors, so the saved seasonal trace cannot reproduce the comparison without rerunning custom reconstruction ([fit lines 41–56](/Users/angelinale/voi-remsen/scripts/02_fit_seasonal.py:41)).
- Unlike the stationary script, it prints convergence but exports regardless of failure. The current run converged; future failures would still produce an authoritative file.
- `P_woy_mean` is computed by transforming mean coefficients, not by averaging posterior transition matrices. The difference is small in the present fit (maximum about 0.0019) but the name is mathematically inaccurate.
- Only one harmonic is fitted. There is no sensitivity to harmonic order, shrinkage by order, physical covariates, year effects, or site-specific seasonality.

### 6.11 Export boundary and decision-layer incompleteness

The stationary export includes \(K\), mean/draws of \(Q\) and weekly \(P\), cell means/scales, tissue rows, and rate units ([export lines 22–38](/Users/angelinale/voi-remsen/src/voi_remsen/worldmodel/export.py:22)). The seasonal export includes harmonic coefficient draws, a plug-in weekly transition cycle, emissions, and the incorrect frozen occupancy curve.

Neither export includes:

- the fitted initial distribution \(\pi_0\) or a current filtered belief;
- state names and defensible decision semantics;
- tissue category and cell-transform definitions;
- site/protocol structure;
- training dates and data hashes;
- code/package versions;
- convergence and validation diagnostics; or
- model/version provenance.

The two schemas differ, and no consumer selects between them. `worldmodel.npz` remains the stationary model even though a separate seasonal comparison is presented. Most importantly, there is no consumer at all: no decision-layer module loads either export.

---

## 7. Deep audit of the Arcieri reference notebooks

### 7.1 How these notebooks should be understood

The three notebooks contain 171 cells and are historical tutorials, not the executed code behind all three Arcieri papers. The [reference README, lines 1–7](/Users/angelinale/voi-remsen/notebooks/reference/README_arcieri.md:1) accurately says the deterioration notebook is merely “close” to the 2023 paper. Key differences are larger than that phrase may suggest:

- the notebooks use simulated data;
- central examples use three states and two actions, not four states and three actions;
- they use Gaussian or truncated-Normal emissions, not the paper’s truncated Student-\(t\);
- they use old PyMC3/Theano APIs and explicitly sample discrete trajectories;
- they contain no reward, belief-policy solver, QMDP, posterior Q-value averaging, PPO, domain randomization, or DBMM; and
- the current HAB model does not import or execute notebook code.

The notebook metadata targets Python 3.8.8 for the first two and 3.7.7 for the third. Saved paths in the third reference a Python 3.8 environment, already showing statefulness. The current project targets Python 3.11 or later and modern PyMC. These notebooks are not runnable in the current environment without a deliberate archival environment and API migration.

Notebook cell numbers below are **zero-based**, matching the JSON cell order.

### 7.2 `HMMs for Markov Chains.ipynb`

**Source:** [local notebook](</Users/angelinale/voi-remsen/notebooks/reference/HMMs for Markov Chains.ipynb>).

#### What it teaches

Cells 4–11 generate a three-state Markov chain and state-dependent autoregressive heteroskedastic Gaussian measurements. Cells 14–20 infer a transition matrix when states are observed. Cells 23–31 jointly infer latent states, transitions, state means/scales, and autoregression. Cells 34–49 extend the explicit-state model to multiple sequences. Cells 50–67 attempt a marginalized forward algorithm, motivated by the inefficiency of sampling categorical trajectories.

The conceptual progression is sensible. The current repository’s analytic forward model is the right descendant of the last section, not the explicit-state sampler.

#### Confirmed defects

1. **Ground-truth equilibrium leakage in cell 17.** The inferred transition random variable is `transition_mat`, but the equilibrium calculation is called with the global generating matrix `p_transition`. The first fully observed fit therefore uses a true initial distribution unavailable in real inference. Later cells pass the inferred matrix correctly.

2. **The marginalized forward likelihood is off by one in cells 55 and 63.** It initializes the forward message with the emission likelihood of \(y_0\), then scans the sequence starting at \(y_0\) for \(T-1\) steps. Consequently \(y_0\) is counted twice and \(y_{T-1}\) is omitted. The multiple-sequence version repeats the error for each sequence. The likelihood also omits the intended equilibrium/initial distribution, effectively using equal unnormalized starting weights.

3. **A markdown distribution is wrong.** Cell 24 describes a state conditional as Dirichlet. A transition row has a Dirichlet prior; the state drawn from that row is categorical.

4. **Posterior mean labels are rounded in cells 30 and 48.** Numeric category labels are not continuous quantities. Rounding their posterior mean is not a posterior mode and can turn a bimodal state posterior into a spurious middle label. It also ignores label-switching uncertainty.

5. **Exact row-sum assertions are fragile.** Cell 5 compares floating-point sums to exactly one. The current constants happen to pass; tolerance-based validation is safer for general inputs.

#### Saved diagnostics contradict narrative claims

The explicit latent single-sequence fit took 1,622 seconds and reports \(\hat R>1.2\) and effective sample size below 200. The multiple-sequence latent fit took 4,286 seconds and reports \(\hat R>1.05\). Cell 31 nevertheless says the HMM “correctly inferred” all parameters and hidden states, and cell 49 makes another broad success claim.

The marginalized single and multiple fits took 2,215 and 4,165 seconds and both report \(\hat R>1.4\), explicitly “did not converge.” Cells 59 and 67 do acknowledge those failures and possible adaptation mistakes. The off-by-one likelihood gives a concrete reason not to use their behavior as evidence against marginalization.

### 7.3 `HMMs for MDP and POMDP.ipynb`

**Source:** [local notebook](</Users/angelinale/voi-remsen/notebooks/reference/HMMs for MDP and POMDP.ipynb>).

#### What it actually does

Cells 4–10 simulate a three-state process with two action-conditioned transition matrices and autoregressive Gaussian observations. Actions are drawn independently from a fixed 90%/10% distribution, not from a state- or observation-dependent policy. Cells 13–22 infer action-specific transition matrices with observed states. Cells 25–34 infer latent states plus emissions. Cells 37–54 repeat both cases with two sequences.

This is a controlled/action-conditioned HMM identification tutorial. It is called an MDP/POMDP because transition probabilities depend on an action and states may be hidden. It does not implement decision making: there is no reward, belief update for a controller, policy, Bellman equation, exploration strategy, or VoI.

#### Confirmed defects and limitations

1. **The multiple-sequence generator ignores its argument.** Cell 37 accepts `p_transition` but calls the simulator using global `p_transitions`. A caller-supplied model is silently ignored.

2. **The markdown repeats the Dirichlet/categorical error.** Cell 26 describes the next state as Dirichlet-distributed.

3. **State estimates are again rounded label means.** Cells 34 and 53 do not provide valid Bayesian decoding.

4. **Rare action support is weak.** The second action occurs only about 10% of the time. Its transition matrix contains generating zeros, but ordinary Dirichlet rows have only interior support and cannot represent exact zeros. No support/identifiability analysis is performed.

5. **Success claims overstate the saved chains.** The single latent fit took 1,808 seconds and the multiple latent fit 5,154 seconds; both report \(\hat R>1.05\) and low effective sample sizes. Cells 31 and 54 still call the estimates satisfactory/correct.

The randomized action design does avoid the behavioral confounding that affects real maintenance logs, but only in simulation. It does not show that an action effect can be recovered causally from operator-chosen actions.

### 7.4 `HMMs for deterioration process Truncated Normal Process.ipynb`

**Source:** [local notebook](</Users/angelinale/voi-remsen/notebooks/reference/HMMs for deterioration process Truncated Normal Process.ipynb>).

This notebook has 46 code cells and essentially no explanatory markdown. Cells 1–7 define three states, two actions, a threshold action rule, action-conditioned transitions, truncated-Normal deterioration/repair emissions, and a simulated sequence. Cells 8–21 fit observed- and hidden-state versions. Cells 22–35 generate and fit two sequences. Cells 36–44 fit 20 shorter sequences.

#### Current visible code and saved outputs are out of order

The visible observation class in cell 13 has execution count 66, after the saved fits and figures in cells 14–35. The visible multiple-sequence generator in cell 22 has execution count 63, also after cells 23–35. State plots in cells 20–21 were executed before the latest fit in cell 14. Only the final 20-sequence fit at execution 70 definitely follows the current visible class definition.

This is a classic stateful-notebook reproducibility failure: the document’s displayed results cannot be assumed to come from the displayed source. A clean restart-and-run-all record is required before any numerical result is interpreted.

#### Deterioration mean is ignored

Cell 3 declares `deterioration_process(prev_em, mu, sigma)` but uses `prev_em` as the Normal mean and never uses `mu`. The fitted `ObservationModel` in cell 13 computes `mu_det` and likewise never uses it. Under the no-action channel, latent state affects scale but not mean deterioration drift.

This materially differs from the 2023 paper’s Eq. 11, where the state-specific \(\mu_{d\mid s}\) is the mean change. The notebook’s state-specific “deterioration means” cannot be identified from the deterioration likelihood because they are absent from it.

#### Initial emission generator and likelihood disagree

Cell 5 generates the initial emission by passing `prev_em=mus[state]` into the deterioration process. The generator therefore samples a Normal centered at \(\mu_s\) truncated above at \(\mu_s\). Cell 13 fits the boundary as a Normal centered at \(\mu_s\) truncated above at zero. Since \(\mu_s<0\), these are different distributions. Parameter recovery is assessed under a misspecified boundary likelihood.

#### Other defects and interpretive problems

- `p_actions` is accepted but unused because random action selection is commented out. Actions follow a deterministic threshold rule.
- The multiple generator accepts `p_transition` but uses global `p_transitions`.
- `raise("No action defined")` tries to raise a string and would itself produce a `TypeError`.
- Strong priors are centered near generating means, scales, and transitions, making simulated recovery easier and helping fix labels.
- The threshold policy means some action-state rows have little or no support; transition posteriors can be prior-dominated.
- Cells 21 and 35 report state “accuracy” from rounded posterior mean labels.
- Saved latent fits report \(\hat R>1.05\) and low effective sample sizes after 2,376 seconds (one sequence), 6,204 seconds (two sequences), and 14,503 seconds (20 short sequences). In the two-sequence fit, one scale has \(\hat R\approx1.12\) and bulk ESS about 26.

### 7.5 Implications for the current project

The notebooks are useful as historical demonstrations of:

- a transition tensor indexed by action;
- shared inference across sequences;
- state-dependent observations;
- the difficulty of discrete latent-state sampling; and
- the need for informative constraints in weakly identified HMMs.

They should not be treated as validated templates. The current repository has already improved on them by marginalizing hidden states and handling irregular gaps. It must retain that improvement while avoiding their stateful execution, global-variable leakage, invalid decoding, convergence overclaims, and likelihood errors.

The notebooks also do not justify a HAB POMDP. They contain no C-HARM observation channel, no management reward, no belief-space policy, and no VoI experiment. Their action-transition code is a building block only.

---

## 8. Crosswalk: literature claims versus the current HAB implementation

### 8.1 What the current model inherits successfully

The repository implements several sound ideas from the supplied sources:

- From Jackson, it uses a continuous-time generator and matrix exponentiation to handle irregular gaps rather than equating an observed gap with one fixed transition.
- From HMM practice, it analytically marginalizes latent states through a scaled forward filter.
- From Arcieri, it keeps posterior parameter draws rather than exporting only a point estimate.
- From ecological nonstationarity concerns, the seasonal extension uses time-varying transition rates and the correct ordered product across missing weeks.
- From latent-state identification practice, it orders cell means to reduce label switching.

These are meaningful foundations. The review’s adverse conclusions do not mean that every component is wrong; they mean the components do not yet form, and are not yet validated as, the intended decision system.

### 8.2 Where the current model departs from Arcieri

| Feature | Arcieri railway model | Current HAB model | Consequence |
|---|---|---|---|
| Time | equally spaced half-year steps | irregular weekly CTMC | Current handling is more appropriate for irregular sampling. |
| Actions | three action-conditioned transition matrices | no action | Current model cannot simulate a management intervention. |
| Observations | one autoregressive heavy-tailed signal | contemporaneous cell Tobit + tissue category | Current model lacks biological lag and uses a restrictive one-axis link. |
| State sampling | paper claims joint NUTS sampling; historical tutorials use hybrid NUTS/Gibbs; executed paper code is unavailable | analytically marginalized | Current approach is computationally preferable. |
| Parameter uncertainty | posterior retained | posterior retained | Good foundation, but export lacks initial belief/provenance. |
| Decision layer | DP/QMDP and later PPO | absent | No policy or VoI can be calculated. |

It would therefore be incorrect to describe the current world model as “Arcieri’s POMDP applied to HAB.” It is a newly designed uncontrolled CT-HMM inspired by one inference component.

### 8.3 Biomass state versus hazard state

The R and Python models provide an informative contrast:

- The tissue-only R model identifies a short, rare toxicity episode: equilibrium occupancy about 9.6%, with 25.8% closure probability inside that episode.
- The joint Python model identifies a long, common high-cell regime: frozen equilibrium about 56.8%, with only 6.5% closure probability and 76.7% non-detect inside that regime.

These are not estimates of the same latent state. The Python state is primarily organized by biomass; the R state by tissue category. Their disagreement does not prove one model is wrong, but it proves that calling the Python state simply “bloom,” “toxic,” or “unsafe” erases a crucial semantic choice.

For a decision model, the state must be defined by what predicts consequences and action value. If the management action is shellfish closure, a state axis driven mostly by cell count may be the wrong sufficient state. If the action is intensified sampling, high biomass may be a useful early-warning state even when current tissue is negative. Those two policies require different rewards and possibly different state representations.

### 8.4 Tissue burden has memory

California evidence reinforces the lag concern. Lane et al. followed approximately weekly monitoring in Monterey Bay for 17 months and found that SPATT detected domoic acid three and seven weeks before traditional bloom recognition and seven and eight weeks before shellfish toxicity in two events (2010, Abstract and pp. 655–656). They discuss temperature- and mussel-size-dependent uptake/depuration and the imperfect coupling between cell abundance, water toxin, and tissue burden. See the [DOI](https://doi.org/10.4319/lom.2010.8.0645) and [open reproduced article](https://escholarship.org/content/qt3hk1d1sf/qt3hk1d1sf_noSplash_30a5d9ae7d63c41fbda65407135fb4e5.pdf).

That result should not be turned into a universal seven- or eight-week lag; it is based on two events at one site. It establishes something more modest and important: material lag and uncoupling can occur, so a same-week common-state emission model needs explicit empirical justification and distributed-lag sensitivity analysis.

### 8.5 A scalar C-HARM “accuracy” is the wrong observation abstraction

C-HARM produces multiple continuous probabilities with different target events. An observation model for a POMDP needs the distribution of those outputs given the decision-relevant latent state or future outcome. The forecast probability is not automatically the likelihood itself, and its three products are not conditionally independent copies of one sensor.

Anderson et al. explicitly report that no single skill score captures forecast usefulness, and that sensitivity, false alarms, bias, ROC/AUC, location, threshold, and decision loss matter differently (§2.6, pp. 5–6). The early product had site/product-specific behavior and high false-positive rates in some settings (§3.2, p. 9). A future analysis should model the full score distribution or calibrated probability, stratified by product, lead, site/region, season, and prevalence.

Sweeping a symmetric matrix with one “accuracy” parameter can still be a pedagogical sensitivity experiment. Its break-even value is a property of that synthetic matrix. It must never be reported as “C-HARM must be at least 70% accurate” without specifying exactly how errors are distributed and which decision/costs generated the threshold.

### 8.6 State count is not state meaning

Pohle et al. show that HMM selection criteria often add states to absorb outliers, unmodeled temporal variation, individual heterogeneity, duration dependence, or conditional-dependence violations even when the true process has fewer genuine states (2017, §§3–4). See the [open manuscript](https://arxiv.org/pdf/1701.08673) and [DOI](https://doi.org/10.1007/s13253-017-0283-8).

This repository contains almost every cited pressure: season, site heterogeneity, protocol change, lag, a broad zero process, and possible channel dependence. A successfully fitted \(K>2\) model might be a better predictive approximation without revealing new ecological states. Conversely, forcing \(K=2\) can make one state pathologically broad. State count should be chosen jointly with emission/process adequacy and the study aim, not as an isolated contest.

---

## 9. Consolidated finding register

### 9.1 Critical and high-severity findings

| ID | Finding | Type | Confidence | Why it matters |
|---|---|---|---:|---|
| C-1 | No C-HARM matchup, sensor likelihood, action, reward, policy, or VoI exists. | Confirmed scope gap | Very high | The intended research question is not yet computed. |
| C-2 | Forward-factor PSIS-LOO is treated as observation/future predictive validation. | Confirmed evaluation defect | Very high | Model/state/seasonality rankings do not target deployment forecasting. |
| C-3 | The low-state cell emission has an astronomical tail and reverses state likelihood above ~857,000 cells/L. | Confirmed model pathology | Very high | Extreme counts can decrease inferred “bloom” probability. |
| C-4 | Nonconverged \(K=3,4\) fits appear in a table ranking \(K=4\) first with weight 1.0. | Confirmed artifact defect | Very high | Results can be misread as evidence for an invalid fit. |
| H-1 | SCW cell “total” is seriata-only while other sites usually sum two groups. | Confirmed data incompatibility | Very high | Site/reporting practice is confounded with latent ecology. |
| H-2 | Same-week cells and tissue DA are modeled as contemporaneous independent emissions. | Strong modeling assumption | High | Tissue accumulation/depuration and within-week lag are ignored. |
| H-3 | A single latent axis is used for biomass and closure risk. | Strong modeling assumption | High | The high state is 76.7% tissue non-detect and is not an unsafe state. |
| H-4 | Variable tissue LODs and DA magnitudes are collapsed into three categories. | Confirmed information loss | High | Assay-era shifts and censoring evidence are lost. |
| H-5 | One transition/emission system pools four sites and decades. | Strong modeling assumption | High | Spatial/protocol heterogeneity can masquerade as states/rates. |
| H-6 | Water-column DA—the closest local target for two C-HARM products—is omitted. | Confirmed design gap | High | The model skips the mechanistic bridge from remote forecast to tissue risk. |
| H-7 | Seasonal \(u/(u+d)\) is labeled actual occupancy and used to derive an onset annotation. | Confirmed interpretation defect | Very high | Relative to a 53-bin propagated plug-in cycle, the saved prevalence phase is about three weeks early and differs by up to 10.4 points; exact calendar occupancy remains history-dependent. |
| H-8 | Seasonal priors reuse stationary estimates from the same data as fixed anchors. | Confirmed empirical-Bayes issue | High | Uncertainty is understated and comparison is not clean. |
| H-9 | “Empirical/nonparametric” seasonal cross-check is another HMM fit to the same data. | Unsupported validation claim | Very high | Agreement is circular rather than independent evidence. |
| H-10 | Observation-weighted mean smoothed probability is labeled “stationary occupancy.” | Confirmed visualization defect | Very high | An empirical posterior average is confused with the CTMC equilibrium distribution. |
| H-11 | Exports omit initial belief and core provenance/semantics. | Confirmed interface defect | Very high | A future controller cannot faithfully reconstruct the fitted filtering problem. |
| H-12 | Monitoring/sampling occurrence is assumed ignorable. | Unchecked assumption | Medium-high | Event-triggered surveillance can bias inferred dynamics/emissions. |
| H-13 | No posterior predictive checking exists. | Confirmed validation gap | Very high | Severe tail and site problems escaped detection. |

### 9.2 Medium and reproducibility findings

| ID | Finding | Type | Confidence | Implication |
|---|---|---|---:|---|
| M-1 | `W-SUN` produces Monday starts and timestamps are timezone-naive. | Confirmed timing/documentation mismatch | High | Future daily satellite alignment can shift observations. |
| M-2 | Weekly water means are paired with weekly tissue maxima. | Modeling choice | High | The joint likelihood mixes different summaries and sample times. |
| M-3 | R model silently fixes every site’s initial state to clean. | Confirmed assumption | High | Early observations and rates can be biased. |
| M-4 | R higher-state emissions are not constrained monotone. | Confirmed identification issue | High | Numeric state labels do not guarantee toxicity ordering. |
| M-5 | Stationary selection uses an arbitrary candidate-dependent \(p_{\text{LOO}}\) gate. | Confirmed methodological weakness | High | Passing/exclusion is not based on a standard criterion. |
| M-6 | Log-likelihood artifacts have fragile xarray structure and no site/date coordinates. | Confirmed interface defect | High | Automatic evaluation fails and influential points cannot be traced. |
| M-7 | Seasonal trace is saved before likelihood reconstruction and always exported. | Confirmed reproducibility/control defect | High | Reported comparison is not self-contained; failed future fits could ship. |
| M-8 | Shared seasonal \(\pi_0\) ignores different starting calendar weeks. | Model misspecification | High | Initial beliefs are inconsistent with a periodic process. |
| M-9 | Forward likelihood has probability-space underflow paths. | Latent numerical risk | High | Extreme proposals can produce NaNs/invalid log likelihood. |
| M-10 | Root README is empty; referenced methodology/tests are absent. | Confirmed documentation gap | Very high | Scope, assumptions, and artifact provenance are opaque. |
| M-11 | Dependencies are unpinned; direct SciPy/xarray imports are undeclared. | Confirmed environment gap | High | Reproduction depends on incidental transitive packages. |
| M-12 | Saved artifacts include stale and development outputs without a manifest. | Confirmed provenance gap | High | It is unclear which files are authoritative. |

### 9.3 Arcieri notebook findings

| ID | Notebook/cell | Finding | Severity |
|---|---|---|---:|
| N-1 | Markov cells 55, 63 | First observation counted twice; last omitted; no initial distribution in marginalized likelihood. | High |
| N-2 | Markov cell 17 | Ground-truth transition matrix used for initial equilibrium instead of inferred matrix. | Medium-high |
| N-3 | All latent examples | Saved \(\hat R\)/ESS warnings conflict with broad “correctly inferred” narrative. | Medium-high |
| N-4 | MDP cell 37; deterioration cell 22 | Function argument ignored in favor of global transition tensor. | Medium |
| N-5 | MDP notebook overall | Random action-conditioned identification is labeled POMDP but contains no decision layer. | High for interpretation |
| N-6 | Deterioration cells 3, 13 | State-specific deterioration mean is accepted/computed but unused. | Medium-high |
| N-7 | Deterioration cells 3, 5, 13 | Initial-emission generator and fitted boundary likelihood differ. | High |
| N-8 | Deterioration execution order | Displayed results precede current visible class/generator definitions. | High |
| N-9 | State plots/accuracy | Rounded posterior means of category labels used as decoded state. | Medium |
| N-10 | Deterioration cell 3 | Attempts to raise a string; unused action-probability argument. | Medium/low |

### 9.4 Issues in the supplied publications that should not be copied

| Source | Issue | Status |
|---|---|---|
| Arcieri 2023 | Algorithm 1 omits prior belief in prediction. | Confirmed published pseudocode error; implementation unknown. |
| Arcieri 2023 | Averaged model-specific optimal Q-functions called robust optimal. | Conceptual overclaim; generally clairvoyant/Bayes-inadaptive. |
| Arcieri 2023/2024 | State-only belief treated as sufficient under autoregressive observations. | Mathematical state-definition gap. |
| Arcieri 2024/2026 | Persistent replay buffer described for on-policy PPO. | Published algorithm incomplete/incompatible as written; code absent. |
| Arcieri 2024 | Belief filter may receive sampled true episode parameters. | Serious unresolved oracle risk; code absent. |
| Arcieri 2026 | Eq. 16 multiplies prior and an already-defined posterior. | Apparent normalization/double-counting error. |
| Arcieri 2026 | Generative transition propagates beliefs rather than realized state. | Apparent failure to represent the state-trajectory joint. |
| Song | Likelihood minimized; failure cost/units and corrosion equations conflict. | Confirmed internal inconsistencies. |
| Kim | High-threshold validation weak; Chain-Bernoulli details and weekly-to-monthly conversion are underspecified; figure refs are wrong. | Confirmed reporting/reproducibility limitations. |
| Jackson manual | Matrix-exponential series omits \(X\); other table/output inconsistencies. | Confirmed manual errors; standard CTMC theory remains sound. |
| Bashar & Torres | Eq. 4 index error; infinite-horizon policy evaluated/interpreted over 25 years. | Confirmed equation/objective inconsistencies. |

---

## 10. A defensible research path from this repository to HAB VoI

This section is a methodological recommendation, not an implementation. The ordering matters: a sophisticated policy solver cannot compensate for an undefined decision or an uncalibrated observation model.

### 10.1 Stage 0: write the decision contract before choosing the model

The project needs a one-page, externally reviewable contract that answers:

1. **Who decides?** A regulator, shellfish grower, public-health officer, or other operator may face different constraints.
2. **What asset/outcome is protected?** Human exposure, a particular harvest area, a shellfish lot, an early-harvest opportunity, a closure/reopening decision, or sampling allocation are not interchangeable.
3. **What is the event?** For example, tissue DA at/above 20 µg/g within one week, water particulate DA above a threshold, or a latent toxic-bloom regime.
4. **What actions are legally and operationally available?** C-HARM may inform confirmatory sampling but may not legally replace tissue testing. Closure may be mandatory after a laboratory result.
5. **What is the exact within-epoch timing?** Daily nowcast/forecast, field collection, lab turnaround, harvest, advisory, closure, and reopening must be ordered.
6. **What is the baseline policy?** “No remote sensing” is not the same as “no monitoring.” The current CDPH/HABMAP schedule and event-triggered behavior must be described.
7. **What is the horizon?** One-step early warning, a seasonal finite horizon, and an infinite-horizon stationary policy answer different questions.
8. **What costs/consequences count?** Illness severity, false negatives, unnecessary closures, lost harvest, early harvest quality/price, lab tests, staff time, fixed data integration, and public trust may matter.
9. **What constraints dominate expected cost?** A chance constraint on unsafe harvest or a risk measure may be more appropriate than average economic value alone.

A useful initial target might be:

> Estimate the value of adding C-HARM nowcast/forecast scores to the current field-monitoring history for deciding whether to order confirmatory tissue sampling and/or take a time-sensitive harvest action over the next one to four weeks at a specified site.

That statement is only an example. It is deliberately narrower and more testable than “the value of C-HARM for HABs.”

### 10.2 Stage 1: repair the data contract

#### Preserve measurement definitions

- Keep the two *Pseudo-nitzschia* size groups separate unless a documented common-denominator series can be constructed.
- Distinguish missing, below detection, true zero, and not-assayed.
- Retain exact collection date/time, laboratory result date, assay method, LOD/LOQ, shellfish species, replicate/sample count, and location coordinates.
- Retain water particulate/total DA and other mechanistically relevant variables rather than reducing the modeling panel to cell total and tissue category.
- Document site-name mapping and test that expected coverage is not silently lost.

#### Define temporal alignment scientifically

- Select an explicit local-time/UTC convention and decision-week boundary.
- Keep within-week dates or daily resolution when lead/lag is the research question.
- If weekly aggregation remains, justify mean, maximum, last value, or time-weighted summary separately for each channel.
- Estimate distributed lags between water cells, water DA, C-HARM scores, and tissue DA rather than prespecifying same-week alignment.

#### Model or audit the observation process

- Quantify sampling frequency by site, season, prior measurement, alert/closure, and year.
- Determine whether field or tissue sampling is triggered by visible conditions, complaints, preliminary results, or external forecasts.
- If missingness is informative, include a sampling/selection model or restrict the estimand to the historical monitoring process.

### 10.3 Stage 2: build a genuine C-HARM matchup and observation model

#### Matchup design

For each product and lead time, specify:

- which pixel/neighborhood/offshore footprint represents each management area;
- how coastline masks and missing/gap-filled fields are handled;
- whether the 0-, 1-, 2-, or 3-day product is available before the decision;
- how a daily score maps to field and tissue collection dates; and
- which time/site blocks are reserved before any calibration choices.

The local overlap begins only in November 2022. With roughly 146 same-week cell+tissue observations in the overlap era across four sites, flexible calibration must be strongly regularized and evaluated honestly. Spatially correlated daily pixels are not thousands of independent validation cases.

#### Preserve the full observation

C-HARM’s three output probabilities should initially remain continuous and distinct. Calibration questions include:

- Does a nominal 0.7 score correspond to a 70% event frequency in the relevant site/season/lead setting?
- What are Brier score, log score, calibration intercept/slope, reliability curve, and precision–recall behavior for rare closure events?
- How do false-negative and false-positive rates change with threshold?
- Does the score add predictive information after current field/lab history and seasonal climatology?

For POMDP use, estimate a distribution such as \(P(C_t\mid X_t)\), where \(C_t\) is the C-HARM score vector and \(X_t\) is the decision-relevant latent state/outcome. Do not insert the score itself as though \(P(C_t\mid X_t)=C_t\). If a discrete observation matrix is necessary, binning and operating thresholds must be learned/tuned inside training data and evaluated on untouched data.

#### Avoid double counting

C-HARM already incorporates satellite reflectance/chlorophyll, SST, salinity, circulation, and calendar-related empirical structure. Adding C-HARM plus the same inputs as conditionally independent sensor channels reuses evidence. A joint likelihood or explicit conditional structure is needed.

### 10.4 Stage 3: choose an ecological model proportionate to the decision

#### Minimum credible near-term model

If the first decision is a one- to four-week test/harvest action, a calibrated dynamic outcome model may be enough. It can predict future tissue exceedance from past tissue, water cells/DA, season, site, and C-HARM scores without claiming that two hidden labels are true ecological states. A one-step or short-horizon decision analysis can calculate information value without a full infinite-horizon POMDP.

This is an important scientific option. A POMDP is valuable when sequential hidden-state dynamics and future information/action tradeoffs materially affect policy; it should not be adopted solely because the reference repository is about POMDPs.

#### Structured state-space model for a full POMDP

If a latent ecological model is needed, a more plausible structure separates at least:

1. **water-column biomass** \(B_t\);
2. **toxigenicity/water toxin** \(W_t\); and
3. **shellfish tissue burden** \(H_t\), with uptake and depuration from past \(W\).

Calendar season and forecast environmental variables are observed context. Site effects/hierarchies account for different baselines and measurement practice. C-HARM bloom probability primarily observes/predicts \(B\); particulate/cellular products relate more closely to \(W\); tissue assays observe \(H\). Conditional independence becomes more plausible after this separation.

A stylized dependency is:

\[
B_t\longrightarrow W_t\longrightarrow H_{t+1:t+L},
\]

with season/environment affecting \(B_t,W_t\), and site/species affecting observations and tissue kinetics. This can be continuous, discrete, factorial, switching, or semi-Markov; the choice should be driven by data and decision needs.

#### Action semantics

Many likely actions do not alter bloom ecology:

- consuming a passive C-HARM product changes information, not \(B_t\);
- ordering a test changes information and cost;
- closing harvest changes exposure/economic state, not ocean dynamics;
- early harvest changes asset/terminal status and payoff.

Those actions need not create action-conditioned ecological transition matrices. The decision state may instead include whether the resource is open, harvested, closed, or terminal. Action-conditioned ecological transitions should be introduced only for interventions that physically affect the hazard and have causal/physical support.

### 10.5 Stage 4: impose validation gates before policy work

#### Gate A: computation

- All chains satisfy prespecified \(\hat R\), ESS, divergence, and energy criteria.
- Multiple initializations reach the same modes.
- Likelihood and belief updates pass independent algebraic tests.
- Tail probabilities and state order are sensible over the whole operational range.

#### Gate B: posterior predictive adequacy

Replicated data must match site/era-specific zeros, quantiles, extremes, tissue categories, lag correlations, episode durations, seasonal patterns, and missingness. Checks should be stratified rather than pooled.

#### Gate C: prospective prediction

- Rolling-origin leave-future-out validation uses only past data.
- Forecast horizon matches the decision horizon.
- A final chronological block is untouched until model lock.
- Held-out-site validation is used if geographic transfer is claimed.
- Baselines include seasonal climatology, last-observation/persistence, and simple regression/state-space alternatives.

#### Gate D: rare-event calibration

Closure events need sensitivity, precision/positive predictive value, false-negative rates, calibration, and proper scores with uncertainty. Overall accuracy is dominated by non-detects and is not sufficient.

#### Gate E: state meaning

If state labels enter rewards, their biological/decision meaning must survive site, season, protocol, and lag checks. If states are merely predictive compression, report them that way and map rewards to predicted consequences rather than names like “safe” or “toxic.”

### 10.6 Stage 5: formulate the POMDP with correct timing

A useful decision epoch can be represented as two substeps:

1. **Information substep:** passively receive C-HARM; optionally commission field/lab sampling if that is a choice.
2. **Management substep:** after updating belief, harvest, wait, issue an advisory, close, maintain closure, reopen, or take another feasible action.

The action set must respect legal rules. If C-HARM cannot authorize reopening, the solver must not be allowed to do so.

The information state must include every observed context needed for the future distribution. If emissions/transitions depend on recent observations, calendar week, site, open/closed status, or test pipeline, the policy input is not just a vector over ecological labels.

For a small discrete model, exact finite-horizon dynamic programming or a transparent point-based solver should be the primary method. QMDP is a benchmark, not the main VoI solver, because it assumes uncertainty resolves after one step and cannot properly value tests/information actions. Deep RL should be considered only after a small interpretable solver and simulation tests establish the objective.

Optimize and evaluate the same horizon. If the real season has a terminal date, use a finite-horizon policy with calendar time in state rather than an infinite-horizon policy evaluated for one season.

### 10.7 Stage 6: compute VoI as a nested policy comparison

Recommended regimes include:

1. **Current-practice baseline:** actual scheduled/event-triggered field and tissue monitoring, without C-HARM.
2. **Baseline + C-HARM nowcast:** passive current remote information.
3. **Baseline + C-HARM forecasts:** lead-specific future products.
4. **C-HARM-triggered confirmatory sampling:** remote information changes testing allocation, then tests change management.
5. **Sparse-monitoring setting:** relevant only if an actual management area lacks the baseline service.

The with-information action set should contain the baseline actions. Gross optional VoI is then nonnegative under exact optimization. Report separately:

- gross expected decision value;
- per-use data/test cost;
- unavoidable fixed integration/adoption cost;
- net value after those costs; and
- break-even fixed cost.

Use paired/common-random-number policy simulations across the same posterior world draws to reduce Monte Carlo noise. Propagate uncertainty in ecological parameters, C-HARM calibration, rewards, baseline policy, and structural alternatives. Report the distribution of incremental value, not only its mean.

If public-health downside is asymmetric, supplement expected value with unsafe-harvest probability, regret, lower quantiles, or CVaR. Calling posterior-mean optimization “robust” is not enough.

### 10.8 Stage 7: decision-focused sensitivity analysis

The most informative sweeps are not a single generic accuracy variable. Vary:

- sensitivity and specificity separately;
- calibration slope/intercept and prevalence shift;
- lead time;
- spatial mismatch;
- tissue accumulation/depuration lag;
- assay LOD and turnaround time;
- false-negative health cost and false-positive closure cost;
- harvest value and timing;
- baseline sampling frequency;
- fixed versus marginal information cost;
- model/state count and site hierarchy; and
- risk tolerance.

Each sensitivity should state exactly which conditional distribution or cost changes. A scalar “accuracy” sweep is acceptable only as a deliberately synthetic channel experiment with no claim of universal threshold.

### 10.9 Stage 8: reproducibility and governance

Before results are used externally, the project needs:

- a current scope README and mathematical methodology;
- a data dictionary and decision-timing diagram;
- pinned Python/R environments and an artifact manifest;
- immutable raw data plus versioned derived datasets;
- restart-and-run workflows rather than stateful notebook outputs;
- tests for data mappings, week/lead alignment, row-stochastic transitions, matrix exponentials, forward likelihood, belief update, periodic occupancy, and nested nonnegative gross VoI;
- posterior predictive and validation reports attached to each export;
- provenance inside exports: model version, data hash, dates, transforms, category definitions, state semantics, diagnostics, and initial/current belief; and
- an explicit label on exploratory/nonconverged/stale artifacts.

### 10.10 Suggested go/no-go sequence

| Gate | Required evidence before proceeding |
|---|---|
| 1. Decision definition | Stakeholders agree on action, outcome, timing, baseline, horizon, and loss. |
| 2. Data comparability | PN-group, LOD, site, date, species, coordinates, and monitoring-process issues are resolved or modeled. |
| 3. World-model validity | Prospective/site validation and posterior predictive checks pass; hazard ordering is not pathological. |
| 4. C-HARM incremental skill | Held-out C-HARM scores add calibrated predictive information beyond current monitoring/climatology. |
| 5. Policy validity | Transparent solver, timing, constraints, and baselines are verified in simulation across posterior/structural uncertainty. |
| 6. Decision value | Net value remains favorable under cost, calibration, lag, and risk sensitivities. |
| 7. Operational evidence | A prospective shadow-mode pilot confirms data availability, calibration, and action feasibility before live reliance. |

---

## 11. What can and cannot be concluded today

### Defensible conclusions

1. The saved interim panels are reproducible from the current data builders and have unique site/week keys.
2. A pooled two-state stationary CT-HMM can be sampled reliably under the current parameterization.
3. A pooled two-state one-harmonic seasonal CT-HMM can also be sampled reliably.
4. The custom matrix exponential is numerically accurate over the saved posterior range, and the 52-week gap cap is negligible for the current two-state rate posterior.
5. The pooled observations exhibit an in-sample seasonal pattern that the fitted onset-rate harmonic represents strongly; site/protocol composition and future predictive skill remain unresolved.
6. A tissue-only two-state CT-HMM provides a plausible descriptive split between a long non-detect regime and short toxic episodes, with AIC/BIC favoring it over the fitted three/four-state alternatives.
7. Explicit belief filtering and posterior parameter uncertainty remain appropriate conceptual foundations for a later decision model.

### Conclusions that are not yet defensible

1. That two is the true or predictively optimal number of ecological states.
2. That state 1 in the Python model is a monotone toxic-bloom or closure-risk state.
3. That the plotted seasonal curve is actual bloom occupancy or that its week-17 peak is ecological prevalence onset.
4. That the seasonal model predicts future held-out HAB conditions better than the stationary model.
5. That the model generalizes across sites, protocol eras, or unobserved future years.
6. That cells and tissue DA are conditionally independent same-week signals of one state.
7. That C-HARM is calibrated against these sites/outcomes, adds incremental information, or has any estimated break-even accuracy.
8. That any management policy, early-harvest rule, closure strategy, or sampling allocation is optimal.
9. That remote sensing has positive, negative, or zero economic/public-health value in this application.
10. That the Arcieri tutorial notebook outputs validate their visible code, and especially not the 2024 RL or 2026 DBMM methods that are absent here.

### Overall judgment

This is a promising exploratory **world-model prototype** with a clear intended research direction. It is not a failed POMDP; it is a POMDP that has not yet been built. The immediate scientific priority is not adding a more sophisticated solver or neural policy. It is making the hazard state, data channels, lag structure, site/protocol effects, C-HARM observation model, and prospective validation trustworthy.

The two-state stationary output should not currently cross the boundary into a policy simulator. Its extreme-tail reversal alone can make more alarming cell counts reduce inferred high-state probability. The seasonal output should likewise not be used to initialize a calendar policy until calendar-conditioned beliefs—or a clearly defined periodic approximation—are propagated correctly and held-out seasonal skill is shown.

Once those foundations are repaired, the literature supports a strong design: a full calibrated belief, an observe-update-act epoch, a transparent finite-horizon policy, nested monitoring regimes, and posterior/structural uncertainty carried into the distribution of incremental decision value.

---

## Appendix A. End-to-end execution narrative

### A.1 Panel construction

Running the panel entry point calls five builders in sequence:

1. calHABMAP rows are parsed, two PN columns are summed, and all available hazard/physical fields are averaged by site/week.
2. CDPH rows are mapped to four sites, classified by censoring/20-µg/g threshold, and reduced to weekly maxima.
3. A four-site cell panel retains only weeks with constructed PN total and transforms it to log10(cells+1).
4. A tissue panel assigns one integer subject per site and elapsed years for R.
5. The cell and tissue panels are outer-joined, preserving weeks in either channel.

The script prints counts and writes a coverage heatmap. It does not touch C-HARM.

### A.2 Tissue reference fit

The R script reads the tissue panel, fits \(K=2,3,4\) CT-HMMs independently from one starting configuration each, prints Q/sojourn summaries, saves R objects, and writes AIC/BIC. It does not export its transition/emission estimates into the Python world model or decision layer.

### A.3 Stationary joint fit

The stationary script reads the joint panel, creates four sequences, and for each requested \(K\):

1. constructs a birth-death generator and ordered cell/tissue emissions;
2. uses NUTS/nutpie to sample continuous parameters while the forward algorithm sums over states;
3. reconstructs forward likelihood factors;
4. writes a NetCDF trace;
5. calculates the disputed PSIS score and sampling diagnostics; and
6. adds the fit to an all-model comparison.

It then excludes fits based on an ad hoc \(\hat R/p_{\text{LOO}}\) gate and exports the highest-ELPD survivor. In the saved run, \(K=2\) is the only survivor.

### A.4 Seasonal fit

The seasonal script expands each site from its first to last observed week, fits \(K=2,H=1\), saves the trace, reconstructs likelihood factors in memory, compares with the stationary trace, and exports coefficient draws plus a plug-in weekly transition cycle. It computes the frozen-rate ratio and calls it occupancy. The seasonal figures live in visualization functions and are not invoked by the seasonal fit script itself; they were generated separately.

### A.5 Intended but absent continuation

The export docstring calls the NPZ file the “world-building boundary” consumed by a decision model. No current file implements that consumer. A future continuation would need to load a versioned model, establish a current belief from historical observations, ingest a new C-HARM observation through a calibrated likelihood, select an action from a policy, and update both ecological and management status.

---

## Appendix B. Numerical audit summary

### B.1 Data

| Quantity | Value |
|---|---:|
| Raw calHABMAP rows | 7,203 |
| Raw calHABMAP sites | 17 |
| calHABMAP weekly rows | 7,004 |
| Selected raw CDPH samples | 1,159 |
| Selected CDPH site-weeks | 1,110 |
| Cell panel weeks | 2,996 |
| Joint panel weeks | 3,548 |
| Joint same-week pairs | 558 |
| Closure-category weeks | 34 |
| Closure weeks with same-week cells | 20 |
| Exact cell zeros | 914 |
| Maximum cell count | ~3.25 million/L |
| Gaps over 52 weeks | 11 |
| Maximum gap | 255 weeks |
| C-HARM-era joint weeks | 683 |
| C-HARM-era cell weeks | 674 |
| C-HARM-era tissue weeks | 155 |
| C-HARM-era same-week cell+tissue | 146 |

### B.2 Saved stationary fits

| Model | Computation status | Scientific status |
|---|---|---|
| \(K=2\) | Converged; max \(\hat R=1.013\), no divergences | Usable exploratory posterior, but emission/state meaning fails key predictive checks. |
| \(K=3\) | Not converged; max \(\hat R=1.544\) | Invalid for inference/comparison. |
| \(K=3\), longer six-chain run | Not converged; max \(\hat R=1.623\) | Confirms structural multimodality/identification difficulty. |
| \(K=4\) | Not converged; max \(\hat R=1.536\) | Invalid for inference/comparison. |

### B.3 Saved seasonal fit

| Quantity | Value/interpretation |
|---|---|
| \(K,H\) | 2 states, 1 annual harmonic |
| max \(\hat R\) | 1.011 |
| divergences | 0 |
| posterior-mean onset range | ~0.054–0.243 per week |
| posterior-mean recovery range | ~0.081–0.099 per week |
| exported frozen ratio range | ~35.6–74.8% |
| propagated 53-bin plug-in high-state probability | ~39.6–73.6% |
| exported/correct peak | week 17 / week 20 |

---

## Appendix C. Reproducibility observations

- The root README is empty and there is no current methodology document.
- No tests are present.
- Python and R dependencies are unpinned; there is no lockfile.
- Direct SciPy/xarray dependencies are undeclared.
- The local base virtual environment lacks the optional inference stack. Saved posterior metadata indicate a different environment with Python 3.13.7, PyMC 6.1.0, nutpie 0.16.11, and ArviZ 1.2.0.
- The Arcieri notebooks require old PyMC3/Theano environments and have no fixed random seeds.
- The deterioration notebook’s visible cell source and saved output execution order are inconsistent.
- Output directories contain exploratory, stale, nonconverged, and current artifacts without a manifest.
- The local C-HARM products are large external data assets but are not version-linked to any derived matchup.
- The review preserved all pre-existing working-tree changes. The uncommitted seasonal implementation was treated as part of the reviewed snapshot, not as an established released component.

---

## Appendix D. Glossary

**Action-conditioned transition:** A transition distribution that changes with the chosen action. Association in logged data is not automatically a causal action effect.

**Belief:** A probability distribution over the hidden state given the available history and chosen model.

**Calibration:** Agreement between predicted probabilities and observed event frequencies. A calibrated 70% prediction should occur roughly 70% of the time in the relevant setting.

**Conditional independence:** Two observations are independent after the latent state is known. This is a modeling assumption; biological lag can violate it.

**CTMC:** Continuous-time Markov chain, parameterized by an intensity/generator matrix \(Q\).

**Emission/observation model:** The distribution of a measurement or forecast output conditional on the latent state.

**ELPD:** Expected log predictive density. Its interpretation depends entirely on which data are genuinely held out.

**HMM:** Hidden Markov model: latent Markov states plus an observation model.

**Informative sampling:** The chance/timing of observing data depends on unmodeled state or risk.

**LOD:** Limit of detection. A result below 1 and below 2.5 convey different censored information.

**MDP:** Markov decision process with observed state, actions, transitions, and rewards.

**POMDP:** Partially observable MDP. A policy usually acts on a belief rather than the unobserved true state.

**Posterior predictive check:** Comparing simulated replicated observations from the fitted posterior with real observations.

**PSIS-LOO:** Pareto-smoothed importance-sampling approximation to leave-one-out cross-validation. Stable weights do not rescue a wrongly defined holdout unit.

**QMDP:** A heuristic that acts under current uncertainty but assumes the state becomes fully observed after one step; it undervalues information-gathering actions.

**Semi-Markov model:** A model in which transition behavior can depend on time already spent in the current state.

**Structural/model-form uncertainty:** Uncertainty about equations, state dimensions, lag form, emission family, or omitted mechanisms—not merely parameter uncertainty inside one model.

**Tobit model:** A censored continuous model that assigns a probability mass to observations at a censoring boundary.

**Value of information:** Improvement in optimized decision value when an information source is available, under a specified action set, cost/reward model, baseline, horizon, and observation process.

**Viterbi path:** The single joint hidden-state path with highest posterior/probability under an HMM; it is not the sequence of marginal posterior modes and does not preserve state uncertainty.

---

## Appendix E. Sources and precise use in this review

### Supplied sources

1. Arcieri, G., et al. (2023). “Bridging POMDPs and Bayesian decision making for robust maintenance planning under model uncertainty.” *Reliability Engineering & System Safety* 239, 109496. Used for the Bayesian HMM, posterior-averaged control, QMDP, railway data/model/results, and the Algorithm 1 audit. [Attached PDF](</Users/angelinale/Downloads/Bridging POMDPs and Bayesion Decision Making (Arcieri, 2023).pdf>).
2. Song, C., et al. (2022). “Value of information analysis in non-stationary stochastic decision environments: A reliability-assisted POMDP approach.” *Reliability Engineering & System Safety* 217, 108034. Used for layered nonstationarity, two-substep timing, observation derivation, beam results, and internal-consistency audit. [Attached PDF](</Users/angelinale/Downloads/Non-Stationary POMDP for VoI (Song).pdf>).
3. Arcieri, G., et al. (2026). “Deep belief Markov models for POMDP inference.” *Neural Networks* 196, 108386. Used for DBMM factorization, experiments, RL claims, and equations/appendix audit. [Attached PDF](</Users/angelinale/Downloads/Deep Belief Markov Models for PODMDP (Arcieri, 2026).pdf>).
4. Arcieri, G., et al. (2024). “POMDP inference and robust solution via deep reinforcement learning: an application to railway optimal maintenance.” *Machine Learning* 113, 7967–7995. Used for belief/LSTM/GTrXL PPO comparison, domain randomization, results, and Algorithm 1 audit. [Attached PDF](</Users/angelinale/Downloads/POMDP Inference & Deep RL (Arcieri, 2024).pdf>).
5. Kim, K. B., et al. (2022). “A multivariate Chain-Bernoulli-based prediction model for cyanobacteria algal blooms at multiple stations in South Korea.” *Environmental Pollution* 313, 120078. Used for covariate-driven network HMM, validation results, predictor experiment, and reproducibility audit. [Attached PDF](</Users/angelinale/Downloads/multivariate HAB prediction (Kim, 2022).pdf>).
6. Jackson, C. (2024). *Multi-state modelling with R: the msm package*, v1.8.2. Used for CTMC theory, observation timing, hidden emissions, diagnostics, and the `msm` implementation context. [Attached manual](</Users/angelinale/Downloads/multi-state modelling (Jackson).pdf>). Peer-reviewed foundation: [Jackson (2011)](https://doi.org/10.18637/jss.v038.i08).
7. Bashar, M., & Torres-Machi, C. (2023). “Quantifying the Value of Satellite-Based Pavement Monitoring in Partially Observable Stochastic Environments.” *Journal of Computing in Civil Engineering* 37(3), 04023004. Used for empirical satellite observation matrices, nested monitoring comparisons, accuracy sensitivity, policy results, and timing/equation audit. [Attached PDF](</Users/angelinale/Downloads/VoI in POMDP (Bashar & Torres).pdf>).

### Focused external sources

8. Bürkner, P.-C., Gabry, J., & Vehtari, A. (2020). “Approximate leave-future-out cross-validation for Bayesian time series models.” *Journal of Statistical Computation and Simulation* 90(14), 2499–2523. Used specifically for the forecasting estimand and why future-informed LOO is optimistic (§2, §3.1, §4.1, §5). [Open manuscript](https://arxiv.org/pdf/1902.06281); [DOI](https://doi.org/10.1080/00949655.2020.1783262).
9. Pohle, J., Langrock, R., van Beest, F. M., & Schmidt, N. M. (2017). “Selecting the Number of States in Hidden Markov Models.” *Journal of Agricultural, Biological and Environmental Statistics* 22, 270–293. Used specifically for extra states absorbing misspecification and for the pragmatic selection workflow (§§2.2–2.3, §3, §4, §6). [Open manuscript](https://arxiv.org/pdf/1701.08673); [DOI](https://doi.org/10.1007/s13253-017-0283-8).
10. Anderson, C. R., et al. (2016). “Initial skill assessment of the California Harmful Algae Risk Mapping (C-HARM) system.” *Harmful Algae* 59, 1–18. Used for C-HARM product construction/targets, spatial and temporal matchup, multi-metric skill, site dependence, and water-to-tissue lead (§1, §§2.1–2.6, §3.2, §§4.1 and 4.3). [NOAA full text](https://repository.library.noaa.gov/view/noaa/33076/noaa_33076_DS1.pdf); [DOI](https://doi.org/10.1016/j.hal.2016.08.006).
11. Lane, J. Q., Roddam, C. M., Langlois, G. W., & Kudela, R. M. (2010). “Application of Solid Phase Adsorption Toxin Tracking (SPATT) for field detection of the hydrophilic phycotoxins domoic acid and saxitoxin in coastal California.” *Limnology and Oceanography: Methods* 8, 645–660. Used specifically for early-warning lead, water/tissue uncoupling, and uptake/depuration caveats (Abstract, pp. 655–659). [DOI](https://doi.org/10.4319/lom.2010.8.0645); [open reproduced article](https://escholarship.org/content/qt3hk1d1sf/qt3hk1d1sf_noSplash_30a5d9ae7d63c41fbda65407135fb4e5.pdf).

### Review boundaries

The review visually inspected important equation, algorithm, table, and figure pages in every supplied PDF and used full-text extraction for systematic searching. It reviewed notebook source and embedded outputs rather than trusting markdown summaries. Saved posterior artifacts were inspected and checked algebraically/numerically, but the expensive Bayesian and R models were not re-fit. The 2024 PPO and 2026 DBMM implementation code is not in this repository; findings about possible oracle filtering or actual PPO buffers are therefore explicitly labeled unresolved risks rather than confirmed code defects.
