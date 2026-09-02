# TrialEmulationj 1.0

TrialEmulationj is a jamovi module for designing, auditing, estimating, and
reporting causal analyses that emulate treatment strategies with observational
data. Version 1.1 keeps design decisions visible and separates workflows whose
data structures, estimands, and identifying assumptions are not interchangeable.

- A target-trial protocol builder with readiness checks and a SHA-256
  specification fingerprint.
- Prespecified causal DAG validation, minimally sufficient adjustment sets,
  selected-covariate checks, and a static DAG figure.
- Audited 24/48-hour landmark analyses for two treatment components and four
  observed strategies.
- Labelled binary variables as well as 0/1 coding.
- AIPW and Hajek IPTW estimates for risk difference, risk ratio, odds ratio, or
  a continuous mean difference.
- Cluster-robust influence-function inference and optional whole-cluster
  bootstrap intervals.
- Automatic or requested Firth mean bias-reduced logistic regression when
  separation makes ordinary logistic estimates unstable.
- Explicit blank-as-absence recoding, complete-case analysis, or stochastic
  hot-deck multiple imputation, each with an audit trail.
- Propensity overlap, effective sample size, truncation sensitivity, covariate
  balance (Love plots), cohort flow, and combined effect/forest plots.

## Analyses in the Causal Inference menu

### 1. Target Trial Protocol

Define the target population, eligibility, strategies, assignment procedure,
time zero, grace period, follow-up, outcome, competing events, causal contrast,
estimand, confounding-control plan, missingness, and sensitivity analyses before
modelling. A readiness table flags incomplete design decisions. The fingerprint
identifies the entered specification; it does not hash patient data.

### 2. Causal DAG and Adjustment Set

Enter a DAGitty graph, select the exposure and outcome columns, and optionally
select the proposed adjustment variables. The analysis checks acyclicity,
returns minimal sufficient sets, flags exposure descendants and simple collider
patterns, and draws a publication-ready DAG. The graph is prespecified; the
module does not learn causal structure from the dataset.

### 3. Landmark Trial (24/48 h)

Construct aligned cohorts at one or two decision windows. Version 1.1 can:

- classify two treatment components into four strategies at each landmark;
- analyze one or both strata of treatment component A while comparing component
  B exposed versus unexposed;
- exclude treatment, outcome, competing-event, or discharge events occurring
  before the landmark;
- audit missing, negative, or inconsistent treatment and event times;
- estimate AIPW and Hajek IPTW effects with clustered uncertainty; and
- output flow, balance, and effect plots.

For an NGT/WSEC study, component A can represent NGT and component B WSEC, with
24- and 48-hour windows. Supply event-time variables whenever failure can occur
before a landmark; otherwise outcome-free status at time zero cannot be verified.

### 4. Sequential Target Trial

Delegates longitudinal sequential-trial emulation to
[TrialEmulation](https://github.com/Causal-LDA/TrialEmulation). It supports
marginal initiation/assignment analogues and per-protocol estimands, stabilized
adherence and censoring weights, pooled logistic marginal structural models,
cluster-robust standard errors, marginal curves, and absolute contrasts. Input
must contain one row per person-period with unique ID-period combinations.

### 5. Clone-Censor Weighting

Delegates cloning, artificial censoring, long-data expansion, and censoring
weights to [survivalCCW](https://github.com/Genentech/survivalCCW). It supports
two grace-period strategies, weight winsorization, weighted survival contrasts,
RMST differences, and a person-level bootstrap that keeps each person's clones
together. Input must contain one baseline row per person.

### 6. AIPW / Doubly Robust

Estimates marginal point-outcome effects using AIPW and normalized IPTW. Choose
pooled or arm-specific outcome models; independently specify outcome and
propensity covariates; request cluster-aware inference, bootstrap intervals,
Firth fitting, and missing-data handling; and inspect overlap, weights, effective
sample size, balance, and propensity-truncation sensitivity.

AIPW is doubly robust to nuisance-model misspecification, not to violations of
the causal assumptions. Consistency still requires a correct propensity model
or a correct outcome model, plus a well-defined intervention, exchangeability,
positivity, consistency, and correct temporal ordering.

### 7. Doubly Robust Survival

Delegates to `precmed::atefitsurv()` for treatment, outcome, and censoring
models. It reports treatment-specific restricted mean survival time, a
restricted-mean-time-lost ratio, an adjusted hazard ratio, bootstrap intervals,
and diagnostic plots from one-row-per-person survival data.

## Install in jamovi

1. Open jamovi.
2. Open the module library and choose **Sideload**.
3. Select the version 1.1 `.jmo` file.
4. Open **Causal Inference** and begin with **1. Target Trial Protocol**.

For development builds, install R, jamovi, and
[jmvtools](https://github.com/jamovi/jmvtools), then prepare and install this
package from its source directory.

## Example and validation data

`example-data/data_censored.xlsx` is a synthetic longitudinal example for the
Sequential Target Trial analysis. A starting mapping is:

- ID: `id`
- period: `period`
- treatment: `treatment`
- outcome: `outcome`
- eligible: `eligible`
- outcome covariates: `age_s`, `x1`, `x2`
- informative censoring: `censored`

Automated tests generate separate synthetic datasets for AIPW, Firth
separation, DAGs, and clustered 24/48-hour landmarks. No patient-level SnapSBO
data are distributed with the module.

## Interpretation and safeguards

The module cannot establish that a causal model is correct. Prespecify the
target trial and DAG; align eligibility, treatment classification, and follow-up
at time zero; distinguish structural infeasibility from sparse empirical
positivity; and use baseline causes rather than mediators, colliders, or
post-treatment variables. Record active jamovi filters and review timing audits,
overlap, effective sample size, extreme weights, balance, model warnings, and
sensitivity analyses before interpreting an effect.

Blank-as-absence recoding is appropriate only when the data-collection process
defines a blank as absence. Hot-deck imputation is a transparent pragmatic
sensitivity method, not a substitute for a protocol-specific imputation model
in a confirmatory analysis. An observational treatment-initiation contrast is
not literally an intention-to-treat effect unless randomized assignment is
observed.

## Method references

- Hernan MA, Robins JM. *Causal Inference: What If*.
- Maringe C, et al. Target trial emulation methodology and the TrialEmulation R
  package.
- Textor J, et al. DAGitty for causal diagrams and covariate adjustment.
- Huntington-Klein's *The Effect*, its replication repository, and the
  `causaldata` teaching examples are useful worked examples; they are not
  runtime dependencies.

Before public distribution, replace the placeholder maintainer email in
`DESCRIPTION` with the project's support address.
