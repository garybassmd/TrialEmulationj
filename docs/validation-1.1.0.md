# TrialEmulationj 1.1 validation report

## Release

- Package: `trialemulationj`
- Version: 1.1.0
- Jamovi build: macOS Apple Silicon (`arm64`)
- Jamovi application used for compilation: 28.2.0
- Jamovi R runtime used for smoke testing: R 4.6.0
- Patient-level data bundled: none

## Included analyses

1. Target Trial Protocol Builder
2. Causal DAG and Adjustment Sets
3. Landmark Target-Trial Analysis (one or two decision windows)
4. Sequential Target Trial Emulation
5. Clone-Censor-Weight Survival Analysis
6. AIPW and Hajek IPTW for Point Outcomes
7. Doubly Robust Survival Effect

## Verification completed

- Jamovi definitions compiled for all seven analyses.
- R package source parsed and installed successfully.
- `R CMD check --no-manual` result: **Status: OK**.
- Automated testthat result: **18 passed, 0 failed, 0 warnings, 0 skipped**.
- The final `.jmo` was unpacked and loaded using Jamovi's own R 4.6.0 runtime.
- Runtime smoke tests passed for the protocol builder, DAG, AIPW/Firth,
  clustered 24/48-hour landmarks, sequential target trials, clone-censor
  weighting, and doubly robust survival.
- Bundled dependencies were confirmed for `brglm2`, `dagitty`, `digest`,
  `TrialEmulation`, `survivalCCW`, and `precmed`.

## Version 1.1 safeguards exercised

- Labelled binary levels and 0/1 variables.
- Automatic/requested Firth mean bias reduction under logistic separation.
- Cluster-robust influence-function standard errors.
- Whole-cluster bootstrap infrastructure.
- Landmark alignment and exclusions for treatment, failure, competing events,
  and discharge occurring at or before time zero.
- Missing, unparseable, negative, and inconsistent time audits.
- Explicit blank-as-absence recoding audit.
- Complete-case and hot-deck multiple-imputation paths.
- Propensity truncation, overlap, inverse-weight, effective-sample-size, and
  standardized-mean-difference diagnostics.
- Love, overlap, flow, effect/forest, and DAG plots.
- SHA-256 specification fingerprints and copyable R specifications.

## Important limitations

- This binary is for macOS Apple Silicon. Windows, Linux, and Intel macOS need
  platform-specific Jamovi builds.
- A causal DAG must be prespecified from subject-matter knowledge; the module
  does not discover causal structure from observed associations.
- Hot-deck imputation is a transparent pragmatic sensitivity method. A
  confirmatory study should use a protocol-specific imputation model and
  sensitivity analysis.
- Short 20-sample bootstraps were used only for runtime testing. Applied work
  should prespecify a substantially larger number appropriate to the study.
- Observational treatment initiation is not literally an intention-to-treat
  effect unless randomized assignment is observed.
- Causal interpretation still requires a well-defined intervention,
  exchangeability, positivity, consistency, correct temporal ordering, and
  appropriate measurement.
- The upstream `survivalCCW` engine should be treated as evolving software and
  results should be checked against protocol-specific code for consequential
  analyses.
- The placeholder maintainer email in `DESCRIPTION` should be replaced before
  public distribution.

## Install

Open Jamovi, open the module library, choose **Sideload**, and select
`trialemulationj-1.1.0-macos-arm64.jmo`. The analyses appear under
**Causal Inference**.

## Reference material used in design

- Hernan and Robins, *Causal Inference: What If*.
- The TrialEmulation methodology and R package documentation.
- DAGitty documentation for graph validation and adjustment sets.
- `brglm2` documentation for bias-reduced logistic regression.
- *The Effect*, its causalbook replication repository, and `causaldata` as
  worked teaching references rather than runtime dependencies.
