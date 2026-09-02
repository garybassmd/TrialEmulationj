# TrialEmulationj 1.1.1 validation report

## Release

- Package: `trialemulationj`
- Version: 1.1.1
- Jamovi build: macOS Apple Silicon (`arm64`)
- Jamovi application used for compilation: 28.2.0
- Jamovi R runtime used for bundle testing: R 4.6.0
- Patient-level data bundled: none

## Verification completed

- All seven analysis definitions and guided user interfaces compiled with the
  jamovi compiler.
- `R CMD check --no-manual` completed with **Status: OK**.
- All 11 automated test cases passed across the protocol, DAG, sequential
  target-trial, clone-censor weighting, AIPW, landmark, and doubly robust
  survival workflows.
- The final `.jmo` was unpacked and loaded from its own bundled library using
  jamovi's R 4.6.0 runtime.
- A private SnapSBO compatibility run completed directly from the final `.jmo`:
  two landmarks, two NGT strata, AIPW and Hajek IPTW, eight finite estimates,
  and the intended NGT/WSEC display labels.
- The analysis-ready SnapSBO OMV was independently checked for 785 unique
  records and 37 analysis/audit variables. It is not part of the repository or
  release assets.

## Version 1.1.1 behaviours exercised

- Numbered analysis menu and guided eight-step landmark workflow.
- Eligibility-first binary validation, including extra coding levels confined
  to ineligible records.
- Labelled binary variables, custom treatment display names, and blank-as-
  absence recoding that preserves an existing negative label.
- Explicit stop/retain policies for missing cluster identifiers and missing,
  unparseable, or negative discharge times.
- Separate flow accounting for strategy classification, observation at the
  landmark, event-free status, and outcome-time usability.
- 24- and 48-hour landmark construction with centre-clustered inference.
- AIPW and Hajek IPTW estimates with automatic Firth fallback under separation.
- Missingness, timing, propensity overlap, effective sample size, weights, and
  covariate-balance diagnostics.

## Included analyses

1. Target-Trial Protocol
2. Causal DAG
3. Landmark Trial (24/48 h)
4. Sequential Target Trial
5. Clone-Censor Weighting
6. AIPW / Doubly Robust
7. Doubly Robust Survival

## Important limitations

- The supplied binary is for macOS Apple Silicon. Windows, Linux, and Intel
  macOS require platform-specific jamovi builds.
- A DAG must be prespecified from subject-matter knowledge; the module does not
  learn causal structure from observed associations.
- AIPW is doubly robust to one nuisance-model misspecification, not to violations
  of exchangeability, positivity, consistency, temporal alignment, measurement,
  or a well-defined intervention.
- The SnapSBO source lacks a verified NGT-placement time. Treating recorded NGT
  status as landmark stratum membership therefore reproduces the historical
  subgroup interpretation but does not establish 24/48-hour NGT timing.
- Missing and invalid dates should be investigated and prespecified sensitivity
  analyses reported; an audit flag is not a data repair.
- Hot-deck imputation is a transparent pragmatic sensitivity method, not a
  substitute for a protocol-specific imputation model in confirmatory work.

## Installation

Open jamovi, open the module library, choose **Sideload**, and select
`trialemulationj-1.1.1-macos-arm64.jmo`. The numbered analyses appear under
**Causal Inference**.
