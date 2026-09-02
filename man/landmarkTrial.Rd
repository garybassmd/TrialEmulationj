\name{landmarkTrial}
\alias{landmarkTrial}
\alias{landmarkTrialBase}
\alias{landmarkTrialClass}
\alias{landmarkTrialOptions}
\alias{landmarkTrialResults}
\title{Landmark Target-Trial Analysis}
\description{Constructs audited landmark cohorts for one or two decision windows and estimates within-stratum AIPW and IPTW effects.}
\usage{
landmarkTrial(
  data, eligibility = NULL, treatmentA, treatmentB, outcome = NULL,
  failureEvent1 = NULL, failureEvent2 = NULL, competingEvent = NULL,
  cluster = NULL, indexTime = NULL, treatmentATime = NULL,
  treatmentBTime = NULL, failureTime1 = NULL, failureTime2 = NULL,
  dischargeTime = NULL, competingTime = NULL, covariates = NULL,
  absenceVars = NULL, eligibilityPositive = "", treatmentAPositive = "",
  treatmentBPositive = "", outcomePositive = "", outcomeSource = "recorded",
  competingHandling = "noFailure", timeMode = "elapsed", timeUnit = "hours",
  window1 = 24, runSecondWindow = TRUE, window2 = 48,
  unknownTreatment = "exclude", unknownEvent = "exclude",
  runBothStrata = TRUE, primaryStratum = "present", primaryWindow = "first",
  missingMethod = "fail", imputations = 20, outcomeMethod = "auto",
  propensityMethod = "auto", minPS = 0.01, maxPS = 0.99,
  bootstrapSamples = 0, randomSeed = 24680, confidenceLevel = 95,
  effectMeasure = "riskDifference", showFlowPlot = TRUE,
  showBalancePlot = TRUE, showEffectPlot = TRUE, showSyntax = TRUE
)
}
\arguments{
  \item{data}{One row per person before the landmark cohort is constructed.}
  \item{eligibility}{Optional baseline eligibility indicator.}
  \item{treatmentA, treatmentB}{Binary treatment components; A defines strata and B is compared within strata.}
  \item{outcome}{Recorded binary failure outcome.}
  \item{failureEvent1, failureEvent2}{Optional binary components for a composite outcome and timing checks.}
  \item{competingEvent}{Optional competing-event indicator.}
  \item{cluster}{Optional cluster identifier.}
  \item{indexTime}{Index date-time when \code{timeMode = "dateTime"}.}
  \item{treatmentATime, treatmentBTime, failureTime1, failureTime2, dischargeTime, competingTime}{Elapsed or date-time variables used to align each landmark.}
  \item{covariates}{Baseline adjustment covariates.}
  \item{absenceVars}{Variables explicitly recoded from missing/blank to absence.}
  \item{eligibilityPositive, treatmentAPositive, treatmentBPositive, outcomePositive}{Optional positive-level labels.}
  \item{outcomeSource}{Use the recorded outcome or combine selected failure components.}
  \item{competingHandling}{Treat as no failure, include in a composite, or exclude as a sensitivity analysis.}
  \item{timeMode, timeUnit}{Elapsed/date-time representation and hours/days scale.}
  \item{window1, runSecondWindow, window2}{Landmark decision-window settings.}
  \item{unknownTreatment, unknownEvent}{Policies for recorded treatment/events with unusable timing.}
  \item{runBothStrata, primaryStratum, primaryWindow}{Stratum and primary diagnostic-plot settings.}
  \item{missingMethod, imputations}{Baseline-covariate missingness policy.}
  \item{outcomeMethod, propensityMethod}{Automatic, maximum-likelihood, or Firth binary models.}
  \item{minPS, maxPS}{Propensity-score truncation limits.}
  \item{bootstrapSamples, randomSeed, confidenceLevel}{Inference settings.}
  \item{effectMeasure}{Risk difference, risk ratio, or odds ratio.}
  \item{showFlowPlot, showBalancePlot, showEffectPlot, showSyntax}{Output display options.}
}
\value{A jamovi results object containing cohort flow, timing audit, four-strategy counts, effects, balance, plots, diagnostics, and reproducible syntax.}
\details{Eligibility, treatment classification, and follow-up are aligned at each landmark. Causal interpretation requires a prespecified protocol, exchangeability, positivity, consistency, and baseline-only adjustment.}
