\name{targetTrial}
\alias{targetTrial}
\alias{targetTrialBase}
\alias{targetTrialClass}
\alias{targetTrialOptions}
\alias{targetTrialResults}
\title{Target Trial Emulation}
\description{
Runs a sequential target-trial emulation from longitudinal person-period data
using the \pkg{TrialEmulation} statistical engine. The analysis supports
marginal intention-to-treat and per-protocol estimands, inverse-probability
weights, a pooled logistic marginal structural model, and marginal cumulative
incidence or survival contrasts.
}
\usage{
targetTrial(
  data, id, period, treatment, outcome, eligible,
  estimand = "itt", outcomeCovs = NULL,
  useCensor = FALSE, censor = NULL,
  censorNumCovs = NULL, censorDenCovs = NULL,
  poolModels = "auto", switchNumCovs = NULL,
  switchDenCovs = NULL, followupModel = "quadratic",
  trialPeriodModel = "linear", targetPopulation = "all",
  targetTrialPeriod = 0, maxFollowUp = 10,
  confInt = TRUE, ciSamples = 500,
  predictionType = "cum_inc", truncateWeights = TRUE,
  truncatePercentile = 99, pControl = 1, chunkSize = 500,
  seed = 20240219, showSyntax = TRUE, showWeightPlot = TRUE
)
}
\arguments{
  \item{data}{A data frame in longitudinal person-period format.}
  \item{id}{Name of the individual identifier variable.}
  \item{period}{Name of the non-negative integer period variable.}
  \item{treatment}{Name of the binary treatment variable.}
  \item{outcome}{Name of the binary outcome variable.}
  \item{eligible}{Name of the binary trial-eligibility variable.}
  \item{estimand}{Either \code{"itt"} or \code{"pp"}.}
  \item{outcomeCovs}{Baseline covariates for the outcome model.}
  \item{useCensor}{Whether to estimate informative-censoring weights.}
  \item{censor}{Name of the binary censoring-event variable.}
  \item{censorNumCovs, censorDenCovs}{Covariates for the numerator and denominator censoring models.}
  \item{poolModels}{Pooling rule for censoring models: \code{"auto"}, \code{"none"}, \code{"numerator"}, or \code{"both"}.}
  \item{switchNumCovs, switchDenCovs}{Covariates for per-protocol adherence-weight models.}
  \item{followupModel}{Follow-up-time specification: \code{"linear"}, \code{"quadratic"}, or \code{"spline3"}.}
  \item{trialPeriodModel}{Trial-period specification: \code{"none"}, \code{"linear"}, or \code{"quadratic"}.}
  \item{targetPopulation}{Baseline population: \code{"all"}, \code{"earliest"}, \code{"latest"}, or \code{"specific"}.}
  \item{targetTrialPeriod}{Trial period used when \code{targetPopulation = "specific"}.}
  \item{maxFollowUp}{Maximum follow-up period for marginal predictions.}
  \item{confInt}{Whether to calculate pointwise 95 percent confidence intervals.}
  \item{ciSamples}{Number of simulation draws for confidence intervals.}
  \item{predictionType}{Either \code{"cum_inc"} or \code{"survival"}.}
  \item{truncateWeights}{Whether to cap final analysis weights.}
  \item{truncatePercentile}{Upper percentile used as the weight cap.}
  \item{pControl}{Sampling probability for non-event person-period rows.}
  \item{chunkSize}{Number of persons processed per expansion chunk.}
  \item{seed}{Random seed for sampling and simulation-based intervals.}
  \item{showSyntax}{Whether to show reproducible R syntax.}
  \item{showWeightPlot}{Whether to show the final-weight histogram.}
}
\value{
A jamovi results object containing summaries, model coefficients, marginal
predictions, treatment contrasts, diagnostics, plots, assumptions, and R code.
}
\details{
Causal interpretation requires a prespecified target-trial protocol and the
identification assumptions described by the TrialEmulation methodology.
}
