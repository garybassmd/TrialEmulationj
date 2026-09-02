\name{cloneCensor}
\alias{cloneCensor}
\alias{cloneCensorBase}
\alias{cloneCensorClass}
\alias{cloneCensorOptions}
\alias{cloneCensorResults}
\title{Clone-Censor-Weight Survival Analysis}
\description{
Clones each person into two grace-period exposure strategies, artificially
censors clones when they deviate, estimates inverse-probability-of-censoring
weights, and reports weighted survival effects with person-level bootstrap
inference.
}
\usage{
cloneCensor(
  data, id, event, timeToEvent, exposure, timeToExposure,
  weightCovs = NULL, gracePeriod = 100, horizon = 365,
  winsorize = TRUE, lowerPercentile = 1, upperPercentile = 99,
  winsorizeByClone = FALSE, bootstrapSamples = 500,
  seed = 20240219, showWeightPlot = TRUE, showSyntax = TRUE
)
}
\arguments{
  \item{data}{A baseline data frame with one row per person.}
  \item{id}{Name of the unique person identifier.}
  \item{event}{Name of the binary outcome-event variable.}
  \item{timeToEvent}{Name of time to the outcome event or censoring.}
  \item{exposure}{Name of exposure during the grace period, coded 0/1.}
  \item{timeToExposure}{Name of time to exposure; may be missing for unexposed people.}
  \item{weightCovs}{Numeric baseline predictors of artificial censoring.}
  \item{gracePeriod}{Time at which the exposure-assignment window closes.}
  \item{horizon}{Time horizon for survival and RMST contrasts.}
  \item{winsorize}{Whether to winsorize clone-censor weights.}
  \item{lowerPercentile, upperPercentile}{Weight winsorization percentiles.}
  \item{winsorizeByClone}{Whether to winsorize separately by assigned strategy.}
  \item{bootstrapSamples}{Number of person-level bootstrap samples.}
  \item{seed}{Random seed for bootstrap resampling.}
  \item{showWeightPlot}{Whether to show the clone-censor-weight histogram.}
  \item{showSyntax}{Whether to show reproducible R syntax.}
}
\value{
A jamovi results object containing weighted hazard, survival, and RMST
contrasts, percentile bootstrap intervals, survival curves, weight diagnostics,
assumptions, and R syntax.
}
\details{
The workflow delegates cloning and weighting to \pkg{survivalCCW}. People are
resampled before cloning so bootstrap inference preserves within-person
dependence between clones.
}
