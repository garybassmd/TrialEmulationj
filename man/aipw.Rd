\name{aipw}
\alias{aipw}
\alias{aipwBase}
\alias{aipwClass}
\alias{aipwOptions}
\alias{aipwResults}
\title{AIPW and IPTW for Point Outcomes}
\description{
Estimates marginal effects for a binary or continuous point outcome using
augmented inverse probability weighting and Hajek-normalized inverse
probability weighting, with cluster-aware inference and diagnostics.
}
\usage{
aipw(
  data, treatment, outcome, treatmentPositive = "", outcomePositive = "",
  outcomeType = "binary", outcomeCovs = NULL, propensityCovs = NULL,
  cluster = NULL, absenceVars = NULL, missingMethod = "fail",
  imputations = 20, randomSeed = 24680, outcomeSpecification = "pooled",
  outcomeMethod = "auto", propensityMethod = "auto",
  minPS = 0.01, maxPS = 0.99, bootstrapSamples = 0,
  confidenceLevel = 95, plotMeasure = "difference",
  sensitivityGrid = TRUE, showOverlapPlot = TRUE,
  showBalancePlot = TRUE, showEffectPlot = TRUE, showSyntax = TRUE
)
}
\arguments{
  \item{data}{A data frame with one independent row per observational unit.}
  \item{treatment}{Name of the binary treatment variable.}
  \item{outcome}{Name of the point-outcome variable.}
  \item{treatmentPositive, outcomePositive}{Optional labels identifying the treated and event levels.}
  \item{outcomeType}{Either \code{"binary"} or \code{"continuous"}.}
  \item{outcomeCovs}{Covariates in both arm-specific outcome models.}
  \item{propensityCovs}{Covariates in the treatment propensity model.}
  \item{cluster}{Optional cluster identifier for cluster-robust inference and whole-cluster bootstrap sampling.}
  \item{absenceVars}{Variables for which missing/blank is explicitly recoded as absence.}
  \item{missingMethod}{One of \code{"fail"}, \code{"completeCases"}, or \code{"hotDeckMI"}.}
  \item{imputations}{Number of hot-deck imputations.}
  \item{randomSeed}{Random seed for imputation and bootstrap sampling.}
  \item{outcomeSpecification}{Either a \code{"pooled"} or \code{"armSpecific"} outcome model.}
  \item{outcomeMethod, propensityMethod}{One of \code{"auto"}, \code{"maximumLikelihood"}, or \code{"firth"}.}
  \item{minPS, maxPS}{Lower and upper truncation bounds for estimated propensity scores.}
  \item{bootstrapSamples}{Number of row or whole-cluster bootstrap samples; zero disables bootstrap intervals.}
  \item{confidenceLevel}{Confidence level as a percentage.}
  \item{plotMeasure}{Difference, risk ratio, or odds ratio for the effect plot.}
  \item{sensitivityGrid}{Whether to refit a grid of propensity-score truncation limits.}
  \item{showOverlapPlot}{Whether to show propensity-score overlap.}
  \item{showBalancePlot, showEffectPlot}{Whether to show balance and effect plots.}
  \item{showSyntax}{Whether to show reproducible R syntax.}
}
\value{
A jamovi results object with marginal potential-outcome means, the ATE,
influence-function standard errors, overlap and weight diagnostics, balance,
assumptions, plots, and R syntax.
}
\details{
The estimator is doubly robust with respect to the propensity and outcome
nuisance models. Causal interpretation still requires consistency, positivity,
no interference, and conditional exchangeability.
}
