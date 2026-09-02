\name{drSurvival}
\alias{drSurvival}
\alias{drSurvivalBase}
\alias{drSurvivalClass}
\alias{drSurvivalOptions}
\alias{drSurvivalResults}
\title{Doubly Robust Survival Effect}
\description{
Estimates a doubly robust survival treatment effect with outcome-regression,
inverse-probability treatment, and inverse-probability censoring components.
The primary estimand is a restricted-mean-time-lost ratio.
}
\usage{
drSurvival(
  data, time, event, treatment, followupTime = NULL,
  outcomeCovs = NULL, propensityCovs = NULL, censorCovs = NULL,
  horizon = 365, psMethod = "glm", censorMethod = "breslow",
  minPS = 0.01, maxPS = 0.99, minCensorSurvival = 0.025,
  bootstrapSamples = 500, seed = 20240219,
  showBootstrapPlot = TRUE, showSyntax = TRUE
)
}
\arguments{
  \item{data}{A survival data frame with one row per person.}
  \item{time}{Name of observed follow-up time.}
  \item{event}{Name of the binary event indicator.}
  \item{treatment}{Name of the binary treatment variable coded 0/1.}
  \item{followupTime}{Optional name of potential censoring time.}
  \item{outcomeCovs}{Covariates in the survival outcome model.}
  \item{propensityCovs}{Covariates in the treatment propensity model.}
  \item{censorCovs}{Covariates in the censoring model.}
  \item{horizon}{Truncation time defining restricted mean time lost.}
  \item{psMethod}{Either \code{"glm"} or \code{"lasso"}.}
  \item{censorMethod}{One of \code{"breslow"}, \code{"aft_exp"}, \code{"aft_weibull"}, \code{"aft_lognormal"}, or \code{"aft_loglogistic"}.}
  \item{minPS, maxPS}{Propensity-score truncation limits.}
  \item{minCensorSurvival}{Lower bound for estimated censoring survival.}
  \item{bootstrapSamples}{Number of bootstrap samples.}
  \item{seed}{Random seed for bootstrap inference.}
  \item{showBootstrapPlot}{Whether to show bootstrap effect distributions.}
  \item{showSyntax}{Whether to show reproducible R syntax.}
}
\value{
A jamovi results object with treatment-specific RMST, the RMTL ratio, an
adjusted hazard ratio, bootstrap inference, diagnostics, assumptions, and R
syntax.
}
\details{
The estimator delegates to \code{precmed::atefitsurv}. Causal interpretation
requires the usual consistency, positivity, no-interference, and exchangeability
assumptions in addition to adequate nuisance-model specification.
}
