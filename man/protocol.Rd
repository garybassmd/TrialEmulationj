\name{protocol}
\alias{protocol}
\alias{protocolBase}
\alias{protocolClass}
\alias{protocolOptions}
\alias{protocolResults}
\title{Target-Trial Protocol Builder}
\description{
Creates a structured target-trial specification, design-readiness checklist,
copy-ready protocol, and reproducibility fingerprint before outcome modelling.
}
\usage{
protocol(
  studyTitle = "", population = "", eligibility = "", strategies = "",
  assignment = "Random assignment to a treatment strategy at time zero",
  timeZero = "", gracePeriod = "", followUp = "", outcome = "",
  competingEvents = "", causalContrast = "initiation",
  estimand = "ateDifference", confoundingPlan = "", analysisPlan = "",
  missingPlan = "", sensitivityPlan = "", acknowledgesAlignment = FALSE,
  acknowledgesFeasible = FALSE, acknowledgesBaseline = FALSE,
  acknowledgesPositivity = FALSE, showReport = TRUE
)
}
\arguments{
  \item{studyTitle, population, eligibility, strategies}{Core question, population, eligibility, and strategy specifications.}
  \item{assignment, timeZero, gracePeriod, followUp}{Assignment and analysis-time specifications.}
  \item{outcome, competingEvents}{Outcome, horizon, and competing-event specifications.}
  \item{causalContrast}{Initiation/assignment analogue, per-protocol, or observed randomized intention-to-treat contrast.}
  \item{estimand}{Risk/mean difference, risk ratio, survival difference, or RMST difference.}
  \item{confoundingPlan, analysisPlan, missingPlan, sensitivityPlan}{Prespecified analysis plans.}
  \item{acknowledgesAlignment, acknowledgesFeasible, acknowledgesBaseline, acknowledgesPositivity}{Design-readiness acknowledgements.}
  \item{showReport}{Whether to show the copy-ready protocol.}
}
\value{A jamovi results object containing the specification, readiness checks, protocol text, and SHA-256 fingerprint.}
\details{This analysis does not use patient-level data and should be completed before inspecting treatment-effect estimates.}
