\name{causalDAG}
\alias{causalDAG}
\alias{causalDAGBase}
\alias{causalDAGClass}
\alias{causalDAGOptions}
\alias{causalDAGResults}
\title{Causal DAG and Adjustment Sets}
\description{Validates a prespecified DAG, identifies minimally sufficient adjustment sets, checks selected covariates, and draws the graph.}
\usage{
causalDAG(
  data, exposure, outcome, selectedAdjusters = NULL, dagText = "",
  effect = "total", maxSets = 20, showIndependencies = TRUE,
  showPlot = TRUE, showSyntax = TRUE
)
}
\arguments{
  \item{data}{A data frame whose selected column names correspond to DAG nodes.}
  \item{exposure, outcome}{Exposure and outcome column/node names.}
  \item{selectedAdjusters}{Proposed baseline adjustment variables.}
  \item{dagText}{DAGitty syntax; a bare list of arrows is wrapped in \code{dag \{ \}}.}
  \item{effect}{Either \code{"total"} or \code{"direct"}.}
  \item{maxSets}{Maximum minimal adjustment sets to return.}
  \item{showIndependencies, showPlot, showSyntax}{Output display options.}
}
\value{A jamovi results object containing graph checks, adjustment sets, selected-variable review, implied independencies, a DAG plot, and validated syntax.}
\details{The module does not learn a causal structure from data. Valid adjustment depends on the causal graph being correct and required variables being measured adequately.}
