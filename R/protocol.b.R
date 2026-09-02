#' @importFrom jmvcore .
# The compiler-generated no-data R wrapper still passes an object named
# `data` to the analysis class. Keep a private, non-exported sentinel so the
# wrapper also works outside the jamovi server, where no dataset is available.
data <- data.frame(.protocol = 1L)

protocolClass <- R6::R6Class(
    "protocolClass",
    inherit = protocolBase,
    private = list(
        .run = function() {
            labels <- list(
                studyTitle = "Study question",
                population = "Target population and setting",
                eligibility = "Eligibility criteria",
                strategies = "Treatment strategies",
                assignment = "Assignment procedure",
                timeZero = "Time zero and alignment",
                gracePeriod = "Grace/decision window",
                followUp = "Follow-up and censoring",
                outcome = "Outcome and horizon",
                competingEvents = "Competing events",
                causalContrast = "Causal contrast",
                estimand = "Primary estimand",
                confoundingPlan = "Confounding-control plan",
                analysisPlan = "Analysis and variance plan"
            )
            contrast <- switch(
                self$options$causalContrast,
                initiation = "Observational analogue of an assignment/initiation effect",
                perProtocol = "Per-protocol effect under adherence to the specified strategies",
                intentionToTreat = "Intention-to-treat effect based on observed randomized assignment"
            )
            estimand <- switch(
                self$options$estimand,
                ateDifference = "Population-average risk or mean difference",
                ateRatio = "Population-average risk ratio",
                survivalDifference = "Marginal survival probability difference",
                rmstDifference = "Restricted mean survival time difference"
            )
            values <- list(
                studyTitle = self$options$studyTitle,
                population = self$options$population,
                eligibility = self$options$eligibility,
                strategies = self$options$strategies,
                assignment = self$options$assignment,
                timeZero = self$options$timeZero,
                gracePeriod = self$options$gracePeriod,
                followUp = self$options$followUp,
                outcome = self$options$outcome,
                competingEvents = self$options$competingEvents,
                causalContrast = contrast,
                estimand = estimand,
                confoundingPlan = self$options$confoundingPlan,
                analysisPlan = self$options$analysisPlan
            )
            for (i in seq_along(values)) {
                value <- values[[i]]
                if (!nzchar(trimws(value))) value <- "Not specified"
                self$results$specification$addRow(
                    rowKey = paste0("component", i),
                    values = list(component = labels[[names(values)[[i]]]], specification = value)
                )
            }

            checks <- list(
                list("Question, population, and strategies are defined", all(nzchar(trimws(c(self$options$studyTitle, self$options$population, self$options$strategies)))), "Define the decision and all contrasted interventions."),
                list("Eligibility is defined at time zero", nzchar(trimws(self$options$eligibility)) && nzchar(trimws(self$options$timeZero)), "State when every eligibility criterion is assessed."),
                list("Eligibility, assignment, and follow-up are aligned", self$options$acknowledgesAlignment, "Do not start follow-up before or after strategy assignment."),
                list("All strategies are feasible for every eligible person", self$options$acknowledgesFeasible, "Restrict the target population if structural positivity fails."),
                list("Outcome, horizon, and competing events are defined", nzchar(trimws(self$options$outcome)) && nzchar(trimws(self$options$followUp)), "Prespecify the outcome clock, horizon, censoring, and competing-event handling."),
                list("Adjustment variables are baseline causes", self$options$acknowledgesBaseline && nzchar(trimws(self$options$confoundingPlan)), "Use subject-matter knowledge/DAGs; exclude descendants, colliders, and outcome-informed variables."),
                list("Structural and empirical positivity are separated", self$options$acknowledgesPositivity, "Diagnose feasibility separately from sparse-data/extreme-weight problems."),
                list("Estimator and dependence-aware inference are prespecified", nzchar(trimws(self$options$analysisPlan)), "Name the estimator, nuisance models, cluster level, bootstrap unit, and confidence interval."),
                list("Missingness and sensitivity analyses are prespecified", nzchar(trimws(self$options$missingPlan)) && nzchar(trimws(self$options$sensitivityPlan)), "State recoding, imputation, overlap, model, and timing sensitivity analyses.")
            )
            for (i in seq_along(checks)) {
                item <- checks[[i]]
                self$results$readiness$addRow(rowKey = paste0("check", i), values = list(
                    check = item[[1]],
                    status = if (isTRUE(item[[2]])) "Ready" else "Incomplete",
                    action = if (isTRUE(item[[2]])) "None" else item[[3]]
                ))
            }

            record <- c(values, list(
                missingPlan = self$options$missingPlan,
                sensitivityPlan = self$options$sensitivityPlan,
                acknowledgements = c(
                    alignment = self$options$acknowledgesAlignment,
                    feasible = self$options$acknowledgesFeasible,
                    baseline = self$options$acknowledgesBaseline,
                    positivity = self$options$acknowledgesPositivity
                )
            ))
            fingerprint <- .te_specification_fingerprint(record)
            self$results$report$setContent(private$.reportText(values, fingerprint))
            self$results$fingerprint$setContent(paste0(
                "Module: TrialEmulationj 1.1.1\n",
                "Specification SHA-256: ", fingerprint, "\n",
                "Short fingerprint: ", substr(fingerprint, 1, 16), "\n",
                "The fingerprint identifies the entered protocol text and choices; it does not hash the dataset."
            ))
        },

        .reportText = function(values, fingerprint) {
            clean <- function(x) if (nzchar(trimws(x))) x else "[Not specified]"
            paste(c(
                paste0("TARGET-TRIAL PROTOCOL: ", clean(values$studyTitle)),
                "",
                paste0("Target population and setting: ", clean(values$population)),
                paste0("Eligibility criteria: ", clean(values$eligibility)),
                paste0("Treatment strategies: ", clean(values$strategies)),
                paste0("Assignment procedure in the target trial: ", clean(values$assignment)),
                paste0("Time zero and alignment: ", clean(values$timeZero)),
                paste0("Grace/decision window: ", clean(values$gracePeriod)),
                paste0("Follow-up and censoring: ", clean(values$followUp)),
                paste0("Outcome and horizon: ", clean(values$outcome)),
                paste0("Competing events: ", clean(values$competingEvents)),
                paste0("Causal contrast: ", values$causalContrast),
                paste0("Primary estimand: ", values$estimand),
                paste0("Confounding-control plan: ", clean(values$confoundingPlan)),
                paste0("Analysis and variance plan: ", clean(values$analysisPlan)),
                paste0("Missing-data plan: ", clean(self$options$missingPlan)),
                paste0("Sensitivity analyses: ", clean(self$options$sensitivityPlan)),
                "",
                paste0("Specification fingerprint: ", fingerprint)
            ), collapse = "\n")
        }
    )
)
