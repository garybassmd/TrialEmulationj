#' @importFrom jmvcore .
aipwClass <- R6::R6Class(
    "aipwClass",
    inherit = aipwBase,
    private = list(
        .run = function() {
            if (is.null(self$options$treatment) || is.null(self$options$outcome))
                return()
            fitted <- tryCatch(private$.fitAnalysis(), error = function(e) e)
            if (inherits(fitted, "error")) {
                jmvcore::reject(paste0(
                    "The AIPW/IPTW analysis could not be completed: ",
                    conditionMessage(fitted),
                    " Check coding, missing-data choices, model variables, and treatment overlap."
                ))
            }
            private$.populateResults(fitted)
        },

        .prepareData = function() {
            treatmentName <- self$options$treatment
            outcomeName <- self$options$outcome
            if (identical(treatmentName, outcomeName))
                stop("Treatment and outcome must be different variables.")
            clusterName <- self$options$cluster
            covariates <- unique(c(self$options$outcomeCovs, self$options$propensityCovs))
            covariates <- setdiff(covariates, c(treatmentName, outcomeName, clusterName))
            used <- unique(c(treatmentName, outcomeName, clusterName, covariates, self$options$absenceVars))
            used <- used[!is.na(used) & nzchar(used)]
            source <- as.data.frame(self$data[, used, drop = FALSE])
            originalN <- nrow(source)
            if (originalN < 20L)
                stop("At least 20 observations are required.")

            recoded <- .te_recode_absence(source, self$options$absenceVars)
            source <- recoded$data
            audit <- recoded$audit
            treatmentRaw <- .te_trim(source[[treatmentName]])
            outcomeRaw <- .te_trim(source[[outcomeName]])
            roleMissing <- is.na(treatmentRaw) | is.na(outcomeRaw)
            if (any(roleMissing)) {
                if (identical(self$options$missingMethod, "fail"))
                    stop(sum(roleMissing), " observations have missing treatment or outcome values.")
                audit <- rbind(audit, data.frame(
                    variable = "Treatment or outcome",
                    recoded = sum(roleMissing),
                    rule = "Rows excluded; treatment and outcome are not imputed",
                    stringsAsFactors = FALSE
                ))
                source <- source[!roleMissing, , drop = FALSE]
            }
            if (!is.null(clusterName) && anyNA(.te_trim(source[[clusterName]])))
                stop("Cluster ID contains missing values. Assign every included row to a cluster.")

            treatment <- .te_as_binary(source[[treatmentName]], "Treatment", self$options$treatmentPositive)
            if (identical(self$options$outcomeType, "binary")) {
                outcome <- .te_as_binary(source[[outcomeName]], "Binary outcome", self$options$outcomePositive)
                Y <- outcome$values
                outcomeMapping <- paste0(outcome$zero, " = 0; ", outcome$one, " = 1")
            } else {
                Y <- suppressWarnings(as.numeric(.te_trim(source[[outcomeName]])))
                if (any(!is.finite(Y)))
                    stop("The continuous outcome must contain finite numeric values.")
                outcomeMapping <- "Continuous numeric outcome"
            }
            A <- treatment$values

            internal <- if (length(covariates) > 0L) sprintf("x%03d", seq_along(covariates)) else character()
            names(internal) <- covariates
            analysis <- data.frame(A = A, Y = Y)
            for (name in covariates)
                analysis[[internal[[name]]]] <- source[[name]]
            if (!is.null(clusterName))
                analysis$.cluster <- as.character(source[[clusterName]])

            covInternal <- unname(internal)
            covMissing <- if (length(covInternal) == 0L) rep(FALSE, nrow(analysis)) else !stats::complete.cases(analysis[, covInternal, drop = FALSE])
            missingCounts <- if (length(covInternal) == 0L) integer() else vapply(analysis[, covInternal, drop = FALSE], function(x) sum(is.na(.te_trim(x))), integer(1))
            if (any(covMissing)) {
                if (identical(self$options$missingMethod, "fail")) {
                    names(missingCounts) <- covariates
                    stop("Covariates contain missing values: ", paste(names(missingCounts)[missingCounts > 0], collapse = ", "), ".")
                }
                if (identical(self$options$missingMethod, "completeCases")) {
                    audit <- rbind(audit, data.frame(
                        variable = "Model covariates",
                        recoded = sum(covMissing),
                        rule = "Rows excluded by complete-case analysis",
                        stringsAsFactors = FALSE
                    ))
                    analysis <- analysis[!covMissing, , drop = FALSE]
                } else {
                    audit <- rbind(audit, data.frame(
                        variable = "Model covariates",
                        recoded = sum(covMissing),
                        rule = paste0(self$options$imputations, " stochastic hot-deck imputations within treatment groups"),
                        stringsAsFactors = FALSE
                    ))
                }
            }

            if (nrow(analysis) < 20L)
                stop("Fewer than 20 observations remain after missing-data handling.")
            if (length(unique(analysis$A)) != 2L || min(table(analysis$A)) < 5L)
                stop("Each treatment group must contain at least five observations.")

            labelMap <- setNames(covariates, unname(internal))
            list(
                data = analysis,
                originalN = originalN,
                audit = audit,
                outcomeVars = unname(internal[intersect(self$options$outcomeCovs, covariates)]),
                propensityVars = unname(internal[intersect(self$options$propensityCovs, covariates)]),
                imputeVars = covInternal,
                labelMap = labelMap,
                treatmentMapping = paste0(treatment$zero, " = 0; ", treatment$one, " = 1"),
                outcomeMapping = outcomeMapping,
                clusterName = if (is.null(clusterName)) NULL else ".cluster",
                originalTreatment = treatmentName,
                originalOutcome = outcomeName,
                originalCluster = clusterName,
                originalOutcomeCovs = self$options$outcomeCovs,
                originalPropensityCovs = self$options$propensityCovs
            )
        },

        .fitArguments = function(prepared) {
            list(
                outcomeType = self$options$outcomeType,
                outcomeVars = prepared$outcomeVars,
                propensityVars = prepared$propensityVars,
                minPS = self$options$minPS,
                maxPS = self$options$maxPS,
                outcomeMethod = self$options$outcomeMethod,
                propensityMethod = self$options$propensityMethod,
                outcomeSpecification = self$options$outcomeSpecification,
                confidenceLevel = self$options$confidenceLevel
            )
        },

        .fitAnalysis = function() {
            if (self$options$minPS >= self$options$maxPS)
                stop("The minimum propensity score must be below the maximum propensity score.")
            prepared <- private$.prepareData()
            fitArgs <- private$.fitArguments(prepared)
            cluster <- if (is.null(prepared$clusterName)) NULL else prepared$data[[prepared$clusterName]]

            if (identical(self$options$missingMethod, "hotDeckMI") && anyNA(prepared$data[, prepared$imputeVars, drop = FALSE])) {
                if (self$options$bootstrapSamples > 0L)
                    stop("Bootstrap confidence intervals and multiple hot-deck imputation cannot currently be combined. Use clustered influence-function inference or complete-case analysis.")
                imputations <- .te_hotdeck_imputations(
                    prepared$data,
                    prepared$imputeVars,
                    self$options$imputations,
                    self$options$randomSeed,
                    strata = prepared$data$A
                )
                fits <- lapply(imputations, function(item) {
                    args <- c(list(data = item, cluster = if (is.null(prepared$clusterName)) NULL else item[[prepared$clusterName]]), fitArgs)
                    do.call(.te_fit_causal_once, args)
                })
                fitted <- fits[[1]]
                fitted$effects <- .te_pool_mi_effects(fits, self$options$confidenceLevel)
                fitted$miFits <- length(fits)
                analysisData <- imputations[[1]]
            } else {
                args <- c(list(data = prepared$data, cluster = cluster), fitArgs)
                fitted <- do.call(.te_fit_causal_once, args)
                draws <- .te_bootstrap_causal(
                    prepared$data,
                    fitArgs,
                    self$options$bootstrapSamples,
                    self$options$randomSeed,
                    prepared$clusterName
                )
                fitted$effects <- .te_apply_bootstrap_ci(fitted$effects, draws, self$options$confidenceLevel)
                fitted$miFits <- 0L
                analysisData <- prepared$data
            }
            fitted$balance$covariate <- private$.restoreLabels(fitted$balance$covariate, prepared$labelMap)
            fitted$prepared <- prepared
            fitted$analysisData <- analysisData
            fitted$sensitivity <- private$.sensitivityGrid(analysisData, prepared, fitArgs)
            fitted$fingerprint <- .te_specification_fingerprint(list(
                analysis = "aipw",
                version = "1.1.1",
                treatment = prepared$originalTreatment,
                outcome = prepared$originalOutcome,
                cluster = prepared$originalCluster,
                outcomeCovariates = prepared$originalOutcomeCovs,
                propensityCovariates = prepared$originalPropensityCovs,
                options = list(
                    treatmentPositive = self$options$treatmentPositive,
                    outcomePositive = self$options$outcomePositive,
                    outcomeType = self$options$outcomeType,
                    absenceVariables = self$options$absenceVars,
                    missingMethod = self$options$missingMethod,
                    imputations = self$options$imputations,
                    randomSeed = self$options$randomSeed,
                    outcomeSpecification = self$options$outcomeSpecification,
                    outcomeMethod = self$options$outcomeMethod,
                    propensityMethod = self$options$propensityMethod,
                    minPS = self$options$minPS,
                    maxPS = self$options$maxPS,
                    bootstrapSamples = self$options$bootstrapSamples,
                    confidenceLevel = self$options$confidenceLevel,
                    sensitivityGrid = self$options$sensitivityGrid
                )
            ))
            fitted
        },

        .restoreLabels = function(labels, mapping) {
            for (internal in names(mapping))
                labels <- sub(paste0("^", internal), mapping[[internal]], labels)
            labels
        },

        .sensitivityGrid = function(data, prepared, fitArgs) {
            if (!isTRUE(self$options$sensitivityGrid))
                return(data.frame())
            limits <- unique(c(.001, .01, .025, .05, .10, self$options$minPS))
            limits <- sort(limits[limits < .5])
            rows <- list()
            cluster <- if (is.null(prepared$clusterName)) NULL else data[[prepared$clusterName]]
            for (limit in limits) {
                args <- fitArgs
                args$minPS <- limit
                args$maxPS <- 1 - limit
                fit <- tryCatch(do.call(.te_fit_causal_once, c(list(data = data, cluster = cluster), args)), error = function(e) NULL)
                if (is.null(fit))
                    next
                difference <- if (identical(self$options$outcomeType, "binary")) "Risk difference" else "Mean difference"
                keep <- fit$effects$estimand == difference
                selected <- fit$effects[keep, , drop = FALSE]
                for (i in seq_len(nrow(selected))) {
                    rows[[length(rows) + 1L]] <- data.frame(
                        limits = paste0("[", format(limit, trim = TRUE), ", ", format(1 - limit, trim = TRUE), "]"),
                        method = selected$method[[i]],
                        estimate = selected$estimate[[i]],
                        lower = selected$lower[[i]],
                        upper = selected$upper[[i]],
                        truncated = fit$truncated,
                        stringsAsFactors = FALSE
                    )
                }
            }
            if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
        },

        .populateResults = function(fitted) {
            prepared <- fitted$prepared
            n <- length(fitted$A)
            maxBalance <- if (nrow(fitted$balance) == 0L) NA_real_ else max(abs(fitted$balance$weighted), na.rm = TRUE)
            clusterText <- if (is.null(prepared$clusterName)) "Independent observations" else paste0(length(unique(fitted$analysisData[[prepared$clusterName]])), " clusters; cluster-robust SE")
            missingText <- switch(
                self$options$missingMethod,
                fail = "No remaining missing values allowed",
                completeCases = "Complete-case analysis",
                hotDeckMI = paste0(fitted$miFits, " hot-deck imputations")
            )
            rows <- list(
                c("Rows received after jamovi filters", format(prepared$originalN, big.mark = ",")),
                c("Rows analyzed", format(n, big.mark = ",")),
                c("Rows excluded", format(prepared$originalN - n, big.mark = ",")),
                c("Treatment coding", prepared$treatmentMapping),
                c("Outcome coding", prepared$outcomeMapping),
                c("Untreated / treated", paste0(sum(fitted$A == 0), " / ", sum(fitted$A == 1))),
                c("Dependence", clusterText),
                c("Missing-data method", missingText),
                c("Propensity model", fitted$psMethod),
                c("Outcome model", fitted$outcomeMethod),
                c("Outcome specification", if (self$options$outcomeSpecification == "pooled") "Pooled model with treatment term" else "Separate arm-specific models"),
                c("Propensity estimates truncated", paste0(fitted$truncated, " (", round(100 * fitted$truncated / n, 1), "%)")),
                c("Largest weighted absolute SMD", if (is.finite(maxBalance)) format(round(maxBalance, 4), nsmall = 4) else "Not applicable"),
                c("Specification fingerprint", substr(fitted$fingerprint, 1, 16))
            )
            for (i in seq_along(rows))
                self$results$summary$setRow(rowNo = i, values = list(measure = rows[[i]][1], value = rows[[i]][2]))

            for (i in seq_len(nrow(fitted$effects)))
                self$results$effects$addRow(rowKey = paste0("effect", i), values = as.list(fitted$effects[i, ]))
            for (group in 0:1) {
                keep <- fitted$A == group
                self$results$propensity$addRow(
                    rowKey = paste0("group", group),
                    values = list(
                        group = if (group == 0) "No treatment" else "Treatment",
                        n = sum(keep),
                        psMin = min(fitted$rawPS[keep]),
                        psMedian = stats::median(fitted$rawPS[keep]),
                        psMax = max(fitted$rawPS[keep]),
                        weightMean = mean(fitted$ipw[keep]),
                        weightMax = max(fitted$ipw[keep]),
                        ess = fitted$ess[[if (group == 0) "untreated" else "treated"]]
                    )
                )
            }
            if (nrow(fitted$balance) > 0L) {
                for (i in seq_len(nrow(fitted$balance))) {
                    row <- fitted$balance[i, ]
                    self$results$balance$addRow(rowKey = paste0("balance", i), values = list(
                        covariate = row$covariate,
                        unweighted = row$unweighted,
                        weighted = row$weighted,
                        acceptable = if (abs(row$weighted) < .10) "Yes" else "No"
                    ))
                }
            }
            if (nrow(prepared$audit) > 0L) {
                for (i in seq_len(nrow(prepared$audit)))
                    self$results$missingAudit$addRow(rowKey = paste0("audit", i), values = as.list(prepared$audit[i, ]))
            }
            if (nrow(fitted$sensitivity) > 0L) {
                for (i in seq_len(nrow(fitted$sensitivity)))
                    self$results$sensitivity$addRow(rowKey = paste0("sensitivity", i), values = as.list(fitted$sensitivity[i, ]))
            }

            self$results$overlapPlot$setState(data.frame(
                propensity = fitted$rawPS,
                treatment = factor(fitted$A, levels = c(0, 1), labels = c("No treatment", "Treatment"))
            ))
            self$results$balancePlot$setState(fitted$balance)
            self$results$effectPlot$setState(private$.effectPlotData(fitted$effects))
            self$results$notes$setContent(private$.notes(fitted))
            self$results$syntax$setContent(private$.syntaxText(fitted))
        },

        .effectPlotData = function(effects) {
            desired <- switch(
                self$options$plotMeasure,
                difference = if (self$options$outcomeType == "binary") "Risk difference" else "Mean difference",
                riskRatio = "Risk ratio",
                oddsRatio = "Odds ratio"
            )
            data <- effects[effects$estimand == desired, c("method", "estimand", "estimate", "lower", "upper"), drop = FALSE]
            data$null <- if (desired %in% c("Risk ratio", "Odds ratio")) 1 else 0
            data
        },

        .notes = function(fitted) {
            notes <- c(
                "The AIPW estimate is doubly robust to nuisance-model misspecification: consistency requires a correct propensity model or a correct outcome model, not necessarily both. It still requires a well-defined intervention, conditional exchangeability, positivity, consistency, and appropriate temporal ordering.",
                "IPTW is shown as a model-comparison diagnostic and uses Hajek-normalized weights.",
                if (is.null(fitted$prepared$clusterName)) "Influence-function standard errors treat rows as independent." else "Influence-function standard errors aggregate contributions within the selected cluster. Bootstrap samples resample whole clusters and preserve all rows within a sampled cluster.",
                if (self$options$bootstrapSamples > 0L) paste0("Confidence limits for supported contrasts are percentile limits from ", self$options$bootstrapSamples, " refitted bootstrap samples; the SE column remains influence-function based.") else "Confidence limits use influence-function standard errors.",
                if (grepl("Firth", fitted$outcomeMethod) || grepl("Firth", fitted$psMethod)) "Firth mean bias-reduced logistic regression was used because it was requested or because separation/instability was detected." else "No Firth fallback was required.",
                "Absolute weighted SMD below 0.10 is a common balance diagnostic, not proof of exchangeability. Inspect overlap, effective sample size, extreme weights, and the truncation sensitivity table.",
                if (identical(self$options$missingMethod, "hotDeckMI")) "Hot-deck multiple imputation samples observed donors within treatment groups and combines estimates using Rubin-style within/between-imputation variance. For confirmatory analyses, compare with a protocol-specific model-based imputation analysis." else NULL,
                if (length(self$options$absenceVars) > 0L) "Blank-as-absence recoding was explicitly requested. The audit table reports every changed variable; this rule should only be used when supported by the data-collection process." else "No blank-as-absence recoding was requested.",
                "Only baseline covariates measured before treatment assignment should enter nuisance models. Do not adjust for mediators, colliders, post-baseline variables, or transformations that use the outcome.",
                "An observational comparison of treatment initiation is not literally an intention-to-treat effect unless randomized assignment is observed. State the observational analog and target population explicitly.",
                "Active jamovi filters are applied before the module receives the rows. Record those filters with the protocol and use the specification fingerprint to identify this setup.",
                if (length(fitted$warnings) > 0L) paste0("Model warning: ", paste(unique(fitted$warnings), collapse = "; ")) else NULL
            )
            paste0("- ", notes, collapse = "\n")
        },

        .overlapPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || nrow(image$state) == 0L)
                return(FALSE)
            ggplot2::ggplot(image$state, ggplot2::aes(x = propensity, fill = treatment, colour = treatment)) +
                ggplot2::geom_density(alpha = .18, linewidth = .9) +
                ggplot2::geom_vline(xintercept = c(self$options$minPS, self$options$maxPS), linetype = "dashed", colour = "#666666") +
                ggplot2::scale_fill_manual(values = c("No treatment" = "#4477AA", "Treatment" = "#CC6677")) +
                ggplot2::scale_colour_manual(values = c("No treatment" = "#4477AA", "Treatment" = "#CC6677")) +
                ggplot2::labs(x = "Estimated propensity score", y = "Density", fill = "Group", colour = "Group") +
                ggtheme
        },

        .balancePlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || nrow(image$state) == 0L)
                return(FALSE)
            state <- image$state
            order <- state$covariate[order(abs(state$unweighted))]
            plotData <- rbind(
                data.frame(covariate = state$covariate, smd = state$unweighted, sample = "Unadjusted"),
                data.frame(covariate = state$covariate, smd = state$weighted, sample = "Adjusted")
            )
            plotData$covariate <- factor(plotData$covariate, levels = unique(order))
            ggplot2::ggplot(plotData, ggplot2::aes(x = smd, y = covariate, colour = sample)) +
                ggplot2::geom_vline(xintercept = 0, colour = "#333333", linewidth = .45) +
                ggplot2::geom_vline(xintercept = c(-.10, .10), linetype = "dashed", colour = "#777777", linewidth = .45) +
                ggplot2::geom_point(size = 2.6) +
                ggplot2::scale_colour_manual(values = c("Unadjusted" = "#F8766D", "Adjusted" = "#00AEB3")) +
                ggplot2::labs(x = "Standardized mean difference", y = NULL, colour = "Sample") +
                ggtheme
        },

        .effectPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || nrow(image$state) == 0L)
                return(FALSE)
            data <- image$state
            ggplot2::ggplot(data, ggplot2::aes(y = method, x = estimate, colour = method)) +
                ggplot2::geom_vline(xintercept = data$null[[1]], linetype = "dashed", colour = "#777777") +
                ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = .15, orientation = "y", linewidth = .8) +
                ggplot2::geom_point(size = 3) +
                ggplot2::scale_colour_manual(values = c("AIPW" = "#0072B2", "IPTW (Hajek)" = "#E69F00")) +
                ggplot2::labs(x = unique(data$estimand), y = NULL, colour = NULL) +
                ggtheme +
                ggplot2::theme(legend.position = "none")
        },

        .syntaxText = function(fitted) {
            quoteText <- function(x) paste0('"', gsub('"', '\\\"', x), '"')
            vectorText <- function(x) if (length(x) == 0L) "NULL" else paste0("c(", paste(vapply(x, quoteText, character(1)), collapse = ", "), ")")
            prepared <- fitted$prepared
            paste(c(
                "# Reproduce the jamovi v1.1 point-outcome analysis",
                "result <- trialemulationj::aipw(",
                "  data = data,",
                paste0("  treatment = ", quoteText(prepared$originalTreatment), ","),
                paste0("  outcome = ", quoteText(prepared$originalOutcome), ","),
                paste0("  treatmentPositive = ", quoteText(self$options$treatmentPositive), ","),
                paste0("  outcomePositive = ", quoteText(self$options$outcomePositive), ","),
                paste0("  outcomeType = ", quoteText(self$options$outcomeType), ","),
                paste0("  outcomeCovs = ", vectorText(prepared$originalOutcomeCovs), ","),
                paste0("  propensityCovs = ", vectorText(prepared$originalPropensityCovs), ","),
                if (!is.null(prepared$originalCluster)) paste0("  cluster = ", quoteText(prepared$originalCluster), ",") else NULL,
                paste0("  absenceVars = ", vectorText(self$options$absenceVars), ","),
                paste0("  missingMethod = ", quoteText(self$options$missingMethod), ","),
                paste0("  imputations = ", self$options$imputations, ", randomSeed = ", self$options$randomSeed, ","),
                paste0("  outcomeSpecification = ", quoteText(self$options$outcomeSpecification), ","),
                paste0("  outcomeMethod = ", quoteText(self$options$outcomeMethod), ","),
                paste0("  propensityMethod = ", quoteText(self$options$propensityMethod), ","),
                paste0("  minPS = ", self$options$minPS, ", maxPS = ", self$options$maxPS, ","),
                paste0("  bootstrapSamples = ", self$options$bootstrapSamples, ","),
                paste0("  sensitivityGrid = ", if (isTRUE(self$options$sensitivityGrid)) "TRUE" else "FALSE", ","),
                paste0("  confidenceLevel = ", self$options$confidenceLevel),
                ")",
                paste0("# Specification fingerprint: ", fitted$fingerprint)
            ), collapse = "\n")
        }
    )
)
