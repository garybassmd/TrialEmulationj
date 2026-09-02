#' @importFrom jmvcore .
targetTrialClass <- R6::R6Class(
    "targetTrialClass",
    inherit = targetTrialBase,
    private = list(
        .run = function() {
            required <- c(
                self$options$id,
                self$options$period,
                self$options$treatment,
                self$options$outcome,
                self$options$eligible
            )

            if (any(vapply(required, is.null, logical(1))))
                return()

            prepared <- private$.prepareData()
            private$.fillInputSummary(prepared$data, prepared$roles)

            warningMessages <- character()
            fitted <- tryCatch(
                withCallingHandlers(
                    private$.fitTrial(prepared$data, prepared$roles),
                    warning = function(w) {
                        warningMessages <<- c(warningMessages, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                ),
                error = function(e) e
            )

            if (inherits(fitted, "error")) {
                jmvcore::reject(paste0(
                    "Trial emulation could not be completed: ",
                    conditionMessage(fitted),
                    " Check the variable coding and model settings."
                ))
            }

            private$.populateResults(fitted, warningMessages)
        },

        .prepareData = function() {
            roles <- list(
                id = self$options$id,
                period = self$options$period,
                treatment = self$options$treatment,
                outcome = self$options$outcome,
                eligible = self$options$eligible,
                censor = if (self$options$useCensor) self$options$censor else NULL
            )

            roleNames <- unlist(roles, use.names = FALSE)
            roleNames <- roleNames[!vapply(roleNames, is.null, logical(1))]
            if (anyDuplicated(roleNames))
                jmvcore::reject("Each data-role box must use a different variable.")

            if (self$options$useCensor && is.null(roles$censor))
                jmvcore::reject("Select a censoring event variable or turn off informative-censoring adjustment.")

            covariates <- unique(c(
                self$options$outcomeCovs,
                if (self$options$estimand == "pp") self$options$switchNumCovs else NULL,
                if (self$options$estimand == "pp") self$options$switchDenCovs else NULL,
                if (self$options$useCensor) self$options$censorNumCovs else NULL,
                if (self$options$useCensor) self$options$censorDenCovs else NULL
            ))
            covariates <- setdiff(covariates, roleNames)
            used <- unique(c(roleNames, covariates))

            data <- as.data.frame(self$data[, used, drop = FALSE])
            if (nrow(data) == 0)
                jmvcore::reject("The data set contains no rows.")

            missingByVariable <- vapply(data, function(x) sum(is.na(x)), integer(1))
            if (any(missingByVariable > 0)) {
                affected <- paste(names(missingByVariable)[missingByVariable > 0], collapse = ", ")
                jmvcore::reject(paste0(
                    "Variables used by this analysis contain missing values: ", affected,
                    ". Resolve these values before emulating trials; silently deleting person-period rows can change eligibility and follow-up."
                ))
            }

            data[[roles$treatment]] <- private$.asBinary(data[[roles$treatment]], "Treatment")
            data[[roles$outcome]] <- private$.asBinary(data[[roles$outcome]], "Outcome")
            data[[roles$eligible]] <- private$.asBinary(data[[roles$eligible]], "Eligibility")
            if (self$options$useCensor)
                data[[roles$censor]] <- private$.asBinary(data[[roles$censor]], "Censoring")

            period <- suppressWarnings(as.numeric(as.character(data[[roles$period]])))
            if (any(!is.finite(period)) || any(period < 0) || any(period != floor(period)))
                jmvcore::reject("Period must contain non-negative integers with no missing values.")
            data[[roles$period]] <- as.integer(period)

            if (is.factor(data[[roles$id]]))
                data[[roles$id]] <- as.character(data[[roles$id]])

            duplicateKey <- duplicated(data[c(roles$id, roles$period)])
            if (any(duplicateKey))
                jmvcore::reject("The data contain duplicate person-period rows. Each ID and period combination must be unique.")

            ordering <- order(data[[roles$id]], data[[roles$period]])
            data <- data[ordering, , drop = FALSE]
            rownames(data) <- NULL

            list(data = data, roles = roles)
        },

        .asBinary = function(x, label) {
            if (is.logical(x))
                return(as.integer(x))

            if (is.factor(x))
                x <- as.character(x)

            values <- suppressWarnings(as.numeric(x))
            if (any(!is.finite(values)) || any(!values %in% c(0, 1)))
                jmvcore::reject(paste0(label, " must be coded 0 and 1."))
            as.integer(values)
        },

        .rhsFormula = function(vars) {
            vars <- unique(vars)
            if (length(vars) == 0)
                return(~ 1)
            stats::reformulate(vars)
        },

        .timeFormula = function(choice, variable) {
            if (choice == "none")
                return(~ 1)
            if (choice == "linear")
                return(stats::reformulate(variable))
            if (choice == "quadratic")
                return(stats::as.formula(paste0("~ ", variable, " + I(", variable, "^2)")))
            stats::as.formula(paste0("~ splines::ns(", variable, ", df = 3)"))
        },

        .fitTrial = function(data, roles) {
            estimand <- if (self$options$estimand == "itt") "ITT" else "PP"
            modelPath <- file.path(tempdir(), paste0("trialemulationj_", as.integer(Sys.time())))
            dir.create(modelPath, recursive = TRUE, showWarnings = FALSE)

            trial <- TrialEmulation::trial_sequence(estimand) |>
                TrialEmulation::set_data(
                    data = data,
                    id = roles$id,
                    period = roles$period,
                    treatment = roles$treatment,
                    outcome = roles$outcome,
                    eligible = roles$eligible
                )

            if (estimand == "PP") {
                trial <- TrialEmulation::set_switch_weight_model(
                    trial,
                    numerator = private$.rhsFormula(self$options$switchNumCovs),
                    denominator = private$.rhsFormula(self$options$switchDenCovs),
                    model_fitter = TrialEmulation::stats_glm_logit(
                        save_path = file.path(modelPath, "adherence")
                    )
                )
            }

            if (self$options$useCensor) {
                pool <- self$options$poolModels
                if (pool == "auto")
                    pool <- if (estimand == "ITT") "numerator" else "none"

                trial <- TrialEmulation::set_censor_weight_model(
                    trial,
                    censor_event = roles$censor,
                    numerator = private$.rhsFormula(self$options$censorNumCovs),
                    denominator = private$.rhsFormula(self$options$censorDenCovs),
                    pool_models = pool,
                    model_fitter = TrialEmulation::stats_glm_logit(
                        save_path = file.path(modelPath, "censoring")
                    )
                )
            }

            if (estimand == "PP" || self$options$useCensor)
                trial <- TrialEmulation::calculate_weights(trial, quiet = TRUE)

            trial <- TrialEmulation::set_outcome_model(
                trial,
                adjustment_terms = private$.rhsFormula(self$options$outcomeCovs),
                followup_time_terms = private$.timeFormula(self$options$followupModel, "followup_time"),
                trial_period_terms = private$.timeFormula(self$options$trialPeriodModel, "trial_period"),
                model_fitter = TrialEmulation::stats_glm_logit(save_path = NA)
            ) |>
                TrialEmulation::set_expansion_options(
                    output = TrialEmulation::save_to_datatable(),
                    chunk_size = self$options$chunkSize
                ) |>
                TrialEmulation::expand_trials() |>
                TrialEmulation::load_expanded_data(
                    seed = self$options$seed,
                    p_control = self$options$pControl
                )

            expanded <- TrialEmulation::outcome_data(trial)
            if (nrow(expanded) == 0)
                stop("No person-trial observations remained after trial expansion.")

            ipw <- if ("weight" %in% names(expanded)) expanded$weight else rep(1, nrow(expanded))
            sampling <- if ("sample_weight" %in% names(expanded)) expanded$sample_weight else rep(1, nrow(expanded))
            rawFinal <- ipw * sampling
            final <- rawFinal

            modifier <- NULL
            truncationPoint <- NA_real_
            if (self$options$truncateWeights) {
                probability <- self$options$truncatePercentile / 100
                truncationPoint <- as.numeric(stats::quantile(rawFinal, probability, na.rm = TRUE))
                modifier <- function(w) pmin(w, truncationPoint)
                final <- modifier(rawFinal)
            }

            trial <- TrialEmulation::fit_msm(
                trial,
                weight_cols = c("weight", "sample_weight"),
                modify_weights = modifier
            )

            baseline <- TrialEmulation::outcome_data(trial)
            baseline <- baseline[baseline$followup_time == 0, , drop = FALSE]
            populationLabel <- "All baseline person-trials"
            if (self$options$targetPopulation == "earliest") {
                targetPeriod <- min(baseline$trial_period)
                baseline <- baseline[baseline$trial_period == targetPeriod, , drop = FALSE]
                populationLabel <- paste0("Earliest emulated trial (period ", targetPeriod, ")")
            } else if (self$options$targetPopulation == "latest") {
                targetPeriod <- max(baseline$trial_period)
                baseline <- baseline[baseline$trial_period == targetPeriod, , drop = FALSE]
                populationLabel <- paste0("Latest emulated trial (period ", targetPeriod, ")")
            } else if (self$options$targetPopulation == "specific") {
                targetPeriod <- self$options$targetTrialPeriod
                baseline <- baseline[baseline$trial_period == targetPeriod, , drop = FALSE]
                populationLabel <- paste0("Specific emulated trial (period ", targetPeriod, ")")
            }

            if (nrow(baseline) == 0)
                stop("The selected target population contains no baseline person-trials.")

            observedMax <- max(TrialEmulation::outcome_data(trial)$followup_time)
            requestedMax <- self$options$maxFollowUp
            predictionMax <- min(requestedMax, observedMax)
            predictionTimes <- 0:predictionMax
            set.seed(self$options$seed)
            predictions <- TrialEmulation::predict(
                trial,
                newdata = baseline,
                predict_times = predictionTimes,
                conf_int = self$options$confInt,
                samples = self$options$ciSamples,
                ci_type = "sandwich",
                type = self$options$predictionType
            )

            list(
                trial = trial,
                expanded = expanded,
                predictions = predictions,
                weights = list(ipw = ipw, sampling = sampling, rawFinal = rawFinal, final = final),
                truncationPoint = truncationPoint,
                populationN = nrow(baseline),
                populationLabel = populationLabel,
                predictionMax = predictionMax,
                requestedMax = requestedMax,
                observedMax = observedMax,
                estimand = estimand,
                roles = roles,
                data = data
            )
        },

        .fillInputSummary = function(data, roles) {
            rows <- list(
                c("Person-period rows", format(nrow(data), big.mark = ",")),
                c("Unique persons", format(length(unique(data[[roles$id]])), big.mark = ",")),
                c("Observed periods", format(length(unique(data[[roles$period]])), big.mark = ",")),
                c("Eligible person-periods", format(sum(data[[roles$eligible]] == 1), big.mark = ",")),
                c("Eligible treatment initiations", format(sum(data[[roles$eligible]] == 1 & data[[roles$treatment]] == 1), big.mark = ",")),
                c("Outcome events", format(sum(data[[roles$outcome]] == 1), big.mark = ",")),
                c("Censoring events", if (is.null(roles$censor)) "Not modelled" else format(sum(data[[roles$censor]] == 1), big.mark = ","))
            )
            table <- self$results$inputSummary
            for (i in seq_along(rows))
                table$setRow(rowNo = i, values = list(measure = rows[[i]][1], value = rows[[i]][2]))
        },

        .populateResults = function(fitted, warningMessages) {
            outcomeData <- TrialEmulation::outcome_data(fitted$trial)
            summaryRows <- list(
                c("Estimand", if (fitted$estimand == "ITT") "Marginal intention-to-treat" else "Marginal per-protocol"),
                c("Target population", fitted$populationLabel),
                c("Target-population rows", format(fitted$populationN, big.mark = ",")),
                c("Emulated trial starts", format(length(unique(outcomeData$trial_period)), big.mark = ",")),
                c("Expanded analysis rows", format(nrow(outcomeData), big.mark = ",")),
                c("Expanded outcome events", format(sum(outcomeData$outcome == 1), big.mark = ",")),
                c("Observed maximum follow-up", as.character(fitted$observedMax)),
                c("Weight truncation", if (self$options$truncateWeights) paste0("Upper ", self$options$truncatePercentile, "% (cap ", signif(fitted$truncationPoint, 5), ")") else "None")
            )
            for (i in seq_along(summaryRows))
                self$results$analysisSummary$setRow(rowNo = i, values = list(measure = summaryRows[[i]][1], value = summaryRows[[i]][2]))

            tidy <- fitted$trial@outcome_model@fitted@summary$tidy
            for (i in seq_len(nrow(tidy))) {
                self$results$coefficients$addRow(
                    rowKey = paste0("coef", i),
                    values = list(
                        term = tidy$term[i],
                        estimate = tidy$estimate[i],
                        se = tidy$std.error[i],
                        statistic = tidy$statistic[i],
                        p = tidy$p.value[i],
                        lower = tidy$conf.low[i],
                        upper = tidy$conf.high[i]
                    )
                )
            }

            pred <- fitted$predictions
            valueName <- self$options$predictionType
            diffName <- paste0(valueName, "_diff")
            untreated <- pred$assigned_treatment_0
            treated <- pred$assigned_treatment_1
            difference <- pred$difference

            for (i in seq_len(nrow(difference))) {
                lower <- if (self$options$confInt) difference[["2.5%"]][i] else NaN
                upper <- if (self$options$confInt) difference[["97.5%"]][i] else NaN
                self$results$marginal$addRow(
                    rowKey = paste0("time", difference$followup_time[i]),
                    values = list(
                        time = difference$followup_time[i],
                        untreated = untreated[[valueName]][i],
                        treated = treated[[valueName]][i],
                        difference = difference[[diffName]][i],
                        lower = lower,
                        upper = upper
                    )
                )
            }

            curveData <- rbind(
                private$.predictionFrame(untreated, valueName, "No treatment"),
                private$.predictionFrame(treated, valueName, "Treatment")
            )
            diffData <- data.frame(
                time = difference$followup_time,
                estimate = difference[[diffName]],
                lower = if (self$options$confInt) difference[["2.5%"]] else NA_real_,
                upper = if (self$options$confInt) difference[["97.5%"]] else NA_real_
            )
            plotState <- list(
                curves = curveData,
                difference = diffData,
                type = valueName,
                confInt = self$options$confInt,
                estimand = fitted$estimand
            )
            self$results$riskPlot$setState(plotState)
            self$results$differencePlot$setState(plotState)

            weights <- list(
                "Inverse-probability" = fitted$weights$ipw,
                "Sampling" = fitted$weights$sampling,
                "Raw combined" = fitted$weights$rawFinal,
                "Final analysis" = fitted$weights$final
            )
            for (i in seq_along(weights)) {
                x <- weights[[i]]
                q <- stats::quantile(x, probs = c(0, .01, .05, .5, .95, .99, 1), na.rm = TRUE, names = FALSE)
                self$results$weightSummary$addRow(
                    rowKey = paste0("weight", i),
                    values = list(
                        weight = names(weights)[i],
                        min = q[1], p1 = q[2], p5 = q[3], median = q[4],
                        mean = mean(x, na.rm = TRUE), p95 = q[5], p99 = q[6], max = q[7]
                    )
                )
            }
            self$results$weightPlot$setState(fitted$weights$final)

            notes <- c(
                "Causal interpretation requires consistency, positivity, no interference, and no unmeasured confounding of treatment assignment at each emulated trial baseline.",
                if (fitted$estimand == "PP") "The per-protocol estimand additionally requires positivity and conditional exchangeability of adherence." else NULL,
                if (self$options$useCensor) "Informative-censoring adjustment additionally requires positivity and sequential exchangeability for censoring, plus correctly specified denominator weight models." else NULL,
                "The reported contrast is Treatment minus No treatment on the selected prediction scale.",
                "Inspect the weight diagnostics for extreme values; truncation changes the target estimator and should be justified prospectively.",
                if (fitted$requestedMax > fitted$predictionMax) paste0("Prediction follow-up was capped at ", fitted$predictionMax, ", the maximum observed follow-up in the expanded data.") else NULL,
                if (self$options$pControl < 1) paste0("Non-event person-periods were sampled with probability ", self$options$pControl, "; inverse sampling weights were included.") else NULL
            )
            if (length(warningMessages) > 0)
                notes <- c(notes, paste0("Model warning: ", unique(warningMessages)))
            self$results$notes$setContent(paste0("- ", notes, collapse = "\n"))

            self$results$syntax$setContent(private$.syntaxText(fitted))
        },

        .predictionFrame = function(frame, valueName, group) {
            data.frame(
                time = frame$followup_time,
                estimate = frame[[valueName]],
                lower = if (self$options$confInt) frame[["2.5%"]] else NA_real_,
                upper = if (self$options$confInt) frame[["97.5%"]] else NA_real_,
                strategy = group
            )
        },

        .riskPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            state <- image$state
            yLabel <- if (state$type == "cum_inc") "Cumulative incidence" else "Survival probability"
            plot <- ggplot2::ggplot(
                state$curves,
                ggplot2::aes(x = time, y = estimate, colour = strategy, fill = strategy)
            )
            if (state$confInt) {
                plot <- plot + ggplot2::geom_ribbon(
                    ggplot2::aes(ymin = lower, ymax = upper),
                    alpha = .15,
                    colour = NA
                )
            }
            plot +
                ggplot2::geom_line(linewidth = .9) +
                ggplot2::scale_colour_manual(values = c("No treatment" = "#4477AA", "Treatment" = "#CC6677")) +
                ggplot2::scale_fill_manual(values = c("No treatment" = "#4477AA", "Treatment" = "#CC6677")) +
                ggplot2::labs(x = "Follow-up period", y = yLabel, colour = "Strategy", fill = "Strategy") +
                ggtheme
        },

        .differencePlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            state <- image$state
            yLabel <- if (state$type == "cum_inc") "Cumulative incidence difference" else "Survival probability difference"
            plot <- ggplot2::ggplot(state$difference, ggplot2::aes(x = time, y = estimate))
            if (state$confInt) {
                plot <- plot + ggplot2::geom_ribbon(
                    ggplot2::aes(ymin = lower, ymax = upper),
                    fill = "#4477AA",
                    alpha = .18
                )
            }
            plot +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#666666") +
                ggplot2::geom_line(linewidth = .9, colour = "#4477AA") +
                ggplot2::labs(x = "Follow-up period", y = yLabel) +
                ggtheme
        },

        .weightPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            data <- data.frame(weight = image$state)
            ggplot2::ggplot(data, ggplot2::aes(x = weight)) +
                ggplot2::geom_histogram(bins = 40, fill = "#4477AA", colour = "white") +
                ggplot2::labs(x = "Final analysis weight", y = "Person-period rows") +
                ggtheme
        },

        .syntaxText = function(fitted) {
            quoteName <- function(x) paste0("\"", gsub("\"", "\\\\\"", x), "\"")
            formulaText <- function(vars) paste(deparse(private$.rhsFormula(vars)), collapse = "")
            lines <- c(
                "library(TrialEmulation)",
                "",
                paste0("trial <- trial_sequence(\"", fitted$estimand, "\") |>"),
                "  set_data(",
                "    data = data,",
                paste0("    id = ", quoteName(fitted$roles$id), ","),
                paste0("    period = ", quoteName(fitted$roles$period), ","),
                paste0("    treatment = ", quoteName(fitted$roles$treatment), ","),
                paste0("    outcome = ", quoteName(fitted$roles$outcome), ","),
                paste0("    eligible = ", quoteName(fitted$roles$eligible)),
                "  )"
            )

            if (fitted$estimand == "PP") {
                lines <- c(lines,
                    "",
                    "trial <- trial |>",
                    "  set_switch_weight_model(",
                    paste0("    numerator = ", formulaText(self$options$switchNumCovs), ","),
                    paste0("    denominator = ", formulaText(self$options$switchDenCovs), ","),
                    "    model_fitter = stats_glm_logit(tempdir())",
                    "  )"
                )
            }

            if (self$options$useCensor) {
                pool <- self$options$poolModels
                if (pool == "auto")
                    pool <- if (fitted$estimand == "ITT") "numerator" else "none"
                lines <- c(lines,
                    "",
                    "trial <- trial |>",
                    "  set_censor_weight_model(",
                    paste0("    censor_event = ", quoteName(fitted$roles$censor), ","),
                    paste0("    numerator = ", formulaText(self$options$censorNumCovs), ","),
                    paste0("    denominator = ", formulaText(self$options$censorDenCovs), ","),
                    paste0("    pool_models = \"", pool, "\","),
                    "    model_fitter = stats_glm_logit(tempdir())",
                    "  )"
                )
            }

            if (fitted$estimand == "PP" || self$options$useCensor)
                lines <- c(lines, "", "trial <- calculate_weights(trial, quiet = TRUE)")

            followup <- paste(deparse(private$.timeFormula(self$options$followupModel, "followup_time")), collapse = "")
            trialPeriod <- paste(deparse(private$.timeFormula(self$options$trialPeriodModel, "trial_period")), collapse = "")
            lines <- c(lines,
                "",
                "trial <- trial |>",
                "  set_outcome_model(",
                paste0("    adjustment_terms = ", formulaText(self$options$outcomeCovs), ","),
                paste0("    followup_time_terms = ", followup, ","),
                paste0("    trial_period_terms = ", trialPeriod),
                "  ) |>",
                paste0("  set_expansion_options(save_to_datatable(), chunk_size = ", self$options$chunkSize, ") |>"),
                "  expand_trials() |>",
                paste0("  load_expanded_data(seed = ", self$options$seed, ", p_control = ", self$options$pControl, ")")
            )

            if (self$options$truncateWeights) {
                lines <- c(lines,
                    "",
                    paste0("cap <- quantile(outcome_data(trial)$weight * outcome_data(trial)$sample_weight, ", self$options$truncatePercentile / 100, ")"),
                    "trial <- fit_msm(trial, modify_weights = function(w) pmin(w, cap))"
                )
            } else {
                lines <- c(lines, "", "trial <- fit_msm(trial)")
            }

            lines <- c(lines,
                "",
                "baseline <- outcome_data(trial)",
                "baseline <- baseline[baseline$followup_time == 0, ]"
            )
            if (self$options$targetPopulation == "earliest")
                lines <- c(lines, "baseline <- baseline[baseline$trial_period == min(baseline$trial_period), ]")
            if (self$options$targetPopulation == "latest")
                lines <- c(lines, "baseline <- baseline[baseline$trial_period == max(baseline$trial_period), ]")
            if (self$options$targetPopulation == "specific")
                lines <- c(lines, paste0("baseline <- baseline[baseline$trial_period == ", self$options$targetTrialPeriod, ", ]"))
            lines <- c(lines,
                paste0("set.seed(", self$options$seed, ")"),
                "marginal <- predict(",
                "  trial,",
                "  newdata = baseline,",
                paste0("  predict_times = 0:", fitted$predictionMax, ","),
                paste0("  conf_int = ", toupper(as.character(self$options$confInt)), ","),
                paste0("  samples = ", self$options$ciSamples, ","),
                paste0("  type = \"", self$options$predictionType, "\""),
                ")"
            )
            paste(lines, collapse = "\n")
        }
    )
)
