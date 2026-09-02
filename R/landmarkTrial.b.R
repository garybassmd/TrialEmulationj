#' @importFrom jmvcore .
landmarkTrialClass <- R6::R6Class(
    "landmarkTrialClass",
    inherit = landmarkTrialBase,
    private = list(
        .run = function() {
            if (is.null(self$options$treatmentA) || is.null(self$options$treatmentB))
                return()
            fitted <- tryCatch(private$.fitLandmarks(), error = function(e) e)
            if (inherits(fitted, "error"))
                jmvcore::reject(paste0("The landmark target-trial analysis could not be completed: ", conditionMessage(fitted)))
            private$.populate(fitted)
        },

        .optionVariables = function() {
            values <- c(
                self$options$eligibility, self$options$treatmentA, self$options$treatmentB,
                self$options$outcome, self$options$failureEvent1, self$options$failureEvent2,
                self$options$competingEvent, self$options$cluster, self$options$indexTime,
                self$options$treatmentATime, self$options$treatmentBTime,
                self$options$failureTime1, self$options$failureTime2,
                self$options$dischargeTime, self$options$competingTime,
                self$options$covariates, self$options$absenceVars
            )
            unique(values[!is.na(values) & nzchar(values)])
        },

        .binaryOrZero = function(source, name, label, positive = "", mask = NULL) {
            if (is.null(name))
                return(list(values = rep.int(0L, nrow(source)), zero = "No", one = "Yes", supplied = FALSE))
            result <- .te_as_binary(source[[name]], label, positive, allowMissing = TRUE, mask = mask)
            result$supplied <- TRUE
            result
        },

        .elapsedVariable = function(source, name) {
            if (is.null(name))
                return(rep(NA_real_, nrow(source)))
            index <- if (identical(self$options$timeMode, "dateTime")) {
                if (is.null(self$options$indexTime))
                    stop("Select an index/start date-time when date-time mode is used.")
                source[[self$options$indexTime]]
            } else NULL
            .te_elapsed(source[[name]], index, self$options$timeMode, self$options$timeUnit)
        },

        .prepareSource = function() {
            if (identical(self$options$outcomeSource, "recorded") && is.null(self$options$outcome))
                stop("Select a recorded failure outcome or choose composite outcome construction.")
            if (identical(self$options$outcomeSource, "composite") && is.null(self$options$failureEvent1) && is.null(self$options$failureEvent2))
                stop("Composite outcome construction requires at least one failure component.")
            if (self$options$minPS >= self$options$maxPS)
                stop("The minimum propensity score must be below the maximum propensity score.")
            if (isTRUE(self$options$runSecondWindow) && self$options$window2 <= self$options$window1)
                stop("The second decision window must be later than the primary window.")
            if (identical(self$options$missingMethod, "hotDeckMI") && self$options$bootstrapSamples > 0L)
                stop("Bootstrap confidence intervals and multiple hot-deck imputation cannot currently be combined in the landmark workflow.")

            variables <- private$.optionVariables()
            source <- as.data.frame(self$data[, variables, drop = FALSE])
            if (nrow(source) < 20L)
                stop("At least 20 rows are required.")
            originalN <- nrow(source)

            eligibility <- if (is.null(self$options$eligibility)) {
                list(values = rep.int(1L, originalN), zero = "Ineligible", one = "Eligible")
            } else .te_as_binary(source[[self$options$eligibility]], "Eligibility", self$options$eligibilityPositive, allowMissing = TRUE)
            baselineEligible <- eligibility$values == 1L & !is.na(eligibility$values)
            recoded <- .te_recode_absence(source, self$options$absenceVars, mask = baselineEligible)
            source <- recoded$data
            treatmentA <- .te_as_binary(source[[self$options$treatmentA]], "Treatment component A", self$options$treatmentAPositive, allowMissing = TRUE, mask = baselineEligible)
            treatmentB <- .te_as_binary(source[[self$options$treatmentB]], "Treatment component B", self$options$treatmentBPositive, allowMissing = TRUE, mask = baselineEligible)
            event1 <- private$.binaryOrZero(source, self$options$failureEvent1, "Failure component 1", self$options$outcomePositive, baselineEligible)
            event2 <- private$.binaryOrZero(source, self$options$failureEvent2, "Failure component 2", self$options$outcomePositive, baselineEligible)
            competing <- private$.binaryOrZero(source, self$options$competingEvent, "Competing event", self$options$outcomePositive, baselineEligible)

            if (identical(self$options$outcomeSource, "recorded")) {
                outcome <- .te_as_binary(source[[self$options$outcome]], "Recorded failure outcome", self$options$outcomePositive, allowMissing = TRUE, mask = baselineEligible)
                Y <- outcome$values
                outcomeMapping <- paste0(outcome$zero, " = 0; ", outcome$one, " = 1")
                if (is.null(self$options$failureEvent1) && !is.null(self$options$failureTime1)) {
                    event1$values <- Y
                    event1$supplied <- TRUE
                }
            } else {
                Y <- ifelse(event1$values == 1L | event2$values == 1L, 1L,
                            ifelse(!is.na(event1$values) & !is.na(event2$values), 0L, NA_integer_))
                outcomeMapping <- "Composite of selected failure components"
            }
            if (!is.null(self$options$competingEvent)) {
                if (identical(self$options$competingHandling, "noFailure"))
                    Y[competing$values == 1L] <- 0L
                if (identical(self$options$competingHandling, "composite"))
                    Y[competing$values == 1L] <- 1L
            }

            covariates <- setdiff(unique(self$options$covariates), c(
                self$options$treatmentA, self$options$treatmentB, self$options$outcome,
                self$options$failureEvent1, self$options$failureEvent2, self$options$competingEvent,
                self$options$cluster
            ))
            internal <- if (length(covariates) > 0L) sprintf("x%03d", seq_along(covariates)) else character()
            names(internal) <- covariates
            covData <- data.frame(row.names = seq_len(originalN))
            for (name in covariates)
                covData[[internal[[name]]]] <- source[[name]]
            labels <- setNames(covariates, unname(internal))

            cluster <- if (is.null(self$options$cluster)) NULL else as.character(.te_trim(source[[self$options$cluster]]))
            clusterMissing <- if (is.null(cluster)) rep(FALSE, originalN) else is.na(cluster) & baselineEligible
            if (any(clusterMissing)) {
                if (identical(self$options$missingCluster, "stop"))
                    stop(sum(clusterMissing), " baseline-eligible rows have missing cluster IDs. Choose 'Assign Unknown cluster' only when this matches the analysis plan.")
                cluster[clusterMissing] <- "Unknown"
            }

            list(
                source = source,
                originalN = originalN,
                recodeAudit = recoded$audit,
                eligibility = eligibility$values,
                treatmentA = treatmentA$values,
                treatmentB = treatmentB$values,
                event1 = event1$values,
                event2 = event2$values,
                competing = competing$values,
                Y = Y,
                treatmentAMapping = paste0(treatmentA$zero, " = absent; ", treatmentA$one, " = present"),
                treatmentBMapping = paste0(treatmentB$zero, " = unexposed; ", treatmentB$one, " = exposed"),
                outcomeMapping = outcomeMapping,
                covData = covData,
                covariates = covariates,
                internalCovs = unname(internal),
                labels = labels,
                cluster = cluster,
                clusterMissing = sum(clusterMissing),
                times = list(
                    treatmentA = private$.elapsedVariable(source, self$options$treatmentATime),
                    treatmentB = private$.elapsedVariable(source, self$options$treatmentBTime),
                    failure1 = private$.elapsedVariable(source, self$options$failureTime1),
                    failure2 = private$.elapsedVariable(source, self$options$failureTime2),
                    discharge = private$.elapsedVariable(source, self$options$dischargeTime),
                    competing = private$.elapsedVariable(source, self$options$competingTime)
                )
            )
        },

        .classifyTreatment = function(recorded, time, supplied, window) {
            classified <- recorded
            unknown <- is.na(recorded)
            if (!supplied)
                return(list(value = classified, unknown = unknown))
            treated <- recorded == 1L & !is.na(recorded)
            valid <- is.finite(time) & time >= 0
            classified[treated & valid & time <= window] <- 1L
            classified[treated & valid & time > window] <- 0L
            ambiguous <- treated & !valid
            if (identical(self$options$unknownTreatment, "assumeRecorded"))
                classified[ambiguous] <- 1L
            else {
                classified[ambiguous] <- NA_integer_
                unknown <- unknown | ambiguous
            }
            list(value = classified, unknown = unknown)
        },

        .eventTiming = function(flag, time, supplied, window) {
            if (!supplied)
                return(list(early = rep(FALSE, length(flag)), unknown = rep(FALSE, length(flag))))
            event <- flag == 1L & !is.na(flag)
            valid <- is.finite(time) & time >= 0
            list(
                early = event & valid & time <= window,
                unknown = event & !valid
            )
        },

        .auditTimeline = function(prepared) {
            rows <- list()
            add <- function(variable, issue, count, action) {
                if (count > 0L)
                    rows[[length(rows) + 1L]] <<- data.frame(variable = variable, issue = issue, n = count, action = action, stringsAsFactors = FALSE)
            }
            base <- prepared$eligibility == 1L & !is.na(prepared$eligibility)
            if (prepared$clusterMissing > 0L)
                add(self$options$cluster, "Missing/blank cluster ID", prepared$clusterMissing, "Assigned to one shared Unknown cluster")
            if (nrow(prepared$recodeAudit) > 0L) {
                for (i in seq_len(nrow(prepared$recodeAudit))) {
                    item <- prepared$recodeAudit[i, ]
                    add(item$variable, "Missing/blank value explicitly recoded as absence", item$recoded, item$rule)
                }
            }
            treatmentSpecs <- list(
                list(self$options$treatmentATime, prepared$treatmentA, prepared$times$treatmentA),
                list(self$options$treatmentBTime, prepared$treatmentB, prepared$times$treatmentB)
            )
            for (spec in treatmentSpecs) {
                if (is.null(spec[[1]])) next
                name <- spec[[1]]
                recorded <- spec[[2]]
                time <- spec[[3]]
                add(name, "Recorded treatment has missing/unparseable time", sum(base & recorded == 1L & !is.finite(time), na.rm = TRUE), if (self$options$unknownTreatment == "exclude") "Exclude at each landmark" else "Assume within window")
                add(name, "Negative time", sum(base & is.finite(time) & time < 0, na.rm = TRUE), if (self$options$unknownTreatment == "exclude") "Exclude if treatment recorded" else "Assume recorded classification")
                add(name, "Time present while treatment recorded absent", sum(base & recorded == 0L & is.finite(time) & time >= 0, na.rm = TRUE), "Retain recorded absence; review source data")
            }
            eventSpecs <- list(
                list(self$options$failureTime1, prepared$event1, prepared$times$failure1),
                list(self$options$failureTime2, prepared$event2, prepared$times$failure2),
                list(self$options$competingTime, prepared$competing, prepared$times$competing)
            )
            for (spec in eventSpecs) {
                if (is.null(spec[[1]])) next
                name <- spec[[1]]
                flag <- spec[[2]]
                time <- spec[[3]]
                add(name, "Recorded event has missing/unparseable time", sum(base & flag == 1L & !is.finite(time), na.rm = TRUE), if (self$options$unknownEvent == "exclude") "Exclude at each landmark" else "Retain recorded outcome")
                add(name, "Negative event time", sum(base & flag == 1L & is.finite(time) & time < 0, na.rm = TRUE), if (self$options$unknownEvent == "exclude") "Exclude at each landmark" else "Retain recorded outcome")
            }
            if (!is.null(self$options$dischargeTime)) {
                time <- prepared$times$discharge
                dischargeAction <- if (identical(self$options$unknownDischarge, "exclude")) "Exclude at each landmark" else "Retain; discharge before the landmark cannot be verified"
                add(self$options$dischargeTime, "Missing/unparseable discharge time", sum(base & !is.finite(time), na.rm = TRUE), dischargeAction)
                add(self$options$dischargeTime, "Negative discharge time", sum(base & is.finite(time) & time < 0, na.rm = TRUE), dischargeAction)
            }
            if (length(rows) == 0L)
                data.frame(variable = "All selected time variables", issue = "No flagged timing conflicts", n = 0L, action = "None", stringsAsFactors = FALSE)
            else do.call(rbind, rows)
        },

        .constructWindow = function(prepared, window) {
            base <- prepared$eligibility == 1L & !is.na(prepared$eligibility)
            classA <- private$.classifyTreatment(prepared$treatmentA, prepared$times$treatmentA, !is.null(self$options$treatmentATime), window)
            classB <- private$.classifyTreatment(prepared$treatmentB, prepared$times$treatmentB, !is.null(self$options$treatmentBTime), window)
            treatmentUnknown <- classA$unknown | classB$unknown

            event1 <- private$.eventTiming(prepared$event1, prepared$times$failure1, !is.null(self$options$failureTime1), window)
            event2 <- private$.eventTiming(prepared$event2, prepared$times$failure2, !is.null(self$options$failureTime2), window)
            competing <- private$.eventTiming(prepared$competing, prepared$times$competing, !is.null(self$options$competingTime), window)
            dischargeSupplied <- !is.null(self$options$dischargeTime)
            dischargeEarly <- if (!dischargeSupplied) rep(FALSE, prepared$originalN) else is.finite(prepared$times$discharge) & prepared$times$discharge >= 0 & prepared$times$discharge <= window
            dischargeUnknown <- if (!dischargeSupplied) rep(FALSE, prepared$originalN) else !is.finite(prepared$times$discharge) | prepared$times$discharge < 0
            early <- event1$early | event2$early | competing$early | dischargeEarly
            eventUnknown <- event1$unknown | event2$unknown | competing$unknown
            outcomeUnknown <- is.na(prepared$Y)
            competingExcluded <- identical(self$options$competingHandling, "exclude") & prepared$competing == 1L & !is.na(prepared$competing)

            sequential <- base
            flow <- list(data.frame(step = "Rows received after jamovi filters", n = prepared$originalN, excluded = 0L))
            previous <- sum(sequential)
            flow[[2]] <- data.frame(step = "Meet baseline eligibility", n = previous, excluded = prepared$originalN - previous)
            sequential <- sequential & !treatmentUnknown
            current <- sum(sequential)
            flow[[3]] <- data.frame(step = "Treatment strategy classifiable", n = current, excluded = previous - current)
            previous <- current
            sequential <- sequential & !early
            if (identical(self$options$unknownDischarge, "exclude"))
                sequential <- sequential & !dischargeUnknown
            current <- sum(sequential)
            flow[[4]] <- data.frame(step = "Remain observed and event-free at landmark", n = current, excluded = previous - current)
            previous <- current
            if (identical(self$options$unknownEvent, "exclude"))
                sequential <- sequential & !eventUnknown
            sequential <- sequential & !outcomeUnknown
            current <- sum(sequential)
            flow[[5]] <- data.frame(step = "Outcome and event timing usable", n = current, excluded = previous - current)
            previous <- current
            sequential <- sequential & !competingExcluded
            current <- sum(sequential)
            flow[[6]] <- data.frame(step = "Final landmark cohort", n = current, excluded = previous - current)
            flow <- do.call(rbind, flow)
            flow$window <- private$.windowLabel(window)
            flow <- flow[, c("window", "step", "n", "excluded")]

            list(
                window = window,
                label = private$.windowLabel(window),
                keep = sequential,
                treatmentA = classA$value,
                treatmentB = classB$value,
                Y = prepared$Y,
                flow = flow
            )
        },

        .windowLabel = function(window) {
            paste0(format(window, trim = TRUE), " ", self$options$timeUnit)
        },

        .strata = function() {
            if (isTRUE(self$options$runBothStrata)) c(1L, 0L) else if (self$options$primaryStratum == "present") 1L else 0L
        },

        .componentLabel = function(value, fallback) {
            value <- trimws(as.character(value))
            if (length(value) == 0L || is.na(value) || !nzchar(value)) fallback else value
        },

        .labelA = function() private$.componentLabel(self$options$treatmentALabel, "Treatment A"),

        .labelB = function() private$.componentLabel(self$options$treatmentBLabel, "Treatment B"),

        .stratumLabel = function(stratum) {
            paste(private$.labelA(), if (stratum == 1L) "present" else "absent")
        },

        .fitSpecification = function(prepared, cohort, stratum) {
            keep <- cohort$keep & cohort$treatmentA == stratum & !is.na(cohort$treatmentA)
            data <- data.frame(A = cohort$treatmentB[keep], Y = cohort$Y[keep])
            if (length(prepared$internalCovs) > 0L)
                data <- cbind(data, prepared$covData[keep, prepared$internalCovs, drop = FALSE])
            if (!is.null(prepared$cluster))
                data$.cluster <- prepared$cluster[keep]
            covMissing <- if (length(prepared$internalCovs) == 0L) rep(FALSE, nrow(data)) else !stats::complete.cases(data[, prepared$internalCovs, drop = FALSE])
            excludedMissing <- sum(covMissing)
            if (any(covMissing) && identical(self$options$missingMethod, "fail"))
                stop(excludedMissing, " rows in the ", cohort$label, ", ", private$.stratumLabel(stratum), " stratum have missing baseline covariates.")
            if (any(covMissing) && identical(self$options$missingMethod, "completeCases"))
                data <- data[!covMissing, , drop = FALSE]
            if (nrow(data) < 20L || length(unique(data$A)) < 2L || min(table(data$A)) < 5L)
                return(list(error = "Insufficient observations or treatment overlap", data = data, excludedMissing = excludedMissing))

            fitArgs <- list(
                outcomeType = "binary",
                outcomeVars = prepared$internalCovs,
                propensityVars = prepared$internalCovs,
                minPS = self$options$minPS,
                maxPS = self$options$maxPS,
                outcomeMethod = self$options$outcomeMethod,
                propensityMethod = self$options$propensityMethod,
                outcomeSpecification = "pooled",
                confidenceLevel = self$options$confidenceLevel
            )
            if (any(covMissing) && identical(self$options$missingMethod, "hotDeckMI")) {
                imputations <- .te_hotdeck_imputations(data, prepared$internalCovs, self$options$imputations, self$options$randomSeed + round(cohort$window) + stratum, strata = data$A)
                fits <- lapply(imputations, function(item) do.call(.te_fit_causal_once, c(list(data = item, cluster = if (is.null(prepared$cluster)) NULL else item$.cluster), fitArgs)))
                fitted <- fits[[1]]
                fitted$effects <- .te_pool_mi_effects(fits, self$options$confidenceLevel)
                analysisData <- imputations[[1]]
            } else {
                fitted <- do.call(.te_fit_causal_once, c(list(data = data, cluster = if (is.null(prepared$cluster)) NULL else data$.cluster), fitArgs))
                draws <- .te_bootstrap_causal(data, fitArgs, self$options$bootstrapSamples, self$options$randomSeed + round(cohort$window) + stratum, if (is.null(prepared$cluster)) NULL else ".cluster")
                fitted$effects <- .te_apply_bootstrap_ci(fitted$effects, draws, self$options$confidenceLevel)
                analysisData <- data
            }
            fitted$balance$covariate <- private$.restoreLabels(fitted$balance$covariate, prepared$labels)
            list(fitted = fitted, data = analysisData, excludedMissing = excludedMissing)
        },

        .restoreLabels = function(labels, mapping) {
            for (internal in names(mapping))
                labels <- sub(paste0("^", internal), mapping[[internal]], labels)
            labels
        },

        .selectedEstimand = function() {
            switch(self$options$effectMeasure, riskDifference = "Risk difference", riskRatio = "Risk ratio", oddsRatio = "Odds ratio")
        },

        .fitLandmarks = function() {
            prepared <- private$.prepareSource()
            windows <- c(self$options$window1, if (isTRUE(self$options$runSecondWindow)) self$options$window2 else numeric())
            cohorts <- lapply(windows, function(window) private$.constructWindow(prepared, window))
            timeline <- private$.auditTimeline(prepared)
            flow <- do.call(rbind, lapply(cohorts, `[[`, "flow"))

            strategyRows <- list()
            effectRows <- list()
            specFits <- list()
            warnings <- character()
            for (cohort in cohorts) {
                for (a in 0:1) {
                    for (b in 0:1) {
                        keep <- cohort$keep & cohort$treatmentA == a & cohort$treatmentB == b
                        strategyRows[[length(strategyRows) + 1L]] <- data.frame(
                            window = cohort$label,
                            treatmentA = paste(private$.labelA(), if (a == 1) "present" else "absent"),
                            treatmentB = paste(private$.labelB(), if (b == 1) "exposed" else "unexposed"),
                            n = sum(keep, na.rm = TRUE),
                            failures = sum(cohort$Y[keep] == 1L, na.rm = TRUE),
                            stringsAsFactors = FALSE
                        )
                    }
                }
                for (stratum in private$.strata()) {
                    spec <- private$.fitSpecification(prepared, cohort, stratum)
                    key <- paste(cohort$window, stratum, sep = "|")
                    specFits[[key]] <- spec
                    if (!is.null(spec$error)) {
                        warnings <- c(warnings, paste0(cohort$label, ", ", private$.stratumLabel(stratum), ": ", spec$error, "."))
                        next
                    }
                    desired <- private$.selectedEstimand()
                    selected <- spec$fitted$effects[spec$fitted$effects$estimand == desired, , drop = FALSE]
                    if (!"bootstrapSuccess" %in% names(selected))
                        selected$bootstrapSuccess <- NA_integer_
                    for (i in seq_len(nrow(selected))) {
                        effectRows[[length(effectRows) + 1L]] <- data.frame(
                            window = cohort$label,
                            stratum = private$.stratumLabel(stratum),
                            method = selected$method[[i]],
                            n = nrow(spec$data),
                            treated = sum(spec$data$A == 1L),
                            failures = sum(spec$data$Y == 1L),
                            estimate = selected$estimate[[i]],
                            lower = selected$lower[[i]],
                            upper = selected$upper[[i]],
                            outcomeModel = spec$fitted$outcomeMethod,
                            bootstrapSuccess = selected$bootstrapSuccess[[i]],
                            stringsAsFactors = FALSE
                        )
                    }
                    warnings <- c(warnings, spec$fitted$warnings)
                }
            }
            strategies <- do.call(rbind, strategyRows)
            effects <- if (length(effectRows) == 0L) data.frame() else do.call(rbind, effectRows)
            primaryWindowValue <- if (self$options$primaryWindow == "second" && length(windows) > 1L) windows[[2]] else windows[[1]]
            primaryStratumValue <- if (self$options$primaryStratum == "present") 1L else 0L
            primaryKey <- paste(primaryWindowValue, primaryStratumValue, sep = "|")
            primaryFit <- specFits[[primaryKey]]
            balance <- if (is.null(primaryFit) || !is.null(primaryFit$error)) data.frame() else primaryFit$fitted$balance
            fingerprint <- .te_specification_fingerprint(list(
                analysis = "landmarkTrial", version = "1.1.1",
                variables = private$.optionVariables(), windows = windows,
                positiveLevels = list(
                    eligibility = self$options$eligibilityPositive,
                    treatmentA = self$options$treatmentAPositive,
                    treatmentB = self$options$treatmentBPositive,
                    outcome = self$options$outcomePositive
                ),
                displayLabels = list(treatmentA = private$.labelA(), treatmentB = private$.labelB()),
                outcomeSource = self$options$outcomeSource,
                competingHandling = self$options$competingHandling,
                timeMode = self$options$timeMode, timeUnit = self$options$timeUnit,
                unknownTreatment = self$options$unknownTreatment,
                unknownEvent = self$options$unknownEvent,
                unknownDischarge = self$options$unknownDischarge,
                missingCluster = self$options$missingCluster,
                strata = private$.strata(), covariates = self$options$covariates,
                absenceVariables = self$options$absenceVars,
                missing = list(method = self$options$missingMethod, imputations = self$options$imputations),
                models = list(outcome = self$options$outcomeMethod, propensity = self$options$propensityMethod),
                propensityLimits = c(self$options$minPS, self$options$maxPS),
                effect = self$options$effectMeasure, bootstrap = self$options$bootstrapSamples,
                randomSeed = self$options$randomSeed, confidenceLevel = self$options$confidenceLevel,
                cluster = self$options$cluster
            ))
            list(
                prepared = prepared, windows = windows, cohorts = cohorts,
                flow = flow, timeline = timeline, strategies = strategies,
                effects = effects, balance = balance, warnings = unique(warnings[nzchar(warnings)]),
                fingerprint = fingerprint
            )
        },

        .populate = function(fitted) {
            prepared <- fitted$prepared
            finalCounts <- vapply(fitted$cohorts, function(x) sum(x$keep), integer(1))
            rows <- list(
                c("Rows received after jamovi filters", format(prepared$originalN, big.mark = ",")),
                c("Decision windows", paste(vapply(fitted$windows, private$.windowLabel, character(1)), collapse = ", ")),
                c(paste(private$.labelA(), "coding"), prepared$treatmentAMapping),
                c(paste(private$.labelB(), "coding"), prepared$treatmentBMapping),
                c("Outcome", prepared$outcomeMapping),
                c("Final cohort size(s)", paste(finalCounts, collapse = ", ")),
                c("Analyzed strata", paste(vapply(private$.strata(), private$.stratumLabel, character(1)), collapse = "; ")),
                c("Primary effect", private$.selectedEstimand()),
                c("Inference", if (is.null(prepared$cluster)) {
                    if (self$options$bootstrapSamples > 0L) "Independent-row influence-function SE and row bootstrap" else "Independent-row influence-function SE"
                } else {
                    paste0(
                        length(unique(prepared$cluster[prepared$eligibility == 1L])),
                        " clusters; cluster-robust influence-function SE",
                        if (self$options$bootstrapSamples > 0L) " and whole-cluster bootstrap" else ""
                    )
                }),
                c("Missing cluster ID policy", if (is.null(self$options$cluster)) "Not applicable" else if (self$options$missingCluster == "unknown") "Assign one shared Unknown cluster" else "Stop and report"),
                c("Missing/invalid discharge policy", if (is.null(self$options$dischargeTime)) "No discharge variable supplied" else if (self$options$unknownDischarge == "exclude") "Exclude when landmark observation cannot be verified" else "Retain and flag"),
                c("Missing covariates", switch(self$options$missingMethod, fail = "Stop and report", completeCases = "Complete-case", hotDeckMI = paste0(self$options$imputations, " hot-deck imputations"))),
                c("Audit flags across variables (counts may overlap)", sum(fitted$timeline$n)),
                c("Specification fingerprint", substr(fitted$fingerprint, 1, 16))
            )
            for (i in seq_along(rows))
                self$results$summary$setRow(rowNo = i, values = list(measure = rows[[i]][1], value = rows[[i]][2]))
            for (i in seq_len(nrow(fitted$flow)))
                self$results$flow$addRow(rowKey = paste0("flow", i), values = as.list(fitted$flow[i, ]))
            for (i in seq_len(nrow(fitted$timeline)))
                self$results$timeline$addRow(rowKey = paste0("timeline", i), values = as.list(fitted$timeline[i, ]))
            for (i in seq_len(nrow(fitted$strategies)))
                self$results$strategies$addRow(rowKey = paste0("strategy", i), values = as.list(fitted$strategies[i, ]))
            if (nrow(fitted$effects) > 0L) {
                for (i in seq_len(nrow(fitted$effects)))
                    self$results$effects$addRow(rowKey = paste0("effect", i), values = as.list(fitted$effects[i, ]))
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
            self$results$flowPlot$setState(fitted$flow)
            self$results$balancePlot$setState(fitted$balance)
            effectPlot <- if (nrow(fitted$effects) == 0L) fitted$effects else fitted$effects[fitted$effects$method == "AIPW", , drop = FALSE]
            if (nrow(effectPlot) > 0L)
                effectPlot$null <- if (self$options$effectMeasure == "riskDifference") 0 else 1
            self$results$effectPlot$setState(effectPlot)
            self$results$notes$setContent(private$.notes(fitted))
            self$results$syntax$setContent(private$.syntax(fitted))
        },

        .notes = function(fitted) {
            notes <- c(
                "Eligibility, treatment classification, and start of outcome follow-up are aligned at each landmark. Events or discharge at/before the landmark are excluded before treatment effects are estimated.",
                "The treatment-initiation comparison is an observational analogue of an assignment effect, not a literal intention-to-treat effect unless randomized assignment is observed.",
                if (is.null(self$options$treatmentATime)) paste0("No timing variable was supplied for ", private$.labelA(), "; recorded status is assumed to describe strategy membership by each landmark.") else NULL,
                if (is.null(self$options$treatmentBTime)) paste0("No timing variable was supplied for ", private$.labelB(), "; recorded status is assumed to describe strategy membership by each landmark.") else NULL,
                if (is.null(self$options$failureTime1) && is.null(self$options$failureTime2)) "No failure-time variable was supplied. Outcome-free status at each landmark cannot be verified; this is acceptable only when failure cannot occur before the landmark by design." else NULL,
                paste0("Treatment with missing/invalid timing is handled by: ", self$options$unknownTreatment, ". Failure with missing/invalid timing is handled by: ", self$options$unknownEvent, ". Missing/invalid discharge timing is handled by: ", self$options$unknownDischarge, "."),
                "Eligibility is applied before treatment and outcome coding are validated, so values observed only among ineligible rows do not invalidate the eligible cohort.",
                "Competing-event handling is a prespecified estimand choice. Excluding competing events conditions on post-baseline information and is provided only as a sensitivity analysis, not as a default causal estimand.",
                if (is.null(fitted$prepared$cluster)) {
                    "Rows are treated as independent."
                } else if (self$options$bootstrapSamples > 0L) {
                    "Cluster-robust influence functions and whole-cluster bootstrap preserve within-cluster dependence."
                } else {
                    "Cluster-robust influence functions preserve within-cluster dependence."
                },
                if (self$options$bootstrapSamples > 0L) paste0("Percentile confidence limits use up to ", self$options$bootstrapSamples, " refitted bootstrap samples per landmark/stratum specification.") else "Confidence limits use clustered or independent influence-function standard errors.",
                "Firth mean bias-reduced logistic regression is used automatically when outcome or propensity separation/instability is detected.",
                "Inspect the timeline audit, cohort flow, overlap, effective sample size, and balance before interpreting effects. A balance threshold does not establish exchangeability.",
                "Only baseline causes measured before treatment classification belong in the adjustment set. Do not use outcome-informed transformations, mediators, colliders, or post-landmark information.",
                if (length(fitted$warnings) > 0L) paste0("Analysis warning: ", paste(fitted$warnings, collapse = "; ")) else NULL
            )
            paste0("- ", notes, collapse = "\n")
        },

        .flowPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || nrow(image$state) == 0L)
                return(FALSE)
            data <- image$state
            data$step <- factor(data$step, levels = unique(data$step))
            ggplot2::ggplot(data, ggplot2::aes(x = step, y = n, group = window, colour = window)) +
                ggplot2::geom_line(linewidth = .8) +
                ggplot2::geom_point(size = 2.8) +
                ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.7, size = 3) +
                ggplot2::scale_colour_manual(values = c("#E69F00", "#0072B2")) +
                ggplot2::labs(x = NULL, y = "Patients remaining", colour = "Decision window") +
                ggtheme +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
        },

        .balancePlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || nrow(image$state) == 0L)
                return(FALSE)
            state <- image$state
            order <- state$covariate[order(abs(state$unweighted))]
            data <- rbind(
                data.frame(covariate = state$covariate, smd = state$unweighted, sample = "Unadjusted"),
                data.frame(covariate = state$covariate, smd = state$weighted, sample = "Adjusted")
            )
            data$covariate <- factor(data$covariate, levels = unique(order))
            ggplot2::ggplot(data, ggplot2::aes(x = smd, y = covariate, colour = sample)) +
                ggplot2::geom_vline(xintercept = 0, colour = "#333333", linewidth = .45) +
                ggplot2::geom_vline(xintercept = c(-.10, .10), linetype = "dashed", colour = "#777777", linewidth = .45) +
                ggplot2::geom_point(size = 2.7) +
                ggplot2::scale_colour_manual(values = c("Unadjusted" = "#F8766D", "Adjusted" = "#00AEB3")) +
                ggplot2::labs(x = "Standardized mean difference", y = NULL, colour = "Sample") +
                ggtheme
        },

        .effectPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state) || nrow(image$state) == 0L)
                return(FALSE)
            data <- image$state
            data$specification <- paste(data$window, data$stratum, sep = " - ")
            data$specification <- factor(data$specification, levels = rev(unique(data$specification)))
            colours <- stats::setNames(c("#E69F00", "#0072B2"), c(private$.stratumLabel(1L), private$.stratumLabel(0L)))
            ggplot2::ggplot(data, ggplot2::aes(x = estimate, y = specification, colour = stratum)) +
                ggplot2::geom_vline(xintercept = data$null[[1]], linetype = "dashed", colour = "#777777") +
                ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = .18, orientation = "y", linewidth = .85) +
                ggplot2::geom_point(size = 3.2) +
                ggplot2::scale_colour_manual(values = colours) +
                ggplot2::labs(x = private$.selectedEstimand(), y = NULL, colour = "Subgroup") +
                ggtheme
        },

        .syntax = function(fitted) {
            value <- function(x) if (is.null(x)) "NULL" else paste0('"', gsub('"', '\\\"', x), '"')
            vector <- function(x) if (length(x) == 0L) "NULL" else paste0("c(", paste(vapply(x, value, character(1)), collapse = ", "), ")")
            paste(c(
                "# Reproduce the TrialEmulationj 1.1.1 landmark specification",
                "result <- trialemulationj::landmarkTrial(",
                "  data = data,",
                paste0("  eligibility = ", value(self$options$eligibility), ","),
                paste0("  treatmentA = ", value(self$options$treatmentA), ","),
                paste0("  treatmentB = ", value(self$options$treatmentB), ","),
                paste0("  outcome = ", value(self$options$outcome), ","),
                paste0("  failureEvent1 = ", value(self$options$failureEvent1), ","),
                paste0("  failureEvent2 = ", value(self$options$failureEvent2), ","),
                paste0("  competingEvent = ", value(self$options$competingEvent), ","),
                paste0("  cluster = ", value(self$options$cluster), ", missingCluster = ", value(self$options$missingCluster), ","),
                paste0("  indexTime = ", value(self$options$indexTime), ","),
                paste0("  treatmentATime = ", value(self$options$treatmentATime), ","),
                paste0("  treatmentBTime = ", value(self$options$treatmentBTime), ","),
                paste0("  failureTime1 = ", value(self$options$failureTime1), ","),
                paste0("  failureTime2 = ", value(self$options$failureTime2), ","),
                paste0("  dischargeTime = ", value(self$options$dischargeTime), ","),
                paste0("  competingTime = ", value(self$options$competingTime), ","),
                paste0("  covariates = ", vector(self$options$covariates), ","),
                paste0("  absenceVars = ", vector(self$options$absenceVars), ","),
                paste0("  eligibilityPositive = ", value(self$options$eligibilityPositive), ","),
                paste0("  treatmentAPositive = ", value(self$options$treatmentAPositive), ","),
                paste0("  treatmentBPositive = ", value(self$options$treatmentBPositive), ","),
                paste0("  outcomePositive = ", value(self$options$outcomePositive), ","),
                paste0("  treatmentALabel = ", value(self$options$treatmentALabel), ", treatmentBLabel = ", value(self$options$treatmentBLabel), ","),
                paste0("  outcomeSource = ", value(self$options$outcomeSource), ","),
                paste0("  competingHandling = ", value(self$options$competingHandling), ","),
                paste0("  timeMode = ", value(self$options$timeMode), ", timeUnit = ", value(self$options$timeUnit), ","),
                paste0("  window1 = ", self$options$window1, ", runSecondWindow = ", if (isTRUE(self$options$runSecondWindow)) "TRUE" else "FALSE", ", window2 = ", self$options$window2, ","),
                paste0("  unknownTreatment = ", value(self$options$unknownTreatment), ", unknownEvent = ", value(self$options$unknownEvent), ", unknownDischarge = ", value(self$options$unknownDischarge), ","),
                paste0("  runBothStrata = ", if (isTRUE(self$options$runBothStrata)) "TRUE" else "FALSE", ", primaryStratum = ", value(self$options$primaryStratum), ","),
                paste0("  primaryWindow = ", value(self$options$primaryWindow), ","),
                paste0("  missingMethod = ", value(self$options$missingMethod), ", imputations = ", self$options$imputations, ","),
                paste0("  outcomeMethod = ", value(self$options$outcomeMethod), ", propensityMethod = ", value(self$options$propensityMethod), ","),
                paste0("  minPS = ", self$options$minPS, ", maxPS = ", self$options$maxPS, ","),
                paste0("  effectMeasure = ", value(self$options$effectMeasure), ","),
                paste0("  bootstrapSamples = ", self$options$bootstrapSamples, ", randomSeed = ", self$options$randomSeed, ","),
                paste0("  confidenceLevel = ", self$options$confidenceLevel),
                ")",
                paste0("# Specification SHA-256: ", fitted$fingerprint)
            ), collapse = "\n")
        }
    )
)
