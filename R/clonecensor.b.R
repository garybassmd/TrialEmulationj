#' @importFrom jmvcore .
cloneCensorClass <- R6::R6Class(
    "cloneCensorClass",
    inherit = cloneCensorBase,
    private = list(
        .run = function() {
            required <- c(self$options$id, self$options$event, self$options$timeToEvent,
                          self$options$exposure, self$options$timeToExposure)
            if (any(vapply(required, is.null, logical(1))))
                return()

            warningMessages <- character()
            fitted <- tryCatch(
                withCallingHandlers(
                    private$.fitAnalysis(),
                    warning = function(w) {
                        warningMessages <<- c(warningMessages, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                ),
                error = function(e) e
            )
            if (inherits(fitted, "error")) {
                jmvcore::reject(paste0(
                    "Clone-censor-weight analysis could not be completed: ",
                    conditionMessage(fitted),
                    " Check the baseline data layout, grace period, and censor-weight model."
                ))
            }
            private$.populateResults(fitted, warningMessages)
        },

        .asBinary = function(x, label) {
            if (is.logical(x))
                return(as.integer(x))
            if (is.factor(x))
                x <- as.character(x)
            values <- suppressWarnings(as.numeric(x))
            if (any(!is.finite(values)) || any(!values %in% c(0, 1)))
                stop(label, " must be coded 0 and 1.")
            as.integer(values)
        },

        .prepareData = function() {
            roles <- list(
                id = self$options$id,
                event = self$options$event,
                timeToEvent = self$options$timeToEvent,
                exposure = self$options$exposure,
                timeToExposure = self$options$timeToExposure
            )
            roleNames <- unlist(roles, use.names = FALSE)
            if (anyDuplicated(roleNames))
                stop("Each data-role box must use a different variable.")
            originalCovs <- setdiff(unique(self$options$weightCovs), roleNames)
            if (length(originalCovs) == 0)
                stop("Select at least one numeric baseline covariate for the artificial-censoring weights.")

            source <- as.data.frame(self$data[, unique(c(roleNames, originalCovs)), drop = FALSE])
            if (nrow(source) < 20)
                stop("At least 20 people are required.")
            if (anyDuplicated(source[[roles$id]]))
                stop("Input must contain one row per person; the person ID contains duplicates.")

            requiredComplete <- source[, setdiff(names(source), roles$timeToExposure), drop = FALSE]
            missing <- vapply(requiredComplete, function(x) sum(is.na(x)), integer(1))
            if (any(missing > 0))
                stop("Variables used by this analysis contain missing values: ", paste(names(missing)[missing > 0], collapse = ", "), ".")

            event <- private$.asBinary(source[[roles$event]], "Event")
            exposure <- private$.asBinary(source[[roles$exposure]], "Exposure")
            timeEvent <- suppressWarnings(as.numeric(source[[roles$timeToEvent]]))
            timeExposure <- suppressWarnings(as.numeric(source[[roles$timeToExposure]]))
            if (any(!is.finite(timeEvent)) || any(timeEvent <= 0))
                stop("Time to event or censoring must be finite and greater than zero.")
            if (any(exposure == 1 & !is.finite(timeExposure)))
                stop("Exposed people must have a finite time to exposure.")
            if (any(is.finite(timeExposure) & timeExposure < 0))
                stop("Time to exposure cannot be negative.")
            if (self$options$lowerPercentile >= self$options$upperPercentile)
                stop("The lower winsorization percentile must be below the upper percentile.")
            if (self$options$horizon > max(timeEvent))
                stop("The effect-estimation horizon cannot exceed the maximum observed follow-up time (", signif(max(timeEvent), 6), ").")

            data <- data.frame(
                .id = as.character(source[[roles$id]]),
                .event = event,
                .time_event = timeEvent,
                .exposure = exposure,
                .time_exposure = timeExposure,
                stringsAsFactors = FALSE
            )
            covNames <- sprintf(".w%03d", seq_along(originalCovs))
            names(covNames) <- originalCovs
            for (name in originalCovs) {
                if (!is.numeric(source[[name]]) && !is.logical(source[[name]]))
                    stop("Censor-weight covariates must be numeric. Create indicator variables for categorical covariates: ", name, ".")
                values <- as.numeric(source[[name]])
                if (any(!is.finite(values)))
                    stop("Censor-weight covariate contains non-finite values: ", name, ".")
                data[[covNames[[name]]]] <- values
            }

            list(data = data, covariates = unname(covNames), originalCovariates = originalCovs, roles = roles)
        },

        .fitOne = function(data, covariates, includeCurve = TRUE) {
            clones <- suppressMessages(survivalCCW::create_clones(
                df = data,
                id = ".id",
                event = ".event",
                time_to_event = ".time_event",
                exposure = ".exposure",
                time_to_exposure = ".time_exposure",
                ced_window = self$options$gracePeriod
            ))
            long <- survivalCCW::cast_clones_to_long(clones)
            long <- survivalCCW::generate_ccw(long, predvars = covariates)
            if (self$options$winsorize) {
                long <- survivalCCW::winsorize_ccw_weights(
                    long,
                    quantiles = c(self$options$lowerPercentile, self$options$upperPercentile) / 100,
                    per_clone = self$options$winsorizeByClone
                )
            }

            cox <- survival::coxph(
                survival::Surv(t_start, t_stop, outcome) ~ clone,
                data = long,
                weights = weight_cox,
                ties = "efron"
            )
            logHR <- unname(stats::coef(cox)[["clone"]])
            if (!is.finite(logHR))
                stop("The weighted Cox model did not produce a finite strategy coefficient.")

            fit0 <- survival::survfit(
                survival::Surv(t_start, t_stop, outcome) ~ 1,
                data = long[long$clone == 0, , drop = FALSE],
                weights = weight_cox
            )
            fit1 <- survival::survfit(
                survival::Surv(t_start, t_stop, outcome) ~ 1,
                data = long[long$clone == 1, , drop = FALSE],
                weights = weight_cox
            )
            horizon <- self$options$horizon
            surv0 <- unname(summary(fit0, times = horizon, extend = TRUE)$surv[[1]])
            surv1 <- unname(summary(fit1, times = horizon, extend = TRUE)$surv[[1]])
            rmst0 <- unname(summary(fit0, rmean = horizon, extend = TRUE)$table[["rmean"]])
            rmst1 <- unname(summary(fit1, rmean = horizon, extend = TRUE)$table[["rmean"]])

            point <- c(
                hr = exp(logHR),
                survival0 = surv0,
                survival1 = surv1,
                survivalDifference = surv1 - surv0,
                rmst0 = rmst0,
                rmst1 = rmst1,
                rmstDifference = rmst1 - rmst0
            )

            curve <- NULL
            if (includeCurve) {
                grid <- unique(c(0, seq(0, horizon, length.out = 121)))
                curve <- rbind(
                    data.frame(time = grid, survival = summary(fit0, times = grid, extend = TRUE)$surv, strategy = "No exposure during grace period"),
                    data.frame(time = grid, survival = summary(fit1, times = grid, extend = TRUE)$surv, strategy = "Exposure during grace period")
                )
            }

            list(point = point, logHR = logHR, long = long, curve = curve)
        },

        .bootstrap = function(data, covariates) {
            set.seed(self$options$seed)
            n <- nrow(data)
            B <- self$options$bootstrapSamples
            draws <- matrix(NA_real_, nrow = B, ncol = 7)
            colnames(draws) <- c("hr", "survival0", "survival1", "survivalDifference", "rmst0", "rmst1", "rmstDifference")
            for (b in seq_len(B)) {
                index <- sample.int(n, n, replace = TRUE)
                sampleData <- data[index, , drop = FALSE]
                sampleData$.id <- paste0("bootstrap_", b, "_", seq_len(n))
                estimate <- tryCatch(
                    suppressWarnings(private$.fitOne(sampleData, covariates, includeCurve = FALSE)$point),
                    error = function(e) NULL
                )
                if (!is.null(estimate))
                    draws[b, ] <- estimate[colnames(draws)]
            }
            draws
        },

        .fitAnalysis = function() {
            prepared <- private$.prepareData()
            fitted <- private$.fitOne(prepared$data, prepared$covariates, includeCurve = TRUE)
            draws <- private$.bootstrap(prepared$data, prepared$covariates)
            valid <- stats::complete.cases(draws)
            minimum <- max(10, ceiling(.5 * self$options$bootstrapSamples))
            if (sum(valid) < minimum)
                stop("Only ", sum(valid), " of ", self$options$bootstrapSamples, " bootstrap samples converged; simplify the weight model or increase the sample size.")
            draws <- draws[valid, , drop = FALSE]

            se <- apply(draws, 2, stats::sd)
            lower <- apply(draws, 2, stats::quantile, probs = .025)
            upper <- apply(draws, 2, stats::quantile, probs = .975)
            p <- rep(NaN, ncol(draws)); names(p) <- colnames(draws)
            p[["hr"]] <- 2 * stats::pnorm(-abs(fitted$logHR / stats::sd(log(draws[, "hr"]))))
            for (name in c("survivalDifference", "rmstDifference"))
                p[[name]] <- 2 * stats::pnorm(-abs(fitted$point[[name]] / se[[name]]))

            list(
                prepared = prepared,
                point = fitted$point,
                se = se,
                lower = lower,
                upper = upper,
                p = p,
                long = fitted$long,
                curve = fitted$curve,
                validBootstrap = nrow(draws)
            )
        },

        .populateResults = function(fitted, warningMessages) {
            long <- fitted$long
            data <- fitted$prepared$data
            rows <- list(
                c("People", format(nrow(data), big.mark = ",")),
                c("Observed outcome events", format(sum(data$.event), big.mark = ",")),
                c("Observed exposure during grace period", format(sum(data$.exposure), big.mark = ",")),
                c("Grace-period end", self$options$gracePeriod),
                c("Effect-estimation horizon", self$options$horizon),
                c("Clones", format(length(unique(interaction(long$.id, long$clone))), big.mark = ",")),
                c("Expanded clone-period rows", format(nrow(long), big.mark = ",")),
                c("Artificial censoring events", format(sum(long$censor), big.mark = ",")),
                c("Weight handling", if (self$options$winsorize) paste0("Winsorized at ", self$options$lowerPercentile, "th and ", self$options$upperPercentile, "th percentiles") else "No winsorization"),
                c("Successful bootstrap samples", paste0(fitted$validBootstrap, " / ", self$options$bootstrapSamples))
            )
            for (i in seq_along(rows))
                self$results$summary$setRow(rowNo = i, values = list(measure = rows[[i]][1], value = rows[[i]][2]))

            labels <- c(
                hr = "Weighted hazard ratio: exposure strategy / no-exposure strategy",
                survival0 = paste0("Survival at ", self$options$horizon, ": no-exposure strategy"),
                survival1 = paste0("Survival at ", self$options$horizon, ": exposure strategy"),
                survivalDifference = paste0("Survival difference at ", self$options$horizon, " (exposure - no exposure)"),
                rmst0 = paste0("RMST to ", self$options$horizon, ": no-exposure strategy"),
                rmst1 = paste0("RMST to ", self$options$horizon, ": exposure strategy"),
                rmstDifference = paste0("RMST difference to ", self$options$horizon, " (exposure - no exposure)")
            )
            for (name in names(labels)) {
                self$results$effects$addRow(
                    rowKey = name,
                    values = list(
                        estimand = labels[[name]], estimate = fitted$point[[name]],
                        se = fitted$se[[name]], lower = fitted$lower[[name]],
                        upper = fitted$upper[[name]], p = fitted$p[[name]]
                    )
                )
            }

            for (clone in 0:1) {
                weights <- long$weight_cox[long$clone == clone]
                q <- stats::quantile(weights, c(0, .01, .5, .99, 1), na.rm = TRUE)
                self$results$weightSummary$addRow(
                    rowKey = paste0("clone", clone),
                    values = list(
                        strategy = if (clone == 0) "No exposure during grace period" else "Exposure during grace period",
                        min = q[[1]], p1 = q[[2]], median = q[[3]], mean = mean(weights), p99 = q[[4]], max = q[[5]]
                    )
                )
            }
            self$results$survivalPlot$setState(fitted$curve)
            self$results$weightPlot$setState(data.frame(
                weight = long$weight_cox,
                strategy = factor(long$clone, levels = c(0, 1), labels = c("No exposure", "Exposure"))
            ))

            notes <- c(
                "Each person is cloned into two sustained strategies. Clones are artificially censored when their observed exposure history becomes incompatible with the assigned strategy.",
                "Percentile bootstrap intervals resample people before cloning; this accounts for dependence between the two clones from the same person.",
                "The hazard ratio is non-collapsible and relies on proportional hazards. Prefer the survival and RMST differences when an absolute marginal contrast is the scientific target.",
                "Causal interpretation requires a well-defined grace-period strategy, consistency, positivity, no interference, and no unmeasured predictors of artificial censoring and the outcome after conditioning on the selected baseline covariates.",
                "The survivalCCW weighting engine currently models baseline covariates with Cox models and is described by its maintainers as work in progress; validate results against protocol-specific code before clinical or regulatory use.",
                if (length(warningMessages) > 0) paste0("Model warning: ", unique(warningMessages)) else NULL
            )
            self$results$notes$setContent(paste0("- ", notes, collapse = "\n"))
            self$results$syntax$setContent(private$.syntaxText(fitted$prepared))
        },

        .survivalPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            ggplot2::ggplot(image$state, ggplot2::aes(x = time, y = survival, colour = strategy)) +
                ggplot2::geom_step(linewidth = .9) +
                ggplot2::scale_colour_manual(values = c("No exposure during grace period" = "#4477AA", "Exposure during grace period" = "#CC6677")) +
                ggplot2::coord_cartesian(ylim = c(0, 1)) +
                ggplot2::labs(x = "Follow-up time", y = "Weighted survival probability", colour = "Strategy") +
                ggtheme
        },

        .weightPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            ggplot2::ggplot(image$state, ggplot2::aes(x = weight, fill = strategy)) +
                ggplot2::geom_histogram(bins = 45, position = "identity", alpha = .55) +
                ggplot2::scale_fill_manual(values = c("No exposure" = "#4477AA", "Exposure" = "#CC6677")) +
                ggplot2::labs(x = "Clone-censor weight", y = "Clone-period rows", fill = "Strategy") +
                ggtheme
        },

        .syntaxText = function(prepared) {
            q <- function(x) paste0("\"", gsub("\"", "\\\\\"", x), "\"")
            covs <- paste(vapply(prepared$originalCovariates, q, character(1)), collapse = ", ")
            paste(c(
                "library(survivalCCW)",
                "library(survival)",
                "",
                "ccw <- data |>",
                "  create_clones(",
                paste0("    id = ", q(prepared$roles$id), ","),
                paste0("    event = ", q(prepared$roles$event), ","),
                paste0("    time_to_event = ", q(prepared$roles$timeToEvent), ","),
                paste0("    exposure = ", q(prepared$roles$exposure), ","),
                paste0("    time_to_exposure = ", q(prepared$roles$timeToExposure), ","),
                paste0("    ced_window = ", self$options$gracePeriod),
                "  ) |>",
                "  cast_clones_to_long() |>",
                paste0("  generate_ccw(c(", covs, "))"),
                if (self$options$winsorize) paste0("ccw <- winsorize_ccw_weights(ccw, c(", self$options$lowerPercentile / 100, ", ", self$options$upperPercentile / 100, "), per_clone = ", toupper(as.character(self$options$winsorizeByClone)), ")") else NULL,
                "fit <- coxph(Surv(t_start, t_stop, outcome) ~ clone, data = ccw, weights = weight_cox)",
                "# Resample people before recreating clones and weights for bootstrap inference."
            ), collapse = "\n")
        }
    )
)
