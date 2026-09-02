#' @importFrom jmvcore .
drSurvivalClass <- R6::R6Class(
    "drSurvivalClass",
    inherit = drSurvivalBase,
    private = list(
        .run = function() {
            required <- c(self$options$time, self$options$event, self$options$treatment)
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
                    "Doubly robust survival estimation could not be completed: ",
                    conditionMessage(fitted),
                    " Check follow-up, treatment overlap, and model specifications."
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

        .rhsText = function(vars) {
            if (length(vars) == 0) "1" else paste(vars, collapse = " + ")
        },

        .prepareData = function() {
            roles <- list(
                time = self$options$time,
                event = self$options$event,
                treatment = self$options$treatment,
                followup = self$options$followupTime
            )
            roleNames <- unlist(roles, use.names = FALSE)
            roleNames <- roleNames[!vapply(roleNames, is.null, logical(1))]
            if (anyDuplicated(roleNames))
                stop("Each data-role box must use a different variable.")

            originalCovs <- unique(c(self$options$outcomeCovs, self$options$propensityCovs, self$options$censorCovs))
            originalCovs <- setdiff(originalCovs, roleNames)
            used <- unique(c(roleNames, originalCovs))
            source <- as.data.frame(self$data[, used, drop = FALSE])
            if (nrow(source) < 30)
                stop("At least 30 complete observations are required.")
            missing <- vapply(source, function(x) sum(is.na(x)), integer(1))
            if (any(missing > 0))
                stop("Variables used by this analysis contain missing values: ", paste(names(missing)[missing > 0], collapse = ", "), ".")

            y <- suppressWarnings(as.numeric(source[[roles$time]]))
            if (any(!is.finite(y)) || any(y <= 0))
                stop("Follow-up time must be numeric, finite, and greater than zero.")
            d <- private$.asBinary(source[[roles$event]], "Event")
            trt <- private$.asBinary(source[[roles$treatment]], "Treatment")
            if (min(table(trt)) < 10)
                stop("Each treatment group must contain at least ten observations.")
            if (sum(d) < 10)
                stop("At least ten outcome events are required.")
            if (self$options$horizon > max(y))
                stop("The RMTL horizon cannot exceed the maximum observed follow-up time (", signif(max(y), 6), ").")
            if (self$options$minPS >= self$options$maxPS)
                stop("The minimum propensity score must be below the maximum propensity score.")

            data <- data.frame(y = y, d = d, trt = trt)
            if (!is.null(roles$followup)) {
                yf <- suppressWarnings(as.numeric(source[[roles$followup]]))
                if (any(!is.finite(yf)) || any(yf < y))
                    stop("Potential censoring time must be finite and at least as large as observed follow-up time.")
                data$yf <- yf
            }

            internal <- if (length(originalCovs) > 0) sprintf("x%03d", seq_along(originalCovs)) else character()
            names(internal) <- originalCovs
            for (name in originalCovs)
                data[[internal[[name]]]] <- source[[name]]

            list(
                data = data,
                roles = roles,
                outcomeVars = unname(internal[intersect(self$options$outcomeCovs, originalCovs)]),
                propensityVars = unname(internal[intersect(self$options$propensityCovs, originalCovs)]),
                censorVars = unname(internal[intersect(self$options$censorCovs, originalCovs)]),
                originalOutcomeCovs = self$options$outcomeCovs,
                originalPropensityCovs = self$options$propensityCovs,
                originalCensorCovs = self$options$censorCovs
            )
        },

        .fitAnalysis = function() {
            prepared <- private$.prepareData()
            cateFormula <- stats::as.formula(paste0("survival::Surv(y, d) ~ ", private$.rhsText(prepared$outcomeVars)))
            psFormula <- stats::as.formula(paste0("trt ~ ", private$.rhsText(prepared$propensityVars)))
            ipcwFormula <- if (length(prepared$censorVars) == 0) NULL else stats::as.formula(paste0("~ trt + ", private$.rhsText(prepared$censorVars)))
            censorMethod <- switch(
                self$options$censorMethod,
                breslow = "breslow",
                aft_exp = "aft (exponential)",
                aft_weibull = "aft (weibull)",
                aft_lognormal = "aft (lognormal)",
                aft_loglogistic = "aft (loglogistic)"
            )

            fit <- precmed::atefitsurv(
                data = prepared$data,
                cate.model = cateFormula,
                ps.model = psFormula,
                ps.method = self$options$psMethod,
                ipcw.model = ipcwFormula,
                ipcw.method = censorMethod,
                minPS = self$options$minPS,
                maxPS = self$options$maxPS,
                followup.time = if (is.null(prepared$roles$followup)) NULL else "yf",
                tau0 = self$options$horizon,
                surv.min = self$options$minCensorSurvival,
                n.boot = self$options$bootstrapSamples,
                seed = self$options$seed,
                verbose = 0
            )

            list(fit = fit, prepared = prepared, censorMethod = censorMethod)
        },

        .rowValues = function(frame, transform = FALSE, test = TRUE) {
            estimate <- frame$estimate[[1]]
            se <- frame$SE[[1]]
            lower <- frame$CI.lower[[1]]
            upper <- frame$CI.upper[[1]]
            if (transform) {
                se <- exp(estimate) * se
                estimate <- exp(estimate)
                lower <- exp(lower)
                upper <- exp(upper)
            }
            list(estimate = estimate, se = se, lower = lower, upper = upper, p = if (test) frame$pvalue[[1]] else NaN)
        },

        .populateResults = function(fitted, warningMessages) {
            data <- fitted$prepared$data
            fit <- fitted$fit
            rows <- list(
                c("People", format(nrow(data), big.mark = ",")),
                c("Untreated", format(sum(data$trt == 0), big.mark = ",")),
                c("Treated", format(sum(data$trt == 1), big.mark = ",")),
                c("Outcome events", format(sum(data$d), big.mark = ",")),
                c("RMTL horizon", self$options$horizon),
                c("Propensity model", if (length(fitted$prepared$originalPropensityCovs) == 0) "Intercept only" else paste(fitted$prepared$originalPropensityCovs, collapse = ", ")),
                c("Outcome model", if (length(fitted$prepared$originalOutcomeCovs) == 0) "Intercept only" else paste(fitted$prepared$originalOutcomeCovs, collapse = ", ")),
                c("Censoring model", paste0(fitted$censorMethod, if (length(fitted$prepared$originalCensorCovs) == 0) " (default covariates)" else paste0(": ", paste(fitted$prepared$originalCensorCovs, collapse = ", ")))),
                c("Bootstrap samples", self$options$bootstrapSamples)
            )
            for (i in seq_along(rows))
                self$results$summary$setRow(rowNo = i, values = list(measure = rows[[i]][1], value = rows[[i]][2]))

            effectFrames <- list(
                rmst0 = list(label = paste0("RMST to ", self$options$horizon, ": no treatment"), frame = fit$rmst0, transform = FALSE, test = FALSE),
                rmst1 = list(label = paste0("RMST to ", self$options$horizon, ": treatment"), frame = fit$rmst1, transform = FALSE, test = FALSE),
                rmtl = list(label = "RMTL ratio: treatment / no treatment", frame = fit$log.rmtl.ratio, transform = TRUE, test = TRUE),
                hr = list(label = "Adjusted hazard ratio: treatment / no treatment", frame = fit$log.hazard.ratio, transform = TRUE, test = TRUE)
            )
            for (name in names(effectFrames)) {
                spec <- effectFrames[[name]]
                values <- private$.rowValues(spec$frame, spec$transform, spec$test)
                self$results$effects$addRow(
                    rowKey = name,
                    values = c(list(estimand = spec$label), values)
                )
            }

            bootstrap <- as.data.frame(fit$trt.boot)
            plotData <- rbind(
                data.frame(estimate = bootstrap$log.rmtl.ratio, estimand = "Log RMTL ratio"),
                data.frame(estimate = bootstrap$log.hazard.ratio, estimand = "Log hazard ratio")
            )
            self$results$bootstrapPlot$setState(plotData)

            notes <- c(
                "The primary doubly robust estimand is the restricted-mean-time-lost (RMTL) ratio: treatment divided by no treatment. Values above 1 indicate more event-free time lost under treatment by the selected horizon.",
                "RMST is restricted mean survival time; larger values indicate longer event-free survival. The adjusted hazard ratio is supplementary and relies on proportional-hazards interpretation.",
                "Bootstrap uncertainty is produced by precmed::atefitsurv. The estimator combines an outcome regression with inverse-probability treatment and censoring adjustment.",
                "Double robustness protects against misspecification of one nuisance-model component, not violations of consistency, positivity, no interference, or conditional exchangeability.",
                "Inspect sensitivity to the RMTL horizon, propensity truncation, censoring model, and covariate sets. Use subject-level rows only.",
                if (!is.null(fit$warning) && nzchar(fit$warning)) paste0("Estimator warning: ", trimws(fit$warning)) else NULL,
                if (length(warningMessages) > 0) paste0("Model warning: ", unique(warningMessages)) else NULL
            )
            self$results$notes$setContent(paste0("- ", notes, collapse = "\n"))
            self$results$syntax$setContent(private$.syntaxText(fitted$prepared, fitted$censorMethod))
        },

        .bootstrapPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            ggplot2::ggplot(image$state, ggplot2::aes(x = estimate, fill = estimand)) +
                ggplot2::geom_histogram(bins = 35, colour = "white", alpha = .75) +
                ggplot2::facet_wrap(~ estimand, scales = "free_x", ncol = 1) +
                ggplot2::scale_fill_manual(values = c("Log RMTL ratio" = "#4477AA", "Log hazard ratio" = "#CC6677")) +
                ggplot2::labs(x = "Bootstrap estimate", y = "Samples") +
                ggplot2::guides(fill = "none") +
                ggtheme
        },

        .syntaxText = function(prepared, censorMethod) {
            bt <- function(x) paste0("`", gsub("`", "", x), "`")
            rhs <- function(vars) if (length(vars) == 0) "1" else paste(bt(vars), collapse = " + ")
            ipcw <- if (length(prepared$originalCensorCovs) == 0) "NULL" else paste0("~ ", bt(prepared$roles$treatment), " + ", rhs(prepared$originalCensorCovs))
            paste(c(
                "library(precmed)",
                "library(survival)",
                "",
                "fit <- atefitsurv(",
                "  data = data,",
                paste0("  cate.model = Surv(", bt(prepared$roles$time), ", ", bt(prepared$roles$event), ") ~ ", rhs(prepared$originalOutcomeCovs), ","),
                paste0("  ps.model = ", bt(prepared$roles$treatment), " ~ ", rhs(prepared$originalPropensityCovs), ","),
                paste0("  ps.method = \"", self$options$psMethod, "\","),
                paste0("  ipcw.model = ", ipcw, ","),
                paste0("  ipcw.method = \"", censorMethod, "\","),
                paste0("  minPS = ", self$options$minPS, ", maxPS = ", self$options$maxPS, ","),
                if (!is.null(prepared$roles$followup)) paste0("  followup.time = ", deparse(prepared$roles$followup), ",") else NULL,
                paste0("  tau0 = ", self$options$horizon, ", surv.min = ", self$options$minCensorSurvival, ","),
                paste0("  n.boot = ", self$options$bootstrapSamples, ", seed = ", self$options$seed),
                ")"
            ), collapse = "\n")
        }
    )
)
