.te_trim <- function(x) {
    if (is.factor(x))
        x <- as.character(x)
    if (is.character(x)) {
        x <- trimws(x)
        x[x == ""] <- NA_character_
    }
    x
}

.te_binary_label_roles <- function(levels) {
    if (length(levels) != 2L)
        return(NULL)
    normal <- tolower(gsub("[^[:alnum:]]+", " ", levels))
    normal <- trimws(normal)
    zeroWords <- c(
        "0", "no", "false", "control", "untreated", "unexposed", "none",
        "alive", "alive at discharge", "successful", "successful nom",
        "no event", "no treatment", "absence", "absent"
    )
    oneWords <- c("1", "yes", "true", "treated", "exposed", "event", "failure", "died", "death")
    zeroHit <- which(normal %in% zeroWords)
    oneHit <- which(normal %in% oneWords)
    if (length(zeroHit) == 1L)
        return(list(zero = levels[[zeroHit]], one = setdiff(levels, levels[[zeroHit]])[[1]]))
    if (length(oneHit) == 1L)
        return(list(one = levels[[oneHit]], zero = setdiff(levels, levels[[oneHit]])[[1]]))
    NULL
}

.te_as_binary <- function(x, label, positive = "", allowMissing = FALSE, mask = NULL) {
    x <- .te_trim(x)
    if (!is.null(mask)) {
        if (length(mask) != length(x))
            stop("Internal error: the eligibility mask for ", label, " has the wrong length.")
        mask <- as.logical(mask)
        mask[is.na(mask)] <- FALSE
        x[!mask] <- NA
    }
    missing <- is.na(x)
    if (any(missing) && !allowMissing)
        stop(label, " contains missing values.")

    observed <- x[!missing]
    if (is.logical(observed)) {
        values <- as.integer(x)
        return(list(values = values, zero = "FALSE", one = "TRUE"))
    }

    if (is.numeric(observed) || is.integer(observed)) {
        uniqueValues <- sort(unique(as.numeric(observed)))
        if (length(uniqueValues) != 2L || !all(uniqueValues %in% c(0, 1)))
            stop(label, " must have exactly two levels coded 0 and 1, or two labelled levels.")
        return(list(values = as.integer(x), zero = "0", one = "1"))
    }

    text <- as.character(x)
    levels <- unique(text[!is.na(text)])
    if (length(levels) != 2L)
        stop(label, " must have exactly two observed levels.")

    if (nzchar(trimws(positive))) {
        hit <- which(tolower(levels) == tolower(trimws(positive)))
        if (length(hit) != 1L)
            stop("The requested positive level for ", label, " was not found. Available levels: ", paste(levels, collapse = ", "), ".")
        one <- levels[[hit]]
        zero <- setdiff(levels, one)[[1]]
    } else {
        roles <- .te_binary_label_roles(levels)
        if (!is.null(roles)) {
            zero <- roles$zero
            one <- roles$one
        } else {
            zero <- sort(levels)[[1]]
            one <- sort(levels)[[2]]
        }
    }

    values <- rep.int(NA_integer_, length(text))
    values[text == zero & !is.na(text)] <- 0L
    values[text == one & !is.na(text)] <- 1L
    list(values = values, zero = zero, one = one)
}

.te_recode_absence <- function(data, variables, mask = NULL) {
    variables <- intersect(unique(variables), names(data))
    audit <- data.frame(variable = character(), recoded = integer(), rule = character(), stringsAsFactors = FALSE)
    if (length(variables) == 0L)
        return(list(data = data, audit = audit))
    if (is.null(mask)) {
        mask <- rep(TRUE, nrow(data))
    } else {
        if (length(mask) != nrow(data))
            stop("Internal error: the absence-recoding mask has the wrong length.")
        mask <- as.logical(mask)
        mask[is.na(mask)] <- FALSE
    }

    for (name in variables) {
        x <- data[[name]]
        blank <- is.na(x) & mask
        if (is.character(x) || is.factor(x))
            blank <- blank | (mask & trimws(as.character(x)) == "")
        count <- sum(blank)
        if (count > 0L) {
            if (is.logical(x)) {
                x[blank] <- FALSE
                rule <- "Missing/blank -> FALSE"
            } else if (is.numeric(x) || is.integer(x)) {
                x[blank] <- 0
                rule <- "Missing -> 0"
            } else {
                x <- as.character(x)
                observed <- unique(x[mask & !blank & !is.na(x)])
                roles <- .te_binary_label_roles(observed)
                replacement <- if (is.null(roles)) "No" else roles$zero
                x[blank] <- replacement
                x <- factor(x)
                rule <- paste0("Missing/blank -> ", replacement)
            }
            data[[name]] <- x
        } else {
            rule <- "No values changed"
        }
        audit <- rbind(audit, data.frame(variable = name, recoded = count, rule = rule, stringsAsFactors = FALSE))
    }
    list(data = data, audit = audit)
}

.te_hotdeck_imputations <- function(data, variables, m = 20L, seed = 24680L, strata = NULL) {
    variables <- intersect(unique(variables), names(data))
    if (length(variables) == 0L || !anyNA(data[, variables, drop = FALSE]))
        return(list(data))
    m <- max(2L, as.integer(m))
    set.seed(as.integer(seed))
    results <- vector("list", m)
    stratum <- if (is.null(strata)) rep("all", nrow(data)) else interaction(strata, drop = TRUE, lex.order = TRUE)

    for (iteration in seq_len(m)) {
        completed <- data
        for (name in variables) {
            x <- completed[[name]]
            missing <- is.na(x)
            if (is.character(x) || is.factor(x))
                missing <- missing | trimws(as.character(x)) == ""
            if (!any(missing))
                next
            if (all(missing))
                stop("Multiple imputation is impossible because ", name, " has no observed donor values.")
            for (level in unique(stratum[missing])) {
                recipients <- which(missing & stratum == level)
                donors <- which(!missing & stratum == level)
                if (length(donors) < 3L)
                    donors <- which(!missing)
                sampled <- sample(donors, length(recipients), replace = TRUE)
                if (is.factor(x)) {
                    text <- as.character(x)
                    text[recipients] <- text[sampled]
                    x <- factor(text, levels = levels(x))
                } else {
                    x[recipients] <- x[sampled]
                }
            }
            completed[[name]] <- x
        }
        results[[iteration]] <- completed
    }
    results
}

.te_rhs <- function(vars, response = NULL) {
    vars <- unique(vars)
    if (length(vars) == 0L) {
        if (is.null(response))
            return(~ 1)
        return(stats::as.formula(paste(response, "~ 1")))
    }
    stats::reformulate(vars, response = response)
}

.te_fit_binomial <- function(formula, data, method = "auto", label = "binomial model") {
    method <- match.arg(method, c("auto", "maximumLikelihood", "firth"))
    warnings <- character()
    fitML <- function() {
        withCallingHandlers(
            stats::glm(formula, data = data, family = stats::binomial()),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )
    }

    fit <- NULL
    needsFirth <- identical(method, "firth")
    if (!needsFirth) {
        fit <- tryCatch(fitML(), error = function(e) e)
        if (inherits(fit, "error")) {
            if (identical(method, "maximumLikelihood"))
                stop(label, " failed: ", conditionMessage(fit))
            needsFirth <- TRUE
        } else {
            fitted <- stats::fitted(fit)
            unstable <- !isTRUE(fit$converged) || any(!is.finite(stats::coef(fit))) ||
                any(abs(stats::coef(fit)[is.finite(stats::coef(fit))]) > 12) ||
                any(fitted < 1e-07 | fitted > 1 - 1e-07) ||
                any(grepl("did not converge|fitted probabilities numerically 0 or 1", warnings, ignore.case = TRUE))
            needsFirth <- identical(method, "auto") && unstable
        }
    }

    used <- "Maximum likelihood logistic regression"
    if (needsFirth) {
        if (!requireNamespace("brglm2", quietly = TRUE))
            stop(label, " requires Firth bias-reduced logistic regression, but the brglm2 dependency is unavailable.")
        warnings <- character()
        fit <- withCallingHandlers(
            stats::glm(
                formula,
                data = data,
                family = stats::binomial(),
                method = brglm2::brglmFit,
                type = "AS_mean"
            ),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )
        used <- "Firth mean bias-reduced logistic regression"
    }
    list(fit = fit, method = used, warnings = unique(warnings), firth = needsFirth)
}

.te_cluster_se <- function(influence, cluster = NULL) {
    influence <- as.numeric(influence)
    n <- length(influence)
    if (is.null(cluster))
        return(stats::sd(influence) / sqrt(n))
    cluster <- as.character(cluster)
    if (anyNA(cluster))
        stop("Cluster identifiers cannot be missing.")
    sums <- rowsum(influence, cluster, reorder = FALSE)[, 1]
    groups <- length(sums)
    if (groups < 2L)
        stop("Cluster-robust inference requires at least two clusters.")
    sqrt((groups / (groups - 1)) * sum(sums^2) / n^2)
}

.te_effect_rows <- function(mu0, mu1, if0, if1, method, confidenceLevel, cluster = NULL, binary = TRUE) {
    alpha <- 1 - confidenceLevel / 100
    critical <- stats::qnorm(1 - alpha / 2)
    make <- function(estimand, estimate, influence, transform = "identity") {
        se <- .te_cluster_se(influence, cluster)
        if (transform == "log") {
            lower <- exp(log(estimate) - critical * se)
            upper <- exp(log(estimate) + critical * se)
            statistic <- log(estimate) / se
            seScale <- "log"
        } else {
            lower <- estimate - critical * se
            upper <- estimate + critical * se
            statistic <- estimate / se
            seScale <- "estimate"
        }
        data.frame(
            method = method,
            estimand = estimand,
            estimate = estimate,
            se = se,
            seScale = seScale,
            lower = lower,
            upper = upper,
            statistic = statistic,
            p = 2 * stats::pnorm(-abs(statistic)),
            stringsAsFactors = FALSE
        )
    }

    rows <- rbind(
        make("Mean outcome: no treatment", mu0, if0),
        make("Mean outcome: treatment", mu1, if1),
        make(if (binary) "Risk difference" else "Mean difference", mu1 - mu0, if1 - if0)
    )
    if (binary && is.finite(mu0) && is.finite(mu1) && mu0 > 0 && mu1 > 0) {
        rows <- rbind(rows, make("Risk ratio", mu1 / mu0, if1 / mu1 - if0 / mu0, "log"))
    }
    if (binary && all(c(mu0, mu1) > 0 & c(mu0, mu1) < 1)) {
        oddsRatio <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
        rows <- rbind(rows, make("Odds ratio", oddsRatio, if1 / (mu1 * (1 - mu1)) - if0 / (mu0 * (1 - mu0)), "log"))
    }
    rows
}

.te_balance_table <- function(data, vars, A, weights, rawPS = NULL, labels = NULL) {
    vars <- unique(vars)
    if (length(vars) == 0L)
        return(data.frame(covariate = character(), unweighted = numeric(), weighted = numeric(), stringsAsFactors = FALSE))
    mm <- stats::model.matrix(.te_rhs(vars), data = data)
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    if (!is.null(rawPS)) {
        mm <- cbind(`Propensity score` = rawPS, mm)
    }
    if (ncol(mm) == 0L)
        return(data.frame(covariate = character(), unweighted = numeric(), weighted = numeric(), stringsAsFactors = FALSE))

    moments <- function(x, group, w) {
        keep <- A == group
        xx <- x[keep]
        ww <- w[keep]
        mean <- sum(ww * xx) / sum(ww)
        variance <- sum(ww * (xx - mean)^2) / sum(ww)
        c(mean = mean, variance = variance)
    }
    smd <- function(x, w) {
        one <- moments(x, 1, w)
        zero <- moments(x, 0, w)
        denominator <- sqrt((one[["variance"]] + zero[["variance"]]) / 2)
        if (!is.finite(denominator) || denominator == 0)
            return(0)
        (one[["mean"]] - zero[["mean"]]) / denominator
    }
    names <- colnames(mm)
    if (!is.null(labels)) {
        for (internal in names(labels))
            names <- sub(paste0("^", internal), labels[[internal]], names)
    }
    data.frame(
        covariate = names,
        unweighted = vapply(seq_len(ncol(mm)), function(j) smd(mm[, j], rep(1, length(A))), numeric(1)),
        weighted = vapply(seq_len(ncol(mm)), function(j) smd(mm[, j], weights), numeric(1)),
        stringsAsFactors = FALSE
    )
}

.te_fit_causal_once <- function(data, outcomeType = "binary", outcomeVars = character(), propensityVars = character(),
                                minPS = .01, maxPS = .99, outcomeMethod = "auto", propensityMethod = "auto",
                                outcomeSpecification = "pooled", confidenceLevel = 95, cluster = NULL) {
    A <- data$A
    Y <- data$Y
    psResult <- .te_fit_binomial(.te_rhs(propensityVars, "A"), data, propensityMethod, "Propensity model")
    rawPS <- as.numeric(stats::predict(psResult$fit, type = "response"))
    if (any(!is.finite(rawPS)))
        stop("The propensity model produced non-finite predictions.")
    ps <- pmin(pmax(rawPS, minPS), maxPS)

    outcomeWarnings <- character()
    outcomeMethods <- character()
    if (identical(outcomeSpecification, "pooled")) {
        formula <- .te_rhs(c("A", outcomeVars), "Y")
        if (identical(outcomeType, "binary")) {
            outcomeResult <- .te_fit_binomial(formula, data, outcomeMethod, "Outcome model")
            outcomeFit <- outcomeResult$fit
            outcomeWarnings <- outcomeResult$warnings
            outcomeMethods <- outcomeResult$method
        } else {
            outcomeFit <- stats::lm(formula, data = data)
            outcomeMethods <- "Ordinary least squares"
        }
        new0 <- data
        new1 <- data
        new0$A <- 0L
        new1$A <- 1L
        m0 <- as.numeric(stats::predict(outcomeFit, newdata = new0, type = "response"))
        m1 <- as.numeric(stats::predict(outcomeFit, newdata = new1, type = "response"))
    } else {
        fitArm <- function(arm) {
            armData <- data[A == arm, , drop = FALSE]
            formula <- .te_rhs(outcomeVars, "Y")
            if (identical(outcomeType, "binary"))
                .te_fit_binomial(formula, armData, outcomeMethod, paste0("Outcome model for treatment ", arm))
            else
                list(fit = stats::lm(formula, data = armData), method = "Ordinary least squares", warnings = character())
        }
        result0 <- fitArm(0L)
        result1 <- fitArm(1L)
        m0 <- as.numeric(stats::predict(result0$fit, newdata = data, type = "response"))
        m1 <- as.numeric(stats::predict(result1$fit, newdata = data, type = "response"))
        outcomeWarnings <- unique(c(result0$warnings, result1$warnings))
        outcomeMethods <- paste(unique(c(result0$method, result1$method)), collapse = "; ")
    }
    if (any(!is.finite(c(m0, m1))))
        stop("An outcome model produced non-finite predictions.")

    aipw1 <- m1 + A / ps * (Y - m1)
    aipw0 <- m0 + (1 - A) / (1 - ps) * (Y - m0)
    muA1 <- mean(aipw1)
    muA0 <- mean(aipw0)
    ifA1 <- aipw1 - muA1
    ifA0 <- aipw0 - muA0

    w1 <- A / ps
    w0 <- (1 - A) / (1 - ps)
    muI1 <- sum(w1 * Y) / sum(w1)
    muI0 <- sum(w0 * Y) / sum(w0)
    ifI1 <- (w1 * (Y - muI1)) / mean(w1)
    ifI0 <- (w0 * (Y - muI0)) / mean(w0)

    binary <- identical(outcomeType, "binary")
    effects <- rbind(
        .te_effect_rows(muA0, muA1, ifA0, ifA1, "AIPW", confidenceLevel, cluster, binary),
        .te_effect_rows(muI0, muI1, ifI0, ifI1, "IPTW (Hajek)", confidenceLevel, cluster, binary)
    )
    ipw <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
    balance <- .te_balance_table(data, propensityVars, A, ipw, rawPS)
    ess <- c(
        untreated = sum(ipw[A == 0])^2 / sum(ipw[A == 0]^2),
        treated = sum(ipw[A == 1])^2 / sum(ipw[A == 1]^2)
    )
    list(
        effects = effects,
        A = A,
        Y = Y,
        rawPS = rawPS,
        ps = ps,
        ipw = ipw,
        balance = balance,
        ess = ess,
        truncated = sum(rawPS < minPS | rawPS > maxPS),
        psMethod = psResult$method,
        outcomeMethod = outcomeMethods,
        warnings = unique(c(psResult$warnings, outcomeWarnings))
    )
}

.te_bootstrap_causal <- function(data, fitArgs, samples = 0L, seed = 24680L, clusterName = NULL) {
    samples <- as.integer(samples)
    if (samples <= 0L)
        return(NULL)
    set.seed(as.integer(seed))
    keys <- c("AIPW|Risk difference", "AIPW|Risk ratio", "AIPW|Odds ratio",
              "IPTW (Hajek)|Risk difference", "IPTW (Hajek)|Risk ratio", "IPTW (Hajek)|Odds ratio")
    draws <- matrix(NA_real_, nrow = samples, ncol = length(keys), dimnames = list(NULL, keys))
    if (!is.null(clusterName)) {
        clusters <- unique(as.character(data[[clusterName]]))
        if (length(clusters) < 2L)
            stop("Cluster bootstrap requires at least two clusters.")
    }
    for (iteration in seq_len(samples)) {
        if (is.null(clusterName)) {
            index <- sample.int(nrow(data), nrow(data), replace = TRUE)
            boot <- data[index, , drop = FALSE]
            cluster <- NULL
        } else {
            sampled <- sample(clusters, length(clusters), replace = TRUE)
            pieces <- lapply(seq_along(sampled), function(j) {
                piece <- data[as.character(data[[clusterName]]) == sampled[[j]], , drop = FALSE]
                piece$.bootstrapCluster <- paste0("b", j)
                piece
            })
            boot <- do.call(rbind, pieces)
            rownames(boot) <- NULL
            cluster <- boot$.bootstrapCluster
        }
        if (length(unique(boot$A)) < 2L)
            next
        args <- c(list(data = boot, cluster = cluster), fitArgs)
        fitted <- tryCatch(do.call(.te_fit_causal_once, args), error = function(e) NULL)
        if (is.null(fitted))
            next
        for (j in seq_len(nrow(fitted$effects))) {
            key <- paste(fitted$effects$method[[j]], fitted$effects$estimand[[j]], sep = "|")
            if (key %in% keys)
                draws[iteration, key] <- fitted$effects$estimate[[j]]
        }
    }
    draws
}

.te_apply_bootstrap_ci <- function(effects, draws, confidenceLevel) {
    if (is.null(draws))
        return(effects)
    probs <- c((1 - confidenceLevel / 100) / 2, 1 - (1 - confidenceLevel / 100) / 2)
    effects$bootstrapSuccess <- NA_integer_
    for (i in seq_len(nrow(effects))) {
        key <- paste(effects$method[[i]], effects$estimand[[i]], sep = "|")
        if (!key %in% colnames(draws))
            next
        values <- draws[, key]
        values <- values[is.finite(values)]
        effects$bootstrapSuccess[[i]] <- length(values)
        if (length(values) >= 20L) {
            interval <- stats::quantile(values, probs = probs, names = FALSE)
            effects$lower[[i]] <- interval[[1]]
            effects$upper[[i]] <- interval[[2]]
        }
    }
    effects
}

.te_pool_mi_effects <- function(fits, confidenceLevel) {
    template <- fits[[1]]$effects
    alpha <- 1 - confidenceLevel / 100
    for (i in seq_len(nrow(template))) {
        estimates <- vapply(fits, function(x) x$effects$estimate[[i]], numeric(1))
        variances <- vapply(fits, function(x) x$effects$se[[i]]^2, numeric(1))
        transform <- template$seScale[[i]]
        if (identical(transform, "log")) {
            estimates <- log(estimates)
            qbar <- mean(estimates)
            total <- mean(variances) + (1 + 1 / length(fits)) * stats::var(estimates)
            critical <- stats::qnorm(1 - alpha / 2)
            template$estimate[[i]] <- exp(qbar)
            template$se[[i]] <- sqrt(total)
            template$lower[[i]] <- exp(qbar - critical * sqrt(total))
            template$upper[[i]] <- exp(qbar + critical * sqrt(total))
            template$statistic[[i]] <- qbar / sqrt(total)
        } else {
            qbar <- mean(estimates)
            total <- mean(variances) + (1 + 1 / length(fits)) * stats::var(estimates)
            critical <- stats::qnorm(1 - alpha / 2)
            template$estimate[[i]] <- qbar
            template$se[[i]] <- sqrt(total)
            template$lower[[i]] <- qbar - critical * sqrt(total)
            template$upper[[i]] <- qbar + critical * sqrt(total)
            template$statistic[[i]] <- qbar / sqrt(total)
        }
        template$p[[i]] <- 2 * stats::pnorm(-abs(template$statistic[[i]]))
    }
    template
}

.te_parse_time <- function(x) {
    if (inherits(x, c("POSIXct", "POSIXlt")))
        return(as.POSIXct(x, tz = "UTC"))
    text <- .te_trim(x)
    suppressWarnings(as.POSIXct(
        as.character(text),
        tz = "UTC",
        tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d", "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M")
    ))
}

.te_elapsed <- function(value, index = NULL, mode = "elapsed", unit = "hours") {
    if (identical(mode, "elapsed")) {
        result <- suppressWarnings(as.numeric(.te_trim(value)))
    } else {
        if (is.null(index))
            stop("An index/start time is required when time mode is date-time.")
        result <- as.numeric(difftime(.te_parse_time(value), .te_parse_time(index), units = "hours"))
        if (identical(unit, "days"))
            result <- result / 24
        return(result)
    }
    result
}

.te_specification_fingerprint <- function(x) {
    digest::digest(x, algo = "sha256", serialize = TRUE)
}
