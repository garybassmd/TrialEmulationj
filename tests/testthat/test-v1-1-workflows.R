make_landmark_data <- function(n = 640L, seed = 1101L) {
    set.seed(seed)
    site <- sprintf("site%02d", sample(seq_len(32L), n, replace = TRUE))
    age <- stats::rnorm(n, 65, 12)
    comorbidity <- stats::rbinom(n, 1, .35)
    treatmentA <- stats::rbinom(n, 1, stats::plogis(-.1 + .01 * (age - 65) + .4 * comorbidity))
    treatmentB <- stats::rbinom(n, 1, stats::plogis(-.4 + .3 * treatmentA + .012 * (age - 65) + .35 * comorbidity))
    treatmentATime <- ifelse(treatmentA == 1, stats::rgamma(n, 2, .12), NA_real_)
    treatmentBTime <- ifelse(treatmentB == 1, stats::rgamma(n, 2, .08), NA_real_)
    failure <- stats::rbinom(n, 1, stats::plogis(-2 + .35 * treatmentB - .3 * treatmentA + .018 * (age - 65) + .5 * comorbidity))
    failureTime <- ifelse(failure == 1, stats::runif(n, 8, 120), NA_real_)
    treatmentBTime[which(treatmentB == 1)[1:2]] <- NA_real_
    failureTime[which(failure == 1)[1]] <- -2
    data.frame(
        eligible = ifelse(stats::runif(n) < .96, "Yes", "No"),
        treatmentA = ifelse(treatmentA == 1, "Yes", "No"),
        treatmentB = ifelse(treatmentB == 1, "Yes", "No"),
        failure = ifelse(failure == 1, "Yes", "No"),
        treatmentATime = treatmentATime,
        treatmentBTime = treatmentBTime,
        failureTime = failureTime,
        dischargeTime = stats::runif(n, 55, 180),
        age = age,
        comorbidity = factor(comorbidity),
        site = site
    )
}

test_that("target-trial protocol builder produces a complete design record", {
    result <- trialemulationj::protocol(
        studyTitle = "Synthetic target trial",
        population = "Eligible admissions",
        eligibility = "Eligibility assessed at admission",
        strategies = "Treatment within 24 hours versus no treatment within 24 hours",
        timeZero = "Admission",
        gracePeriod = "24 hours",
        followUp = "From the 24-hour landmark to discharge",
        outcome = "Failure by discharge",
        competingEvents = "Death treated as a competing event",
        confoundingPlan = "Prespecified baseline DAG adjustment set",
        analysisPlan = "AIPW with whole-site bootstrap",
        missingPlan = "Multiple imputation",
        sensitivityPlan = "Repeat at 48 hours",
        acknowledgesAlignment = TRUE,
        acknowledgesFeasible = TRUE,
        acknowledgesBaseline = TRUE,
        acknowledgesPositivity = TRUE
    )
    expect_equal(nrow(result$specification$asDF), 14L)
    expect_equal(nrow(result$readiness$asDF), 9L)
    expect_true(all(result$readiness$asDF$status == "Ready"))
})

test_that("DAG workflow identifies and checks the baseline adjustment set", {
    data <- data.frame(exposure = 0:1, outcome = 1:0, age = c(50, 70), severity = c(0, 1))
    result <- trialemulationj::causalDAG(
        data = data,
        exposure = "exposure",
        outcome = "outcome",
        selectedAdjusters = c("age", "severity"),
        dagText = paste(
            "dag {",
            "age -> exposure; age -> outcome;",
            "severity -> exposure; severity -> outcome;",
            "exposure -> outcome",
            "}"
        ),
        showPlot = FALSE
    )
    expect_true(any(grepl("age", result$adjustmentSets$asDF$variables)))
    expect_true(any(grepl("severity", result$adjustmentSets$asDF$variables)))
    expect_equal(result$summary$asDF$value[result$summary$asDF$measure == "Selected set valid"], "Yes")
})

test_that("Firth models return finite AIPW estimates under complete separation", {
    set.seed(2202)
    n <- 240L
    data <- data.frame(
        treatment = rep(c("No", "Yes"), each = n / 2),
        outcome = rep(c("No", "Yes"), each = n / 2),
        age = stats::rnorm(n),
        site = rep(sprintf("site%02d", 1:24), each = 10)
    )
    result <- trialemulationj::aipw(
        data = data,
        treatment = "treatment",
        outcome = "outcome",
        outcomeCovs = "age",
        propensityCovs = "age",
        cluster = "site",
        outcomeMethod = "firth",
        propensityMethod = "firth",
        sensitivityGrid = FALSE,
        showOverlapPlot = FALSE,
        showBalancePlot = FALSE,
        showEffectPlot = FALSE
    )
    expect_true(all(is.finite(result$effects$asDF$estimate)))
    expect_true(all(grepl("Firth", result$summary$asDF$value[result$summary$asDF$measure %in% c("Propensity model", "Outcome model")])))
})

test_that("landmark workflow constructs 24- and 48-hour clustered analyses", {
    data <- make_landmark_data()
    result <- trialemulationj::landmarkTrial(
        data = data,
        eligibility = "eligible",
        treatmentA = "treatmentA",
        treatmentB = "treatmentB",
        outcome = "failure",
        cluster = "site",
        treatmentATime = "treatmentATime",
        treatmentBTime = "treatmentBTime",
        failureTime1 = "failureTime",
        dischargeTime = "dischargeTime",
        covariates = c("age", "comorbidity"),
        eligibilityPositive = "Yes",
        treatmentAPositive = "Yes",
        treatmentBPositive = "Yes",
        outcomePositive = "Yes",
        showFlowPlot = FALSE,
        showBalancePlot = FALSE,
        showEffectPlot = FALSE
    )
    expect_equal(length(unique(result$flow$asDF$window)), 2L)
    expect_equal(nrow(result$strategies$asDF), 8L)
    expect_equal(nrow(result$effects$asDF), 8L)
    expect_true(all(is.finite(result$effects$asDF$estimate)))
    expect_true(sum(result$timeline$asDF$n) >= 3L)
})
