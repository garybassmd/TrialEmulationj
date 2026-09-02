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

test_that("blank-as-absence reuses an existing labelled negative level", {
    set.seed(2203)
    n <- 180L
    data <- data.frame(
        treatment = rep(c("No", "Yes"), each = n / 2L),
        death = rep(c("alive at discharge", "alive at discharge", "died"), length.out = n),
        age = stats::rnorm(n, 65, 11),
        stringsAsFactors = FALSE
    )
    data$death[seq(9, n, by = 18)] <- NA_character_
    result <- trialemulationj::aipw(
        data = data,
        treatment = "treatment",
        outcome = "death",
        outcomeCovs = "age",
        propensityCovs = "age",
        absenceVars = "death",
        sensitivityGrid = FALSE,
        showOverlapPlot = FALSE,
        showBalancePlot = FALSE,
        showEffectPlot = FALSE
    )
    outcomeCoding <- result$summary$asDF$value[result$summary$asDF$measure == "Outcome coding"]
    expect_match(outcomeCoding, "alive at discharge = 0")
    expect_equal(result$missingAudit$asDF$recoded[result$missingAudit$asDF$variable == "death"], 10L)
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
        treatmentALabel = "NGT",
        treatmentBLabel = "WSEC",
        showFlowPlot = FALSE,
        showBalancePlot = FALSE,
        showEffectPlot = FALSE
    )
    expect_equal(length(unique(result$flow$asDF$window)), 2L)
    expect_equal(nrow(result$strategies$asDF), 8L)
    expect_equal(nrow(result$effects$asDF), 8L)
    expect_true(all(is.finite(result$effects$asDF$estimate)))
    expect_true(sum(result$timeline$asDF$n) >= 3L)
    expect_true(all(grepl("^NGT ", result$effects$asDF$stratum)))
    expect_true(all(grepl("^(NGT|WSEC) ", c(result$strategies$asDF$treatmentA, result$strategies$asDF$treatmentB))))
})

test_that("landmark workflow validates coding after eligibility and audits observation time", {
    nEligible <- 160L
    nIneligible <- 40L
    data <- data.frame(
        eligible = c(rep("Yes", nEligible), rep("No", nIneligible)),
        treatmentA = c(rep(c("No", "Yes"), each = nEligible / 2L), rep("Not assessed", nIneligible)),
        treatmentB = c(rep(c("No", "Yes"), length.out = nEligible), rep("Not applicable", nIneligible)),
        outcome = c(rep(c("alive at discharge", "alive at discharge", "died", "alive at discharge"), length.out = nEligible), rep("Direct surgery", nIneligible)),
        dischargeTime = c(rep(80, nEligible), rep(NA_real_, nIneligible)),
        age = c(seq(45, 84, length.out = nEligible), rep(60, nIneligible)),
        site = c(rep(sprintf("site%02d", 1:20), length.out = nEligible), rep(NA_character_, nIneligible)),
        stringsAsFactors = FALSE
    )
    data$site[[1]] <- NA_character_
    data$outcome[[1]] <- NA_character_
    data$dischargeTime[[1]] <- NA_real_
    data$dischargeTime[[2]] <- -1
    data$dischargeTime[[3]] <- 12

    excluded <- trialemulationj::landmarkTrial(
        data = data,
        eligibility = "eligible",
        treatmentA = "treatmentA",
        treatmentB = "treatmentB",
        outcome = "outcome",
        cluster = "site",
        dischargeTime = "dischargeTime",
        covariates = "age",
        eligibilityPositive = "Yes",
        treatmentAPositive = "Yes",
        treatmentBPositive = "Yes",
        outcomePositive = "died",
        absenceVars = "outcome",
        missingCluster = "unknown",
        unknownDischarge = "exclude",
        runSecondWindow = FALSE,
        showFlowPlot = FALSE,
        showBalancePlot = FALSE,
        showEffectPlot = FALSE
    )
    expect_equal(excluded$flow$asDF$n[excluded$flow$asDF$step == "Meet baseline eligibility"], nEligible)
    expect_equal(excluded$flow$asDF$n[excluded$flow$asDF$step == "Remain observed and event-free at landmark"], nEligible - 3L)
    expect_equal(excluded$timeline$asDF$n[excluded$timeline$asDF$issue == "Missing/blank cluster ID"], 1L)
    expect_equal(excluded$timeline$asDF$n[excluded$timeline$asDF$issue == "Missing/unparseable discharge time"], 1L)
    expect_equal(excluded$timeline$asDF$n[excluded$timeline$asDF$issue == "Negative discharge time"], 1L)
    expect_equal(excluded$timeline$asDF$n[excluded$timeline$asDF$issue == "Missing/blank value explicitly recoded as absence"], 1L)

    retained <- trialemulationj::landmarkTrial(
        data = data,
        eligibility = "eligible",
        treatmentA = "treatmentA",
        treatmentB = "treatmentB",
        outcome = "outcome",
        cluster = "site",
        dischargeTime = "dischargeTime",
        covariates = "age",
        eligibilityPositive = "Yes",
        treatmentAPositive = "Yes",
        treatmentBPositive = "Yes",
        outcomePositive = "died",
        absenceVars = "outcome",
        missingCluster = "unknown",
        unknownDischarge = "retain",
        runSecondWindow = FALSE,
        showFlowPlot = FALSE,
        showBalancePlot = FALSE,
        showEffectPlot = FALSE
    )
    expect_equal(retained$flow$asDF$n[retained$flow$asDF$step == "Remain observed and event-free at landmark"], nEligible - 1L)
})
