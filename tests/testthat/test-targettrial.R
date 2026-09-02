test_that("ITT analysis runs on the package example", {
    data <- TrialEmulation::data_censored

    result <- trialemulationj::targetTrial(
        data = data,
        id = "id",
        period = "period",
        treatment = "treatment",
        outcome = "outcome",
        eligible = "eligible",
        estimand = "itt",
        outcomeCovs = c("age_s", "x1"),
        useCensor = FALSE,
        censor = NULL,
        maxFollowUp = 3,
        confInt = FALSE,
        truncateWeights = FALSE,
        chunkSize = 100,
        showWeightPlot = FALSE
    )

    expect_s3_class(result, "targetTrialResults")
})

test_that("per-protocol analysis with censoring runs", {
    data <- TrialEmulation::data_censored

    result <- trialemulationj::targetTrial(
        data = data,
        id = "id",
        period = "period",
        treatment = "treatment",
        outcome = "outcome",
        eligible = "eligible",
        estimand = "pp",
        outcomeCovs = c("age_s", "x1"),
        switchNumCovs = c("age_s"),
        switchDenCovs = c("age_s", "x1", "x3"),
        useCensor = TRUE,
        censor = "censored",
        censorNumCovs = c("age_s"),
        censorDenCovs = c("age_s", "x1", "x2"),
        maxFollowUp = 3,
        confInt = FALSE,
        truncateWeights = TRUE,
        truncatePercentile = 99,
        chunkSize = 100,
        showWeightPlot = FALSE
    )

    expect_s3_class(result, "targetTrialResults")
})
