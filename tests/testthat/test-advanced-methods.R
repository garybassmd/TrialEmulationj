test_that("AIPW point-outcome analysis runs", {
    set.seed(14)
    n <- 400
    x <- stats::rnorm(n)
    z <- stats::rbinom(n, 1, .4)
    ps <- stats::plogis(-.2 + .6 * x - .3 * z)
    treatment <- stats::rbinom(n, 1, ps)
    outcome <- stats::rbinom(n, 1, stats::plogis(-1 + .5 * treatment + .4 * x + .2 * z))
    data <- data.frame(treatment, outcome, x, z)

    result <- trialemulationj::aipw(
        data = data,
        treatment = "treatment",
        outcome = "outcome",
        outcomeType = "binary",
        outcomeCovs = c("x", "z"),
        propensityCovs = c("x", "z"),
        showOverlapPlot = FALSE
    )

    expect_s3_class(result, "aipwResults")
})

test_that("clone-censor-weight analysis runs with person-level bootstrap", {
    data("dummy_data", package = "survivalCCW")
    data <- survivalCCW::dummy_data[seq_len(80), ]

    result <- trialemulationj::cloneCensor(
        data = data,
        id = "id",
        event = "event",
        timeToEvent = "timetoevent",
        exposure = "exposure",
        timeToExposure = "timetoexposure",
        weightCovs = c("cov1", "cov2"),
        gracePeriod = 80,
        horizon = 200,
        bootstrapSamples = 20,
        showWeightPlot = FALSE
    )

    expect_s3_class(result, "cloneCensorResults")
})

test_that("doubly robust survival analysis runs", {
    data("survivalExample", package = "precmed")
    data <- precmed::survivalExample[seq_len(500), ]
    data$treatment01 <- as.integer(data$trt == "drug1")
    horizon <- with(data, min(
        stats::quantile(y[treatment01 == 1], .8),
        stats::quantile(y[treatment01 == 0], .8)
    ))

    result <- trialemulationj::drSurvival(
        data = data,
        time = "y",
        event = "d",
        treatment = "treatment01",
        outcomeCovs = c("age", "female"),
        propensityCovs = c("age", "previous_treatment"),
        censorCovs = c("age"),
        horizon = horizon,
        bootstrapSamples = 20,
        showBootstrapPlot = FALSE
    )

    expect_s3_class(result, "drSurvivalResults")
})
