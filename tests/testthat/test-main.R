library(wired)

testthat::test_that("support_wired returns callable samplers and correct shapes", {
  testthat::skip_if_not_installed("MASS")
  testthat::skip_if_not_installed("mc2d")
  testthat::skip_if_not_installed("forecast")
  testthat::skip_if_not_installed("quantreg")

  set.seed(1)

  n <- 120
  ts_set <- data.frame(
    A = 100 + cumsum(rnorm(n, 0, 1)),
    B =  80 + cumsum(rnorm(n, 0, 1))
  )

  fit <- wired:::support_wired(
    ts_set = ts_set,
    future = 1,
    mode = "additive",
    n_testing = 2,
    dep_metric = "spearman",
    corr_adapt = "rolling",
    roll_window = 50,
    copula = "gaussian",
    seed = 123,

    # speed knobs passed down to probabilistic_mixer() via ...
    n_crps_mc = 100,
    q_grid_size = 100
  )

  testthat::expect_true(is.list(fit))
  testthat::expect_true(is.function(fit$rfun_trafo))
  testthat::expect_true(is.function(fit$rfun_level))
  testthat::expect_true(is.function(fit$rfun_both))

  xT <- fit$rfun_trafo(3)
  xL <- fit$rfun_level(3)
  xb <- fit$rfun_both(3)

  testthat::expect_true(is.matrix(xT))
  testthat::expect_true(is.matrix(xL))
  testthat::expect_equal(dim(xT), c(3, 2))
  testthat::expect_equal(dim(xL), c(3, 2))

  testthat::expect_true(is.list(xb))
  testthat::expect_true(all(c("trafo", "level") %in% names(xb)))
  testthat::expect_equal(dim(xb$trafo), c(3, 2))
  testthat::expect_equal(dim(xb$level), c(3, 2))
})

testthat::test_that("wired produces [n x p x H] arrays with expected dimnames", {
  testthat::skip_if_not_installed("MASS")
  testthat::skip_if_not_installed("mc2d")
  testthat::skip_if_not_installed("forecast")
  testthat::skip_if_not_installed("quantreg")

  set.seed(2)

  n <- 150
  ts_set <- data.frame(
    A = 120 + cumsum(rnorm(n, 0, 1)),
    B =  90 + cumsum(rnorm(n, 0, 1)),
    C =  70 + cumsum(rnorm(n, 0, 1))
  )

  fitH <- wired(
    ts_set = ts_set,
    future = 3,
    mode = "additive",
    n_testing = 2,
    dep_metric = "kendall",
    corr_adapt = "static",
    copula = "gaussian",
    seed = 123,

    # speed knobs
    n_crps_mc = 50,
    q_grid_size = 100
  )

  testthat::expect_true(is.list(fitH))
  testthat::expect_true(is.function(fitH$rfun_trafo))
  testthat::expect_true(is.function(fitH$rfun_level))
  testthat::expect_true(is.function(fitH$rfun_both))

  arrT <- fitH$rfun_trafo(4)
  arrL <- fitH$rfun_level(4)

  testthat::expect_true(is.array(arrT))
  testthat::expect_true(is.array(arrL))
  testthat::expect_equal(dim(arrT), c(3, 4, 3))
  testthat::expect_equal(dim(arrL), c(3, 4, 3))

  # dimnames: series names + horizon labels
  testthat::expect_equal(dimnames(arrT)[[1]], c("h1", "h2", "h3"))
  testthat::expect_equal(dimnames(arrT)[[3]],  c("A", "B", "C"))
})

testthat::test_that("wired rfun_both is consistent with rfun_trafo and rfun_level", {
  testthat::skip_if_not_installed("MASS")
  testthat::skip_if_not_installed("mc2d")
  testthat::skip_if_not_installed("forecast")
  testthat::skip_if_not_installed("quantreg")

  set.seed(3)

  n <- 130
  ts_set <- data.frame(
    A = 100 + cumsum(rnorm(n, 0, 1)),
    B =  95 + cumsum(rnorm(n, 0, 1))
  )

  fitH <- wired(
    ts_set = ts_set,
    future = 2,
    mode = "additive",
    n_testing = 2,
    dep_metric = "spearman",
    corr_adapt = "rolling",
    roll_window = 50,
    copula = "t",
    t_df = 7,
    seed = 321,

    # speed knobs
    n_crps_mc = 50,
    q_grid_size = 50
  )

  n_draw <- 5
  both <- fitH$rfun_both(n_draw)
  arrT <- fitH$rfun_trafo(n_draw)
  arrL <- fitH$rfun_level(n_draw)

  testthat::expect_true(is.list(both))
  testthat::expect_true(all(c("trafo", "level") %in% names(both)))
  testthat::expect_equal(dim(both$trafo), dim(arrT))
  testthat::expect_equal(dim(both$level), dim(arrL))

  # With your design (set.seed inside each sampler), these should match exactly.
  testthat::expect_equal(both$trafo, arrT)
  testthat::expect_equal(both$level, arrL)
})
