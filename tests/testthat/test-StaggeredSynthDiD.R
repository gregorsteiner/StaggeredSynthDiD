test_that("Function correctly returns list", {
  # simulate simple example dataset with staggered treatment
  set.seed(1)
  N <- 10 # number of units
  TT <- 10 # number of time periods
  tau <- 0.5 # actual treatment effect
  data <- data.frame("Unit" = rep(1:N, each = TT),
                     "Time" = rep(1:TT, N),
                     "W" = 0)
  data <- within(data, {
    # add staggered treatment (for units 1, 2, and 3 with different start timing)
    W[Unit == 1 & Time >= 3] <- 1
    W[Unit == 2 & Time >= 4] <- 1
    W[Unit == 3 & Time >= 5] <- 1
    # dependent variable is linear function of covariates and treatment
    Y <- tau * W + rnorm(N * TT)
  })

  # use function
  res <- StaggeredSynthDiD(data)

  # test if the output is a list with 3 elements
  expect_equal(length(res), 3)

})
