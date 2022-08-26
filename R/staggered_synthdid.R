

### estimate variance (bootstrap estimator) and construct CI ###
staggered_synthdid <- function(data, unit = "Unit", time = "Time",
                               outcome = "Y", treatment = "W",
                               covariates = NULL,
                               boot_iter = 100, ci_cov = 0.95){


  # select relevant columns and change column names
  data <- data[, c(unit, time, outcome, treatment, covariates)]
  colnames(data)[1:4] <- c("Unit", "Time", "Y", "W")

  # define parameters
  N <- length(unique(data$Unit))
  TT <- length(unique(data$Time))

  # compute the ATT as the weighted average
  tau_hat <- estimate(data, covariates = covariates)

  # estimate variance with the bootstrap estimator
  B <- boot_iter
  taus <- rep(NA, B)

  for (i in 1:B) {
    # sample units
    units_boot <- sample(unique(data$Unit), N, replace = TRUE)

    # only keep these units
    data_boot <- do.call(rbind, lapply(units_boot, function(x) data[data$Unit == x, ]))

    # rename units (such that the panel.matrices function works as intended)
    data_boot$Unit <- rep(1:N, each = TT)

    # if there is at least one treated unit, compute the sdid estimator
    if(any(data_boot$W == 1)) taus[i] <- estimate(data_boot)

  }

  # compute variance
  var_boot <- mean((taus - mean(taus, na.rm = TRUE))^2, na.rm = TRUE)

  # construct CI
  alpha <- 1 - ci_cov
  CI <- tau_hat + c(1, -1) * qnorm(alpha / 2) * sqrt(var_boot)

  # and add names
  names(CI) <- c("Lower", "Upper")

  # return everything
  return(list("ATT" = tau_hat,
              "Bootstrap SE" = sqrt(var_boot),
              "Confidence Interval" = CI))

}
