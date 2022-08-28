

### estimate variance (bootstrap estimator) and construct CI ###
staggered_synthdid <- function(data, unit = "Unit", time = "Time",
                               outcome = "Y", treatment = "W",
                               covariates = NULL, vcov = "boot",
                               boot_iter = 100, ci_cov = 0.95){


  # select relevant columns and change column names
  data <- data[, c(unit, time, outcome, treatment, covariates)]
  colnames(data)[1:4] <- c("Unit", "Time", "Y", "W")

  # define parameters
  N <- length(unique(data$Unit))
  TT <- length(unique(data$Time))

  # compute the ATT as the weighted average
  tau_hat <- estimate(data, covariates = covariates)

  # bootstrap varaince estimator
  if(vcov == "boot"){
    B <- boot_iter
    taus_boot <- rep(NA, B)

    for (i in 1:B) {
      # sample units
      units_boot <- sample(unique(data$Unit), N, replace = TRUE)

      # only keep these units
      data_boot <- do.call(rbind, lapply(units_boot, function(x) data[data$Unit == x, ]))

      # rename units (such that the panel.matrices function works as intended)
      data_boot$Unit <- rep(1:N, each = TT)

      # if there is at least one treated unit, compute the sdid estimator
      if(any(data_boot$W == 1)) taus_boot[i] <- estimate(data_boot)

    }

    # compute variance
    V_hat <- mean((taus_boot - mean(taus_boot, na.rm = TRUE))^2, na.rm = TRUE)
  }

  # jackknife variance estimator
  if(vcov == "jack"){
    # loop over N units and exclude i-th unit
    taus_jack <- sapply(unique(data$Unit), function(i){
      # exclude unit i
      data_jack <- data[data$Unit != i, ]

      # rename units (such that the panel.matrices function works as intended)
      data_jack$Unit <- rep(1:(N-1), each = TT)

      # and compute tau_hat for the reduced dataset
      tau_hat <- estimate(data_jack)
      return(tau_hat)
    })

    # compute variance
    V_hat <- (N - 1) * mean((taus_jack - tau_hat)^2)
  }


  # construct CI
  alpha <- 1 - ci_cov
  CI <- tau_hat + c("Lower" = 1, "Upper" = -1) * qnorm(alpha / 2) * sqrt(V_hat)

  # return everything
  return(list("ATT" = tau_hat,
              "Bootstrap SE" = sqrt(V_hat),
              "Confidence Interval" = CI))

}
