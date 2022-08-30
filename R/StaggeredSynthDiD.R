

### get point estimate, estimate variance and construct CI ###
StaggeredSynthDiD <- function(data, unit = "Unit", time = "Time",
                               outcome = "Y", treatment = "W",
                               covariates = NULL, vcov = "jack",
                               iterations = 100, ci_cov = 0.95){


  # select relevant columns and change column names
  data <- data[, c(unit, time, outcome, treatment, covariates)]
  colnames(data)[1:4] <- c("Unit", "Time", "Y", "W")

  # define parameters
  N <- length(unique(data$Unit))
  TT <- length(unique(data$Time))
  N_tr <- length(unique(data[data$W == 1, "Unit"]))
  N_co <- N - N_tr

  # some checks
  if(!(vcov %in% c("boot", "jack", "placebo"))) stop("Invalid vcov argument")

  # define estimation function (will be reused for the variance estimation)
  estimate <- function(data, covariates = NULL){

    # create a new variable indicating when treatment started (0 indicates never treated)
    # we do this using data.table functionality
    data <- data.table::setDT(data)
    data[, Start := ifelse(any(W == 1), Time[W == 1][1], 0), by = .(Unit)]

    # isolate never treated units
    data_never <- data[data$Start == 0, ]

    # get the total number of treated unit-time combinations (needed for the weights calculation later)
    NT_treated <- sum(data$W == 1)

    # then loop over different adoption dates and compute treatment effect for each
    adop_dates <- unique(data$Start[data$Start != 0])

    result <- sapply(adop_dates, function(adop){

      # merge never treated and treated for this specific adoption date
      data_int <- rbind(data_never, data[data$Start == adop, ])

      # get matrices for this dataset
      input <- suppressWarnings(synthdid::panel.matrices(data_int, unit = unit, time = time,
                                                         outcome = outcome, treatment = treatment))

      # add covariates if desired
      if(!is.null(covariates)){
        # Convert them into 3d array as required by the synthdid function
        X <- simplify2array(lapply(covariates, function(x){
          do.call(rbind, split(data_int[, get(x)], data_int[, get(unit)]))
        }))

        # and compute estimator
        est <- synthdid::synthdid_estimate(input$Y, input$N0, input$T0, X)
      } else{
        # else compute estimator without covariates
        est <- synthdid::synthdid_estimate(input$Y, input$N0, input$T0)
      }

      # and compute weights
      weight <- sum(data_int$W == 1) / NT_treated

      # return estimated Treatment effect and weight
      return(c("Estimate" = as.numeric(est),
               "Weight" = weight))

    })

    # compute the ATT as the weighted average
    tau_hat <- as.numeric(result["Estimate", ] %*% result["Weight", ])
    return(tau_hat)
  }

  # compute the ATE using that function
  tau_hat <- estimate(data = data, covariates = covariates)

  # bootstrap varaince estimator
  if(vcov == "boot"){
    B <- iterations
    taus_boot <- rep(NA, B)

    for (i in 1:B) {
      # sample units
      units_boot <- sample(unique(data$Unit), N, replace = TRUE)

      # only keep these units
      data_boot <- do.call(rbind, lapply(units_boot, function(x) data[data$Unit == x, ]))

      # rename units (such that the panel.matrices function works as intended)
      data_boot$Unit <- rep(1:N, each = TT)

      # if there is at least one treated unit, compute the sdid estimator
      if(any(data_boot$W == 1)) taus_boot[i] <- estimate(data_boot, covariates = covariates)

    }

    # compute variance
    V_hat <- mean((taus_boot - mean(taus_boot, na.rm = TRUE))^2, na.rm = TRUE)
  }

  # jackknife variance estimator
  if(vcov == "jack"){
    # check if there's more than one treated unit (otherwise the jackknife estimator is not defined)
    if(N_tr == 1) stop("The jackknife variance estimator requires at least two treated units")

    # loop over N units and exclude i-th unit
    taus_jack <- sapply(unique(data$Unit), function(i){
      # exclude unit i
      data_jack <- data[data$Unit != i, ]

      # rename units (such that the panel.matrices function works as intended)
      data_jack$Unit <- rep(1:(N-1), each = TT)

      # and compute tau_hat for the reduced dataset
      tau_hat <- estimate(data_jack, covariates = covariates)
      return(tau_hat)
    })

    # compute variance
    V_hat <- (N - 1) * mean((taus_jack - tau_hat)^2)
  }

  # placebo variance estimator
  if(vcov == "placebo"){
    B <- iterations
    taus_plac <- rep(NA, B)

    # isolate control units
    co <- data[data$Time == max(data$Time, na.rm = TRUE) & data$W == 0, "Unit"]

    for (i in 1:B) {
      # sample control units to become treated
      tr <- sample(co, N_tr)

      # and randomize treatment adoption timing
      adop_timing <- sample(sort(unique(data$Time))[3:length(unique(data$Time))], N_tr)

      # create placebo data
      data_tr <- do.call(rbind, Map(function(unit, time){
        # select unit
        data_int <- data[data$Unit == unit, ]

        # add treatment
        data_int[data_int$Time >= time, "W"] <- 1

        return(data_int)

      }, tr, adop_timing))
      # add placebo control data
      data_co <- data[data$Unit %in% co & !(data$Unit %in% tr), ]
      data_plac <- rbind(data_co, data_tr)

      # rename units (such that the panel.matrices function works as intended)
      data_plac$Unit <- rep(1:N_co, each = TT)

      # estimate tau
      taus_plac[i] <- estimate(data_plac, covariates = covariates)
    }

    # compute variance
    V_hat <- mean((taus_plac - mean(taus_plac))^2)
  }

  # construct CI
  alpha <- 1 - ci_cov
  CI <- tau_hat + c("Lower" = 1, "Upper" = -1) * qnorm(alpha / 2) * sqrt(V_hat)

  # return everything
  return(list("Estimate" = tau_hat,
              "SE" = sqrt(V_hat),
              "CI" = CI))

}


