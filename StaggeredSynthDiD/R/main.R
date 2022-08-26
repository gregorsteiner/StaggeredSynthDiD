

### function to estimate the ATT ###
estimate_att <- function(data, unit = "Unit", time = "Time",
                         outcome = "Y", treatment = "W",
                         covariates = NULL){
  
  # select relevant columns and change column names
  data <- data[, c(unit, time, outcome, treatment, covariates)]
  colnames(data)[1:4] <- c("Unit", "Time", "Y", "W")

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
  att <- as.numeric(result["Estimate", ] %*% result["Weight", ])
  return(att)
}

### estimate variance (bootstrap estimator) and construct CI ###
synthdid_staggered <- function(data, unit = "Unit", time = "Time",
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
  tau_hat <- estimate_att(data, covariates = covariates)
  
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
    if(any(data_boot$W == 1)) taus[i] <- estimate_att(data_boot)
    
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

