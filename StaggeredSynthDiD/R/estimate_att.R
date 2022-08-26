
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
