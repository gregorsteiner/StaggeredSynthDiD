# StaggeredSynthDiD: Synthetic Difference-in-Differences Estimation in Staggered Adoption Settings

This package extends the synthetic Difference-in-Differences (DiD) estimator by Arkhangelsky et al. (2021) to staggered adoption settings, building on the [**synthdid**](https://synth-inference.github.io/synthdid/) package.
As proposed in the Appendix of Arkhangelsky et al. (2021), the estimator is applied repeatedly, once for every adoption date, and a weighted average is calculated.

### Installation

The current development version can be installed from Github:

```R
devtools::install_github("gregorsteiner/StaggeredSynthDiD")
```


### Example

```R
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
  # dependent variable
  Y <- tau * W + rnorm(N * TT)
})

# use function
StaggeredSynthDiD(data)
´´´


### References
Arkhangelsky, Dmitry, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. 2021. "Synthetic Difference-in-Differences." American Economic Review, 111 (12): 4088-4118.
