# StaggeredSynthDiD: Synthetic Difference-in-Differences Estimation in Staggered Adoption Settings

This package extends the synthetic Difference-in-Differences (DiD) estimator by Arkhangelsky et al. (2021) to staggered adoption settings, building on the **synthdid** package by David Hirshberg.
As proposed in the Appendix of Arkhangelsky et al. (2021), the estimator is applied repeatedly, once for every adoption date, and a weighted average is calculated.

### Installation

The current development version can be installed from Github:

```R
devtools::install_github("gregorsteiner/StaggeredSynthDiD")
```


### References
Arkhangelsky, Dmitry, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. 2021. "Synthetic Difference-in-Differences." American Economic Review, 111 (12): 4088-4118.
