\name{StaggeredSynthDiD}
\alias{StaggeredSynthDiD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Synthetic DiD estimation for staggered adoption settings
%%  ~~function to do ... ~~
}
\description{Provides a point estimate of and inference on the Average Treatment Effect (ATE) in staggered adoption settings.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
StaggeredSynthDiD(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The dataset used. It Should be in long format.}
  \item{unit}{The name of the column containing the unit id.}
  \item{time}{The name of the column containing the time period.}
  \item{outcome}{The name of the column containing the dependent variable.}
  \item{treatment}{The name of the column containing the treatment indicator.}
  \item{covariates}{A character vector of covariate column names. NULL if no covariates should be included}
  \item{vcov}{The variance estimator to be used. Should be one of "boot", "jack" or "placebo".}
  \item{iterations}{Number of resampling iterations. Only relevant for the bootstrap and placebo variance estimators.}
  \item{ci_cov}{Desired nominal coverage probability of the confidence interval.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Estimate}{Estimate of the ATE.}
  \item{SE}{Standard error of the point estimate.}
  \item{CI}{Confidence Interval with user-specified nominal coverage.}
}
\references{
Arkhangelsky, Dmitry, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. 2021. "Synthetic Difference-in-Differences." American Economic Review, 111 (12): 4088-4118.
}
\author{
Gregor Steiner
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# simulate simple example dataset
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
StaggeredSynthDiD(data, vcov = "jack")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
