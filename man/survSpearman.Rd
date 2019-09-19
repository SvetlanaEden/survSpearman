\name{survSpearman}
\alias{survSpearman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Spearman's Correlation for Bivariate Survival Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes non-parametric estimate of Spearman's rank correlation for bivariate survival data. Two correlations are returned: a highest rank correlation that assigns observations beyond a specified region to have the highest rank value, and a restricted correlation that estimates the correlation within the specified region.
}
\usage{
survSpearman(bivarSurf, tauX = Inf, tauY = Inf)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bivarSurf}{
  A matrix containing the bivariate survival probability and the marginal survival probabilities. The first column of \code{"bivarSurf"} is the survival probability on one axis (we call it X). The first row of \code{"bivarSurf"} is the survival probability on the other axis (we call it Y). The row names of \code{"bivarSurf"} are ordered \code{"X"}'s. The column names of \code{"bivarSurf"} are ordered \code{"Y"}'s. Element \code{"bivarSurf[1,1]"} equals to 1. Its first row and column names are \code{"0"} and \code{"0"} (see examples below).
}
  \item{tauX}{
  A value that defines the restricted region for X coordinate.
}
  \item{tauY}{
  A value that defines the restricted region for Y coordinate.
}
}
\details{
  Returns a 2 x 1 matrix with two Spearman's correlation estimates. Both correlations are computed in the restricted region defined by values \code{"tauX"} and \code{"tauY"}: \code{"[0,tauX)x[0,tauY)"}. The highest rank Spearman's correlation is computed by assigning those events greater than tauX and/or tauY to have highest rank values. The restricted Spearman's correlation is computed as the rank correlation within the restricted region. When argument \code{"bivarSurf"} is computed using Dabrowska's method, the resulting bivariate survival probability can have negative mass. The problem of negative mass can result in zero or negative probability of failure in the restricted region, so in this case the restricted Spearman's correlation cannot be computed and NA value is returned. This only happens when the sample size is small and censoring is heavy.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  Returns a 2x1 matrix with two Spearman's correlation values, highest rank and restricted correlations, in the restricted region.
}
\references{
Eden, S., Li, C., Shepherd B. (2019). Non-parametric Estimation of Spearman's Rank Correlation with Bivariate Survival Data.
}
\author{
%%  ~~who you are~~
Svetlana Eden <svetlana.eden@vanderbilt.edu>
}

\examples{
X = c(0.5, 0.6, 0.7, 0.8)
Y = c(0.44, 0.77, 0.88, 0.99)
deltaX = c(1, 0, 1, 1)
deltaY = c(1, 1, 1, 1)

bivarSurf = survDabrowska(X, Y, deltaX, deltaY)$DabrowskaEst
bivarSurf
#         0  0.44 0.77 0.88 0.99
# 0   1.000 0.750 0.50 0.25    0
# 0.5 0.750 0.750 0.50 0.25    0
# 0.6 0.750 0.750 0.50 0.25    0
# 0.7 0.375 0.375 0.25 0.25    0
# 0.8 0.000 0.000 0.00 0.00    0

### When both last observed events are not censored and
### the restricted region is included the support, both correlations are the same
survSpearman(bivarSurf, tauX = Inf, tauY = Inf)

survSpearman(bivarSurf, tauX = Inf, tauY = 0.8)
#             Correlation
# HighestRank   0.6531973
# Restricted    0.9428090
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~bivariate }
\keyword{ ~survival }
\keyword{ ~Spearman }
\keyword{ ~correlation }
