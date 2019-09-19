\name{survDabrowska}
\alias{survDabrowska}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bivariate survival estimation using Dabrowska's method.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function computes bivariate survival probability for right censored data using Dabrowska's [1988] method.
}
\usage{
survDabrowska(X, Y, deltaX, deltaY)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  Time to event or censoring for variable \code{"X"}. It indicates time to event if argument \code{"deltaX"} is 1 and time to censoring if argument \code{"deltaX"} is 0.}
  \item{Y}{
  Time to event or censoring for variable \code{"Y"}. It indicates time to event if argument \code{"deltaY"} is 1 and time to censoring if argument \code{"deltaY"} is 0.}
  \item{deltaX}{
  Event indicator for variable \code{"X"}. \code{"deltaX"} is 1 if the event was observed and \code{"deltaX"} is 0 if it was censored.}
  \item{deltaY}{
  Event indicator for variable \code{"Y"}. \code{"deltaY"} is 1 if the event was observed and \code{"deltaY"} is 0 if it was censored.}
}
\details{
  The function returns a list with two values: Dabrowska's estimator of bivariate survival probability \code{"DabrowskaEst"} and its corresponding cumulative distribution function (CDF) \code{"DabrowskaCDF"}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{DabrowskaEst}{A matrix containing the bivariate survival probability and the marginal survival probabilities. The first column of \code{"DabrowskaEst"} is the survival probability for right censored event corresponding to \code{"X"}. The first row of \code{"DabrowskaEst"} is the survival probability  for right censored event corresponding to \code{"Y"}. The row names of \code{"DabrowskaEst"} are ordered \code{"X"}'s. The column names of \code{"DabrowskaEst"} are ordered \code{"Y"}'s. Element \code{"DabrowskaEst[1,1]"} equals 1. Its first row and column names are \code{"0"} and \code{"0"} (see examples below).}
  \item{DabrowskaCDF}{A matrix containing the bivariate cumulative distribution function (CDF) and the marginal CDFs. The first column of \code{"DabrowskaCDF"} is the marginal CDF for \code{"X"}. The first row of \code{"DabrowskaCDF"} is the marginal CDF for \code{"Y"}. The row and column names of \code{"DabrowskaCDF"} are the same as for \code{"DabrowskaEst"}.}
}
\references{
Dabrowska, D. M. (1988). Kaplan–Meier estimate on the plane. The Annals of Statistics 16, 1475–1489.
}
\author{
%%  ~~who you are~~
Svetlana Eden <svetlana.eden@vanderbilt.edu>
Joey Sherrill <joey.sherrill@vumc.org>
}

\examples{

X = c(0.5, 0.6, 0.6, 0.8)
Y = c(0.44, 0.77, 0.88, 0.99)
deltaX = c(1, 0, 1, 1)
deltaY = c(1, 1, 1, 1)

survDabrowska(X, Y, deltaX, deltaY)
# $DabrowskaEst
#        0 0.44 0.77 0.88 0.99
# 0   1.00 0.75 0.50 0.25    0
# 0.5 0.75 0.75 0.50 0.25    0
# 0.6 0.50 0.50 0.25 0.25    0
# 0.8 0.00 0.00 0.00 0.00    0
#
# $DabrowskaCDF
#     0 0.44 0.77 0.88 0.99
# 0   0 0.00 0.00 0.00 0.00
# 0.5 0 0.25 0.25 0.25 0.25
# 0.6 0 0.25 0.25 0.50 0.50
# 0.8 0 0.25 0.50 0.75 1.00
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~bivariate }
\keyword{ ~survival }
\keyword{ ~Dabrowska }
