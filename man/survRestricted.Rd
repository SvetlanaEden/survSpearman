\name{survRestricted}
\alias{survRestricted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes conditional bivariate survival probability in a restricted region.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function computes the conditional bivariate survival probability in a restricted region. In bivariate survival analysis we sometimes need to know the survival probability conditionally on surviving in a restricted region. This region is defined as \code{"[0,XLimit)x[0,YLimit)"}. The function takes an already computed survival surface and arguments that define the restricted region \code{"XLimit"} and \code{"YLimit"} and computes a conditional survival surface.
}
\usage{
survRestricted(bivarSurf, tauX = Inf, tauY = Inf)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bivarSurf}{
  A matrix containing the bivariate survival probability and the marginal survival probabilities. The first column of \code{"bivarSurf"} is the survival probability on one axis (let's call it X). The first row of \code{"bivarSurf"} is the survival probability on the other axis (let's call it Y). The row names of \code{"bivarSurf"} are ordered \code{"X"}'s. The column names of \code{"bivarSurf"} are ordered \code{"Y"}'s. Element \code{"bivarSurf[1,1]"} equals to 1. Its first row and column names are \code{"0"} and \code{"0"} (see examples below). The function \code{survDabrowska} returns a matrix in this format.
}
  \item{tauX}{
  A value that defines the restricted region for the X coordinate.
}
  \item{tauY}{
  A value that defines the restricted region for the Y coordinate.
}
}
\details{
  If the probability of surviving in the restricted region is equal or less than zero, \code{"NA"} value is returned. Otherwise, the function returns a list of survival surfaces and their differentials. Element \code{"Sxy"} of this list is the conditional survival probability with its marginals and row/column names (in the same format as argument \code{"bivarSurf"}). The rest of the returned list elements do not contain marginal values or row/column names.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{Sxy}{Conditional survival probability. The first column of \code{"Sxy"} is the survival probability on one axis (let's call it X). The first row of \code{"Sxy"} is the survival probability on the other axis (let's call it Y). The row names of \code{"Sxy"} are ordered times to events (or censoring events) for X. The column names of \code{"Sxy"} are ordered times to events (or censoring events) for Y. (see examples below).
  }
  \item{SxMyM}{ \code{"Sxy"} at point point \code{"(x-, y-)"}, where \code{"x-"} is a left limit of \code{"x"}
  }
  \item{Sx}{Conditional marginal survival probability function for X (in a matrix form)}
  \item{Sy}{Conditional marginal survival probability function for Y (in a matrix form)}
  \item{Sdx}{Conditional marginal probability mass function for X (in a matrix form)}
  \item{Sdy}{Conditional marginal probability mass function for Y (in a matrix form)}
  \item{SxM}{Conditional marginal survival probability function for X (in a matrix form) at point x-}
  \item{SyM}{Conditional marginal survival probability function for Y (in a matrix form) at point y-}
  \item{SxM_y}{Conditional bivariate survival probability function at point (x-, y)}
  \item{Sx_yM}{Conditional bivariate survival probability function at point (x, y-)}
  \item{Sdx_y}{\code{"SxM_y - Sxy"}}
  \item{Sx_dy}{\code{"Sx_yM - Sxy"}}
  \item{Sdx_yM}{\code{"SxMyM - Sx_yM"}}
  \item{SxM_dy}{\code{"SxMyM - SxM_y"}}
  \item{Sdxdy}{Conditional bivariate probability mass function}
}
\references{
Eden, S.K., Li, C., Shepherd B.E. (2019). Non-parametric Estimation of Spearman's Rank Correlation with Bivariate Survival Data.
}
\author{
%%  ~~who you are~~
Svetlana Eden <svetlana.eden@vanderbilt.edu>
}

\examples{
X = c(0.5, 0.6, 0.8)
Y = c(0.44, 0.77, 0.99)
deltaX = c(1, 0, 1)
deltaY = c(1, 1, 1)
bivarSurf = survDabrowska(X, Y, deltaX, deltaY)$DabrowskaEst
bivarSurf
  ###             0      0.44      0.77 0.99
  ### 0   1.0000000 0.6666667 0.3333333    0
  ### 0.5 0.6666667 0.6666667 0.3333333    0
  ### 0.6 0.6666667 0.6666667 0.3333333    0
  ### 0.8 0.0000000 0.0000000 0.0000000    0

condSurf = survRestricted(bivarSurf, tauX = Inf, tauY = 0.88)$Sxy
condSurf
  ### > condSurf
  ###       0 0.44 0.77
  ### 0   1.0  0.5    0
  ### 0.5 0.5  0.5    0
  ### 0.6 0.5  0.5    0
  ### 0.8 0.0  0.0    0
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~conditional }
\keyword{ ~bivariate }
\keyword{ ~survival }
