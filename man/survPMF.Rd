\name{survPMF}
\alias{survPMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes bivariate probability mass function from bivariate survival probability
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function computes bivariate probability mass function from bivariate survival probability.
}
\usage{
survPMF(bivarSurf)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bivarSurf}{
  A matrix containing the bivariate survival probability and the marginal survival probabilities. The first column of \code{"bivarSurf"} is the marginal survival probability on one axis (let's call it X). The first row of \code{"bivarSurf"} is the marginal survival probability on the other axis (let's call it Y). The row names of \code{"bivarSurf"} are ordered \code{"X"}'s. The column names of \code{"bivarSurf"} are ordered \code{"Y"}'s. Element \code{"bivarSurf[1,1]"} equals 1. The rest of the matrix is the joint survival surface.  Its first row and column names are \code{"0"} and \code{"0"} (see examples below). The function \code{survDabrowska} returns a matrix in this format.
}
}
\details{
  The function returns a list of survival surfaces and their differentials. Element \code{"Sdxdy"} of this list is the bivariate probability mass function with its marginals and row/column names (in the same format as argument \code{"bivarSurf"}). The rest of the returned list elements do not contain marginal values or row/column names.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{Sdxdy}{Bivariate probability mass function with its marginals and row/column names (in the same format as argument \code{"bivarSurf"})}
  \item{Sxy}{Bivariate survival probability.
  }
  \item{SxMyM}{ \code{"Sxy"} at point point \code{"(x-, y-)"}, where \code{"x-"} is a left limit of \code{"x"}
  }
  \item{Sx}{Marginal survival probability function for X (in a matrix form)}
  \item{Sy}{Marginal survival probability function for Y (in a matrix form)}
  \item{Sdx}{Marginal probability mass function for X (in a matrix form)}
  \item{Sdy}{Marginal probability mass function for Y (in a matrix form)}
  \item{SxM}{Marginal survival probability function for X (in a matrix form) at point x-}
  \item{SyM}{Marginal survival probability function for Y (in a matrix form) at point y-}
  \item{SxM_y}{Bivariate bivariate survival probability function at point (x-, y)}
  \item{Sx_yM}{Bivariate survival probability function at point (x, y-)}
  \item{Sdx_y}{\code{"SxM_y - Sxy"}}
  \item{Sx_dy}{\code{"Sx_yM - Sxy"}}
  \item{Sdx_yM}{\code{"SxMyM - Sx_yM"}}
  \item{SxM_dy}{\code{"SxMyM - SxM_y"}}
}
\references{
Eden, S., Li, C., Shepherd B. (2019). Non-parametric Estimation of Spearman's Rank Correlation with Bivariate Survival Data.
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

bivarPMF = survPMF(bivarSurf)$Sdxdy
bivarPMF
  ###             0      0.44      0.77      0.99
  ### 0   0.0000000 0.3333333 0.3333333 0.3333333
  ### 0.5 0.3333333 0.3333333 0.0000000 0.0000000
  ### 0.6 0.0000000 0.0000000 0.0000000 0.0000000
  ### 0.8 0.6666667 0.0000000 0.3333333 0.3333333
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~bivariate }
\keyword{ ~survival }
