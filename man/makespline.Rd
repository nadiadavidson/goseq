\name{makespline}
\Rdversion{1.1}
\alias{makespline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monotonic Spline
}
\description{
Fits a monotonic cubic spline to the data provided, using the penalized constrained least squares method from the \code{mgcv} package.
}
\usage{
makespline(x, y, newX=NULL, nKnots = 6, lower_bound = 10^-3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The predictor variable.
}
  \item{y}{
The response variable.  Must be the same length as \code{x}.
}
 \item{newX}{
The points at which to return the value on the fitted spline.  If not specified \code{x} is used.
}
 \item{nKnots}{
The number of knots to use in fitting the spline.
}
 \item{lower_bound}{
The spline cannot drop below this value.
}
}
\details{
This uses the \code{pcls} function from the \pkg{mgcv} package to produce the fit.  The montonicity constraint is enforced using \code{mono.con} from the same package. The \code{lower_bound} argument is only used on the rare occasions when the fitting function becomes negative or arbitrarily close to zero.  If this does occur \code{lower_bound} is added everywhere to ensure that no one length is given essentially infinite weighting.
}
\value{
Returns a vector of values containing the value of the fit at each point \code{newX}.
}
\references{
Package \pkg{mgcv}.  In particular this function is a modification of an example given in the man page for \code{pcls}.
}
\author{
Matthew D. Young \email{myoung@wehi.edu.au}.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

y <- c( rbinom(50,p=0.4,size=1), rbinom(50,p=0.6,size=1) )
x <- 1:100
plot(x,y)
p <- makespline(x,y)
lines(x,p)
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
