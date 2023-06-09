% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simple.linear.regression.R
\name{simple.linear.regression}
\alias{simple.linear.regression}
\title{Simple Linear Regression}
\usage{
simple.linear.regression(y = cars$dist, x = cars$speed, data = cars)
}
\arguments{
\item{y}{vector of length n, with dependent variables. Default is cars data set dist.}

\item{x}{vector of length n, with independent variables. Default is cars data set speed.}

\item{data}{n*p matrix with predixtor variables. Default is cars data set.}
}
\value{
List with the following elements:
        height: coefficient
        intercept: coefficient
        height.se: unit increase in X, expect that Y increases by height.se
        intercept.se: unit increase in Y, expect that X increases by intercept.se
        height.p: significance of change
        intercept.p: significance of change
}
\description{
Simple Linear Regression
}
\examples{
simple.linear.regression(x = c(1,2,3,4,5,6), y = c(1,2,3,4,5,6))
}
\author{
Matthijs Knigge

This is a test function for personal purpose only, reflecting 
template that should be used. It sole purpose is just a reminder.
Simple linear regression: Y = B0 + B1 * X1 + E
}
\keyword{linear}
\keyword{regression}
