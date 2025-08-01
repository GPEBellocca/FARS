% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_stressed_factors.R
\name{compute_stressed_factors}
\alias{compute_stressed_factors}
\title{Compute Stressed Factors}
\usage{
compute_stressed_factors(dep_variable, factors, ellipsoids, h, qtau, min)
}
\arguments{
\item{dep_variable}{A numeric vector representing the dependent variable (e.g., GDP growth, inflation).}

\item{factors}{A matrix or data frame of factor estimates, typically extracted from an MLDFM model.}

\item{ellipsoids}{A list of matrices, where each matrix represents a stressed ellipsoid for a given time period.}

\item{h}{Integer representing the forecast horizon (in time steps). It defines the lag used in regression.}

\item{qtau}{Numeric. The quantile level used in quantile regression (default is 0.05).}

\item{min}{Logical. If \code{TRUE}, the function uses a stepwise minimization method. If \code{FALSE}, it uses a stepwise maximization method.}
}
\value{
A matrix of stressed factors, with each row representing a time period and each column representing a factor.
}
\description{
Computes stressed factors.
}
\keyword{internal}
