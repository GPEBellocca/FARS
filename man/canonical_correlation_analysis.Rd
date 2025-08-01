% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/canonical_correlation_analysis.R
\name{canonical_correlation_analysis}
\alias{canonical_correlation_analysis}
\title{Canonical Correlation Analysis for MLDFM}
\usage{
canonical_correlation_analysis(y, Nregio, r_glob, r_reg)
}
\arguments{
\item{y}{A numeric matrix containing all variables (T × N).}

\item{Nregio}{Integer vector specifying the number of variables in each block.}

\item{r_glob}{Integer. Number of global factors to be extracted.}

\item{r_reg}{Integer vector. Number of regional (local) factors to be extracted for each block.}
}
\value{
A numeric matrix of estimated factors (T × \code{r_glob}).
}
\description{
Performs Canonical Correlation Analysis (CCA)
}
\keyword{internal}
