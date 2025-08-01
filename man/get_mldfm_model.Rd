% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mldfm_model.R
\name{get_mldfm_model}
\alias{get_mldfm_model}
\title{Get a Single \code{mldfm} Object from a \code{mldfm_subsample} Object}
\usage{
get_mldfm_model(object, index, ...)
}
\arguments{
\item{object}{An object of class \code{mldfm_subsample}.}

\item{index}{Integer. The position of the desired model (between 1 and \code{n_samples}).}

\item{...}{Additional arguments (ignored).}
}
\value{
A single \code{mldfm} object.
}
\description{
Returns the \code{mldfm} object at the specified position in a \code{mldfm_subsample} object.
}
