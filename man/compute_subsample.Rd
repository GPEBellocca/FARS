% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_subsample.R
\name{compute_subsample}
\alias{compute_subsample}
\title{Compute Subsample of Data by Block}
\usage{
compute_subsample(data, block_ind, n, sample_size = 1, seed = NULL)
}
\arguments{
\item{data}{A numeric matrix or data frame containing the original data.}

\item{block_ind}{Integer vector specifying the end indices of each block. If \code{NULL}, assumes a single block spanning all columns.}

\item{n}{Integer specifying the number of blocks.}

\item{sample_size}{Numeric between 0 and 1 specifying the proportion of columns to sample within each block. Defaults to 1 (all columns).}

\item{seed}{Optional integer. Seed for reproducibility of the column sampling. If \code{NULL}, sampling is random.}
}
\value{
A list with two elements:
\describe{
\item{sample_data}{A numeric matrix combining the sampled columns from all blocks.}
\item{sample_block_ind}{An integer vector containing the indices corresponding to the sampled data.}
}
}
\description{
Selects a random subset of columns within each block of the input data matrix.
}
\keyword{internal}
