#' @title US GDP Growth Series 
#'
#' @description Quarterly US GDP growth series used as the dependent variable in the replication.
#'
#' @details
#' Derived from the Excel file \emph{Data_IMF.xlsx} included in \code{inst/extdata/}.
#' The original series contains quarterly GDP levels for 63 countries.
#' For replication, all series are converted to log-differenced annualized growth rates
#' (\code{diff(log(x)) * 400}). From this dataset, the U.S. series is extracted
#' and the first observation dropped to obtain 59 observations in total.
#'
#' @format A time series object with 59 quarterly observations.
#' @source Replication materials of González-Rivera et al. (2024).
#' @docType data
#' @name dep_variable
#' @usage data(dep_variable)
#' @keywords datasets
"dep_variable"


#' @title Macro-Financial Database 
#'
#' @description Macro-financial variables used in the replication exercise.
#'
#' @details
#' Derived from the Excel file \emph{DataBase.xlsx} included in \code{inst/extdata/}.
#' The original dataset contains 624 variables. For replication, the first 519 variables
#' are selected, converted to a numeric matrix, and outliers are corrected using
#' the function \code{correct_outliers(..., threshold = 5)} provided in FARS.
#'
#' @format A numeric matrix with 59 rows and 519 columns.
#' @source Replication materials of González-Rivera et al. (2024).
#' @docType data
#' @name mf_data
#' @usage data(mf_data)
#' @keywords datasets
"mf_data"
