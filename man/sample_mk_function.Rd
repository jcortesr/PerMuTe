% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_mk_function.R
\name{sample_mk_function}
\alias{sample_mk_function}
\title{Sample MK trend test function}
\usage{
sample_mk_function(x)
}
\arguments{
\item{x}{time series at a single grid cell}
}
\value{
Mann-Kendall's trend test Z score, corrected for autocorrelation
}
\description{
This function is an example of a valid function as an input to the multiple
testing correction function. Autocorrelation lag 1 , r, is estimated and
removed from each grid cell time series ,x_t, so thatso that the new time
series is y_t = x_t - rx_t-1
}
