% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_stcs.R
\name{get_stcs}
\alias{get_stcs}
\title{Supra-Threshold Cluster Size}
\usage{
get_stcs(data, alpha_local, null_distribution, data_dim)
}
\arguments{
\item{data}{Matrix containing data. Should be a statistical image, i.e., the
test statistic of each grid cell}

\item{data_dim}{dimesions of original data. Used for calculating df for
t distributed test statistics, ignored if the test statistic is normal}

\item{thr}{Threshold of significance.}
}
\value{
A number specifying the size, in number of grid cells, of the largest
exceedance cluster.
}
\description{
From data in a matrix, finds the largest cluster size of grid cells above
the specified threshold.
}
