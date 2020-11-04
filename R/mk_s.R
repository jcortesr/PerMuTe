#' Mann Kendall S statistic
#'
#' Calculates the MK S statistic. Internal function used in mk_z_stat.
#'
#' @param x time series
#' @return Mann-Kendall's S statistic
#' @export mk_s
#'
mk_s<- function(x){
  n<- length(x)
  S0_fn<- function (i,x) -sign(x[1:(n-i)] - x[(i+1):n])
  sum(unlist(lapply(1:(n-1),S0_fn,x)))
}
