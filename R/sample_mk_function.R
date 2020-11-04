#' Sample MK trend test function
#'
#' This function is an example of a valid function as an input to the multiple
#' testing correction function. Autocorrelation lag 1 , r, is estimated and
#' removed from each grid cell time series ,x_t, so thatso that the new time
#' series is y_t = x_t - rx_t-1
#' @param x time series at a single grid cell
#' @return Mann-Kendall's trend test Z score, corrected for autocorrelation
#' @export sample_mk_function

sample_mk_function<- function(x){
  x <- x[is.finite(x)]
  xn <- (x[-1] - (x[-length(x)] * rk_fn(x)))
  z<- mk_z_stat(xn)
  return(z)
}
