#' Lag k serial correlation
#'
#' computes lag-k serial correlation coefficient of vector x , r_k
#'
#' @param x time series
#' @return the lag-k serial correlation coefficient
#' @export rk_fn

rk_fn<- function(x,k=1){
  # computes lag-k serial correlation coefficient of vector x , r_k
  n<- length(x)
  mean_x<- mean(x)
  num<- (1/(n-k))*sum((x[1:n-k] - mean_x)*(x[(1+k):n]-mean_x))
  denom<- (1/n)*sum((x[1:n] - mean_x)^2)
  rk<- num/denom
  return(rk)
}
