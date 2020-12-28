#' Mann Kendall Trend test, corrected for AR(1)
#'
#' Computes the Mann-Kendall trend test, S, and converts it into a Z score.
#'
#' @param x time series of a single grid cell
#' @return Mann-Kendall's S statistic transformed to a Z score
#' @export mk_z_stat

mk_z_stat<- function(x){
  if(any(is.na(x))) x<- x[!is.na(x)]
  if(length(x)<8) return(NA)

  s<- mk_s(x)
  n<- length(x)
  var_s<- n*(n-1)*(2*n+5)/18
  if (length(unique(x)) < n) {
    tmp <- unique(n)
    for (i in 1:length(tmp)) {
      tie <- length(which(x == tmp[i]))
      if (tie > 1) {
        var_s = var_s-tie*(tie-1)*(2*tie+5)/18
      }
    }
  }
  z<- NA
  if(s>0) z<- (s-1)/sqrt(var_s)
  if(s==0) z<- 0
  if(s<0) z<- (s+1)/sqrt(var_s)
  return(z)
}
