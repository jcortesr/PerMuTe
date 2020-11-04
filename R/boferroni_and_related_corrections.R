#' Bonferroni Hochberg correction
#'
#' For given matrix containing p.values, returns a TRUE/FALSE vector where TRUE
#' is significant using the correction.
#'
#' @param pval matrix of p-values
#' @param alpha desired alpha level for the correction.

bh_fn = function(pval, alpha = .05) {
  ppp = order(pval, decreasing=FALSE)
  m<- length(pval)
  thr = ((1:m)/m)*alpha
  bh = pval[ppp] < thr
  for (i in length(bh):1) {
    if (is.na(bh[i])) next
    if (bh[i]) {
      bh[1:i] = TRUE
      break
    }
  }
  res.bh = rep(NA,m)
  res.bh[ppp] = bh
  out<- matrix(res.bh, nrow = dim(pval)[1])
  return(out)
}

#' Bonferroni Yekutieli correction
#'
#' For given matrix containing p.values, returns a TRUE/FALSE vector where TRUE
#' is significant using the correction.
#'
#' @param pval matrix of p-values
#' @param alpha desired alpha level for the correction.

by_fn = function(pval, alpha = .05) {
  # returns a TRUE/FALSE vector where TRUE is significant using the by correction
  ppp = order(pval, decreasing=FALSE)
  m<- length(pval)
  thr = ((1:m)/(m*sum(1/(1:m))))*alpha
  by = pval[ppp] <= thr
  for (i in length(by):1) {
    if (is.na(by[i])) next
    if (by[i]) {
      by[1:i] = TRUE
      break
    }
  }
  res.by = rep(NA,m)
  res.by[ppp] = by
  out<- matrix(res.by, nrow = dim(pval)[1])
  return(out)
}

#' Walker min-p correction
#'
#' For given matrix containing p.values, returns a TRUE/FALSE vector where TRUE
#' is significant using the correction.
#'
#' @param pval matrix of p-values
#' @param alpha desired alpha level for the correction.


walker_fn<- function(pval, alpha = .05){
  pmin<- min(pval)
  m<- length(pval)
  alpha_w<- 1-(1-alpha)^(1/m)
  res.walker<- pval < alpha_w
  out<- matrix(res.walker, nrow = dim(pval)[1])
  return(out)
}

#' Bonferroni correction
#'
#' For given matrix containing p.values, returns a TRUE/FALSE vector where TRUE
#' is significant using the correction.
#'
#' @param pval matrix of p-values
#' @param alpha desired alpha level for the correction.

bonf_fn<- function(pval, alpha = .05){
  tmp<- pval < (alpha/length(pval))
  out<- matrix(tmp, nrow = dim(pval)[1])
  return(out)
}

#' Holmes step-up correction
#'
#' For given matrix containing p.values, returns a TRUE/FALSE vector where TRUE
#' is significant using the correction.
#'
#' @param pval matrix of p-values
#' @param alpha desired alpha level for the correction.

holm_fn<- function(pval, alpha = .05){
  ppp = order(pval, decreasing=FALSE)
  m<- length(pval)
  thr<- rev(alpha/(m-(m:1)+1))
  out = pval[ppp] < thr
  for (i in length(out):1) {
    if (is.na(out[i])) next
    if (out[i]) {
      out[1:i] = TRUE
      break
    }
  }
  res.out= rep(NA,m)
  res.out[ppp] = out
  out<- matrix(res.out, nrow = dim(pval)[1])
  return(out)
}

#' Hochberg step-down correction
#'
#' For given matrix containing p.values, returns a TRUE/FALSE vector where TRUE
#' is significant using the correction.
#'
#' @param pval matrix of p-values
#' @param alpha desired alpha level for the correction.


hochberg_fn<- function(pval, alpha = .05){
  ppp = order(pval, decreasing=FALSE)
  m<- length(pval)
  thr<- rev(alpha/(m-(m:1)+1))
  out<-  pval[ppp] <= thr
  if(sum(out, na.rm = TRUE) > 0){
    for (i in 1:length(out)) {
      if (!out[i]) {
        out[1:(i-1)] = TRUE
        break
      }
    }
  }
  res.out= rep(NA,m)
  res.out[ppp] = out
  out<- matrix(res.out, nrow = dim(pval)[1])
  return(out)
}
