#' Print function for mtc class
#'
#' Print the summary of the test statistic, p values, and detected grid cells
#' for each mutliple testing correction.
#'
#' @param object an object of class mtr

print.mtc<- function(object){
  out1<- t(sapply(object[1:2], function(x) summary(as.vector(x))))
  # test_statistic<- summary(as.vector(object$test_statistic))
  # p_values<- summary(as.vector(object$pval))
  out<- as.data.frame(do.call(rbind, lapply(object[-c(1:2)], function(x) table(as.vector(x)))))
  out$method<- rownames(out)
  rownames(out)<- NULL
  names(out)[1:2]<- c("not significant", "significant")
  out<- out[order(out$significant,decreasing = TRUE),]
  out<- out[, c(3,2,1)]
  cat("Summary of test statistic and p-values:\n\n")
  print(round(out1, digits = 4))
  cat("\n\nNumber of significant grid cells by multiple testing correction:\n\n")
  print(out, row.names = FALSE)
}
