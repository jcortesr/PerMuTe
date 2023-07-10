#' Threshold  data
#' Thresholds the data according to the corresponding multiple testing
#' correction. This function is used internally.
#'
#' @param perm_results object created from perm_dist function. Contains the
#' distribution of maxT and STCS AND the matrix obtained from applying a
#' hypothesis test at each grid cell.
#' @param alpha_local significance threshold at the grid cell level.
#' @param alpha_global desired overall significance level
#' @param data_dim dimesions of original data. Used for calculating df for
#' t distributed test statistics, ignored if the test statistic is normal
#' @export threshold_data


threshold_data<- function(perm_results, alpha_local, alpha_global, data_dim,
                          null_distribution){

  if(null_distribution == "normal") {
    pval<- 2*pnorm(-abs(perm_results$original_results))
    uncorrected_thr<- qnorm(1-alpha_local/2)
  }

  if(null_distribution == "t"){
    pval<- 2*pt(-abs(perm_results$original_results), df = data_dim[3]-2)
    uncorrected_thr<- qt(1-alpha_local/2, df = data_dim[3]-2)
  }

  out<- list()
  out$test_statistic<- perm_results$original_results
  out$pval<- pval
  out$uncorrected<- abs(perm_results$original_results) > uncorrected_thr
  out$bonferroni<- bonf_fn(pval, alpha = alpha_global)
  out$walker<- walker_fn(pval, alpha = alpha_global)
  out$hochberg<- hochberg_fn(pval, alpha = alpha_global)
  out$holm<- holm_fn(pval, alpha = alpha_global)
  out$by<- by_fn(pval, alpha = alpha_global)
  out$bh<- bh_fn(pval, alpha = alpha_global)
  out$maxT<- abs(perm_results$original_results) > quantile(perm_results$maxT, probs = 1-alpha_global, na.rm = TRUE)

  tmp<- get_stcs(data=perm_results$original_results, alpha_local=alpha_local, null_distribution=null_distribution)
  stcs_thr<- quantile(perm_results$stcs, 1-alpha_global, na.rm = TRUE)
  wh_cluster<- which(tmp$clusters$cluster.count > stcs_thr)
  wh_cluster_sel<- tmp$clusters$clusters %in% wh_cluster
  out$stcs<- apply(tmp$clusters$clusters, 1:2, function(x, wh_cluster) wh_cluster %in% x, wh_cluster)

  # using both stcs and maxT in stcs
  stcs_thr<- quantile(perm_results$stcs, 1-alpha_global/2, na.rm = TRUE)
  stcs_maxT_thr<- quantile(perm_results$stcs_maxT, 1-alpha_global/2, na.rm = TRUE)
  wh_cluster_stcs<- which(tmp$clusters$cluster.count > stcs_thr)
  stcs_maxT<- vector(length=length(tmp$clusters$cluster.count))
  for (i in 1:length(tmp$clusters$cluster.count)) stcs_maxT[i]<- max(abs(perm_results$original_results[tmp$clusters$clusters==i]), na.rm = TRUE)
  wh_cluster_stcs_maxT<- which(stcs_maxT > stcs_maxT_thr)
  wh_cluster_sel<- tmp$clusters$clusters %in% c(wh_cluster_stcs, wh_cluster_stcs_maxT)
  out$stcs_maxT<- apply(tmp$clusters$clusters, 1:2, function(x, wh_cluster)  x %in% wh_cluster, wh_cluster=unique(c(wh_cluster_stcs, wh_cluster_stcs_maxT)))

  data_dimnames<- attr(out$test_statistic, "dimnames")
  if(!is.null(data_dimnames)) {
    for(i in 1:length(out)) attr(out[[i]], "dimnames")<- data_dimnames
  }
  return(out)
}

