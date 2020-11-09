#' Main analysis function
#'
#' @param data Array containing data. First two dimensions are assumed to be
#' the spatial dimensions, third dimension must be the variable. Note function fails
#' if data is not as specified.
#' @param fx function to be applied at each grid cell. Should be self
#' sufficient - no extra arguments needed besides the time series.
#' Should return only the test statistic
#' @param method a string of characters indicating which correction for multiple
#' testing to use. Defaults to c("maxT","stcs"). Additional options are
#' "bonferroni", "bh", "by", "holmes", "hochberg". "bh" is the
#' Benjamini-Hochberg method, also known as the false discovery rate, "by" is
#' the Benjamini-Yekutieli method. See references for more information.
#' @param alpha_local Significance level for the hypthesis test performed at
#' each grid cell.
#' @param alpha_global significance level to be applied to the permutation
#' methods. This controls the overall probability of a false positive among all
#' the data at the specified alpha. Recommended to be at the same significance
#' level as the thr. Defaults to alpha_local.
#' @param null_distribution either "normal" or "t". Used to estimate the
#' threshold of significance for the test statistic. Defaults to normal
#' @param seed seed to be fed into set.seed function
#' @param block_size Desired block size for block permutation. Useful for
#' serially correlated data.
#' @param verbose Counter returning when the function is done with 10 function
#' calls
#' @return A named list of matrices: the first two contain the test statistic and the
#' p-values for each grid cell, afterwards each matrix contains true/false for
#' every grid cell in the input data, indicating significance /nonsignificance
#' according to the corresponding method
#' @export multiple_testing_correction
#' @references
#' Cortés, J., Mahecha, M., Reichstein, M. et al. Accounting for multiple testing in the analysis of spatio-temporal environmental data. Environ Ecol Stat 27, 293–318 (2020). https://doi.org/10.1007/s10651-020-00446-4

multiple_testing_correction<- function(data,
                                       fx,
                                       method="all",
                                       nperm=1000,
                                       alpha_local=0.05,
                                       alpha_global=NULL,
                                       null_distribution="normal",
                                       seed=NULL,
                                       block_size=NULL,
                                       verbose=TRUE
                                       ){

  if(is.null(alpha_global)) alpha_global<- alpha_local

  perm_results<- perm_dist(data=data, fx=fx, nperm=nperm, alpha_local=alpha_local,
                  alpha_global=alpha_global, null_distribution=null_distribution,
                  seed=seed, block_size=block_size, verbose=verbose)
  out<- threshold_data(perm_results=perm_results, alpha_local=alpha_local,
                       alpha_global=alpha_global, data_dim=dim(data),
                       null_distribution=null_distribution)
  class(out)<- "mtc"
  summary.mtc(out)
  return(out)
}
