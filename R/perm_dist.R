#' Permutation distribution
#'
#' Obtains the permutation distribution of the maximum statistic and the STCS,
#' named maxT and stcs, respectively. Returns a list containing maxT and stcs:
#' each a vector of length nperm where each entry is the maximum of the nth
#' permutation. Together form the permutation distribution, of which
#' the (1-alpha)th percentile is the threshold of significance adjusted
#' for multiple testing.
#'
#' @param data Data in a 3d array, where the first two dimensions are the
#' physical dimensions (e.g., lon and lat), and the third one is time.
#' @param fx function to be applied at each grid cell. Should be self
#' sufficient - no extra arguments needed. Should return only the test
#' statistic
#' @param nperm number of permutations. Defaults to 1000
#' @param verbose Counter returning when the function is done with 10 function
#' calls
#' @return returns the distribution of maxT and stcs, each a vector in a list
#' @export perm_dist

perm_dist<- function(data, fx, nperm=1000,
                     alpha_local, alpha_global, null_distribution,
                     block_size = NULL, seed, verbose = TRUE){
  perm_matrix<- perm_matrix(nobs = dim(data)[3], nperm = nperm, block_size = block_size, seed = seed)
  maxT<- vector(length = nperm)
  stcs<- vector(length = nperm)
  cat("starting permutations:\n")
  for(i in 1:nperm){
    tmp<- apply(data[,,perm_matrix[i,]], 1:2, fx)
    maxT[i]<- max(abs(as.vector(tmp)), na.rm = TRUE)
    stcs[i]<- get_stcs(tmp, alpha_local, null_distribution)$stcs
    if(verbose) if((i%%10)==0) cat(i,"\n")
  }
  cat("finished!\n\n")
  return(list(maxT = maxT, stcs = stcs, original_results = tmp))
}
