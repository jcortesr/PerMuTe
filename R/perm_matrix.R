#' Permutation Matrix
#'
#' Generates matrix for permutaitons. Each row is a permutation. Makes it easy
#' to do block permutations. Repeated lines are avoided. Original ordering
#' is always the last row.
#' @param nobs number of observations to permute
#' @param nperm number of permutations, including original data ordering.
#' @param block_size Desired block size for block permutation.
#' @param seed seed to be fed into set.seed function
#' @return matrix for use in permutation methods. Each row is one permutation

perm_matrix<- function(nobs, nperm, block_size = NULL, seed = NULL){
  if(!is.null(seed)) set.seed(seed)

  out<- matrix(NA, nrow = nperm, ncol = nobs)
  if (!is.null(block_size)){
    shifter <- function(x, a) {
      if (a == 0) x else c(tail(x, -a), head(x, a))
    }
    for(i in 1:(nperm-1)){
      time_series<- shifter(1:nobs, sample(0:(nobs-1),1))
      blocks<- split(time_series, ceiling(seq_along(time_series)/block_size))
      if(length(blocks[[length(blocks)]]) < block_size) {
        blocks[[(length(blocks)-1)]]<-  c(blocks[[(length(blocks)-1)]],  blocks[[length(blocks)]])
        blocks[[length(blocks)]]<- NULL
      }
      out[i,]<- unname(unlist(blocks[sample.int(length(blocks))]))
    }
    out[nperm,]<- 1:nobs
    while(sum(duplicated(out))>0){
      sel<- which(duplicated(out, fromLast = TRUE))
      for (i in 1:length(sel)) {
        time_series<- shifter(1:nobs, sample(0:(nobs-1),1))
        blocks<- split(time_series, ceiling(seq_along(time_series)/block_size))
        if(length(blocks[[length(blocks)]]) < block_size) {
          blocks[[(length(blocks)-1)]]<-  c(blocks[[(length(blocks)-1)]],  blocks[[length(blocks)]])
          blocks[[length(blocks)]]<- NULL
        }
        tmp<- blocks # lapply(blocks, sample)
        out[sel[i],]<- unname(unlist(tmp[sample.int(length(blocks))]))
      }
    }
    return(out)
  } else{

    for(i in 1:(nperm-1)){
      out[i,]<- sample.int(nobs)
    }
    out[nperm,]<- 1:nobs

    while(sum(duplicated(out))>0){
      sel<- which(duplicated(out, fromLast = TRUE))
      for (i in 1:length(sel)) out[sel[i],]<- sample.int(nobs)
    }
    return(out)
  }
}
