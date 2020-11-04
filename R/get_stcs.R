#' Supra-Threshold Cluster Size
#'
#' From data in a matrix, finds the largest cluster size of grid cells above
#' the specified threshold.
#' @param data Matrix containing data. Should be a statistical image, i.e., the
#' test statistic of each grid cell
#' @param thr Threshold of significance.
#' @return A number specifying the size, in number of grid cells, of the largest
#' exceedance cluster.
#' @param data_dim dimesions of original data. Used for calculating df for
#' t distributed test statistics, ignored if the test statistic is normal
#' @export get_stcs
#' @import osc

get_stcs<- function(data, alpha_local, null_distribution, data_dim){
  if(null_distribution == "normal") thr<- qnorm(1-alpha_local/2)
  if(null_distribution == "t") thr<- qt(1-alpha_local/2, df = data_dim[3]-2)

  pixel_sign<- sign(data)
  pixel_significant<- abs(data)>thr
  pixel_result<- pixel_sign*pixel_significant

  # positive
  pixel_result_pos<- pixel_result
  pixel_result_pos[is.na(pixel_result_pos)]<- -999
  clusters_pos<- osc::cca(pixel_result_pos,count.cells = TRUE, s=1, mode = 2,
                          count.max  = length(pixel_sign))
  stcs_pos<- max(clusters_pos$cluster.count)

  nclust_pos<- length(clusters_pos$cluster.count)
  clusters_sep<- vector(mode = "list", length = 2)

  # negative
  pixel_result_neg<- pixel_result
  pixel_result_neg[pixel_result_neg == -1] = 10
  pixel_result_neg[pixel_result_neg == 1] = -10
  pixel_result_neg[is.na(pixel_result_neg)] = -999
  clusters_neg<- osc::cca(pixel_result_neg,count.cells = TRUE, s=1, mode = 2,
                          count.max = length(pixel_sign))
  stcs_neg<- max(clusters_neg$cluster.count)

  # join
  clusters_neg$clusters[clusters_neg$clusters > 0]<- clusters_neg$clusters[clusters_neg$clusters > 0] + nclust_pos
  clusters_sep[[1]]<- clusters_pos$clusters + clusters_neg$clusters
  clusters_sep[[2]]<- c(clusters_pos$cluster.count, clusters_neg$cluster.count)
  names(clusters_sep)<- c("clusters", "cluster.count")

  stcs<- max(clusters_sep$cluster.count)

  return(list(stcs=stcs, clusters=clusters_sep))
}
