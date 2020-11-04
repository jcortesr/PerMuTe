#' Gridded world temperature data .
#'
#' A dataset containing the temperature at each grid cell in the world.
#'
#'
#' @format An array with dimensions 180, 90, 68
#' \describe{
#'   \item{lon}{longitude, in degrees}
#'   \item{lat}{latitude, in degrees}
#'   \item{t}{time, in years, starting in 1951}
#' }
#' @source \url{https://data.giss.nasa.gov/pub/gistemp/GHCNv3/gistemp1200_ERSSTv5.nc.gz}
"temp_gistemp"
