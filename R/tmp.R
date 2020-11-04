##  NASA GISTEMP EXAMPLE
# setwd("/Users/co27vij/Desktop/PhD/multiple_testing/")
# source("lib/source_all.R")
library(tidyr)
library(ncdf4)
library(lubridate)

# change to own file location
gistemp <- ncdf4::nc_open("/home/jose/multiple_testing/data/gistemp1200_ERSSTv5.nc")

# grab temperature variable, name dimensions appropiately
temp_gistemp<- ncvar_get(gistemp, "tempanomaly")
dimnames(temp_gistemp) <- list(lon = gistemp$dim$lon$vals,
                               lat = gistemp$dim$lat$vals,
                               t = as.character(as.Date(gistemp$dim$time$vals, origin = "1800-01-01")))
# rm(gistemp)

# Filter out data to include only >= 1951
temp_gistemp <- temp_gistemp[,,year(dimnames(temp_gistemp)$t) > 1950]

# aggregate data
make_yearly <- function(x, groups) {
  x<- split(x, groups)
  unlist(unname(lapply(x, mean, na.rm = TRUE)))
}

temp_gistemp<- apply(temp_gistemp, c("lon", "lat"), make_yearly, groups = year(dimnames(temp_gistemp)$t))
temp_gistemp<- aperm(temp_gistemp, c(2,3,1)) # reorder to lon, lat, t
dimnames(temp_gistemp) <- list(lon = gistemp$dim$lon$vals,
                               lat = gistemp$dim$lat$vals,
                               t = 1951:2018) # rename...
saveRDS(temp_gistemp, file = "data/temp_anomalies_1951_2018.rds")
# mk trend test
# pw_fn<- function(x){
#   if((sum(is.na(x))/length(x))>0) return (NA)
#   if (any(is.finite(x) == FALSE)) {
#     x <- x[-c(which(is.finite(x) == FALSE))]
#     warning("The input vector contains non-finite numbers. An attempt was made to remove them")
#   }
#   xn <- (x[-1] - (x[-length(x)] * rk_fn(x)))
#   z<- mk_z_stat_(xn)
#   return(z)
# }
# mk<- apply(temp_gistemp, 1:2, pw_fn)
# x11()
# image.plot(mk)
#
# mk_perm<- perm_dist(temp_gistemp, fx = pw_fn, nperm = 10) # can be slow
#
# robinson <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +over"
# mk_raster<- mk %>% t %>% .[90:1,] %>% raster(crs = CRS("+init=epsg:4326"), xmn=-179, xmx=179, ymn=-89, ymx=89)
# mk_raster<- projectRaster(mk_raster, crs = CRS(robinson), method = "bilinear")
#
# mk_map<- tm_shape(mk_raster) +
#   tm_raster("layer", legend.show = TRUE, palette = "-RdBu", n = 8, title = "",
#             legend.hist = TRUE, legend.is.portrait = FALSE,
#             colorNA = "#e3e3e3") +
#   tm_legend(legend.outside=T, legend.outside.position="bottom") +
#   map_layer
# tmap_arrange(mk_map)
#
#
# sig<- apply(mk, 1:2, function(x, thr) ifelse(abs(x) <= qnorm(1-ALPHA/2), 1, # not significant
#                        ifelse(x > qnorm(1-ALPHA/2) & x <= thr, 2, # significant and positive without correction
#                               ifelse(x < qnorm(ALPHA/2) & x >= -thr, 3, # significant and negative without correction
#                                      ifelse(x > qnorm(1-ALPHA/2) & x > thr, 4, # significant and positive with correction
#                                             ifelse(x < qnorm(ALPHA/2) & x < -thr , 5, NA))))),
#             thr = 2.192464 )
# # 2.192464 for k-FWER , k = 500 (~5% of tentative discoveries)
# # 2.631384 for k = 150
# # 4.026173 for FWER from maxT
# sig_raster<- sig %>% t %>% .[90:1,] %>% raster(crs = CRS("+init=epsg:4326"), xmn=-179, xmx=179, ymn=-89, ymx=89)
# sig_raster<- projectRaster(sig_raster, crs = CRS(robinson), method = "ngb")
# sig_map<- tm_shape(sig_raster) +
#   tm_raster("layer", palette = my_palette,
#             legend.show = TRUE,
#             style = 'fixed', breaks = 0:5+.5,
#             colorNA = "#e3e3e3")  +
#   tm_legend(legend.outside=T, legend.outside.position="bottom") +
#   map_layer
# tmap_arrange(sig_map)
