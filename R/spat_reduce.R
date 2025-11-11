#' Reduce a high resolution image measurement (origin) to a lower resultion
#' measurement (target)
#' @description `spat_reduce()` Creates a single SpatialExperiment object
#' with all measurements mapped to the same coordinate space for plotting
#' and analysis
#' @details add detauls here 
#' @import ggplot2
#' @import ggnewscale
#' @import SpatialExperiment
#' @import dplyr

spat_reduce <- function(spe_target,
                        spe_origin,
                        target_dims = c('spot_width','spot_height'),
                        origin_dims = c('spot_width','spot_height')
                        ){
  
  #first let's get overlapping regions
  target_coords <- spatialCoords(spe_target)
  origin_coords <- spatialCoords(spe_origin)
  target_x <- target_coords[,1]
  origin_x <- origin_coords[,1]
  
  ##get overlapping x values
  
  target_y <- target_coords[,2]
  
  #get overlapping y-values
  #then add them to the target
  
  #return target
}