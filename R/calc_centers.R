#' Calculate centers of spots/samples for distance based analysis
#' @description calc_centers calculates the centers of spots for numerical analysis. Currently not used.
#' @import SpatialExperiment
#'
#' @param spe SpatialExperiment object containing omics data
#' @return Data frame of new coordinates
#'
#' @examples
#' data(pancMeta)
#' data(protMeta)
#' data(smallPancData)
#' img0.spe <- convert_to_spe(smallPancData$Image_0,
#'   pancMeta,
#'   protMeta,
#'   feature_meta_colname = "pancProts",
#'   image_files = system.file("extdata", "Image_0.png", package = "spammR"),
#'   sample_id = "Image0",
#'   spatial_coords_colnames = c("x_pixels", "y_pixels"),
#'   image_sample_ids = "Image0", image_ids = "Image0"
#' )
#' ## can't actually call internal functuon
calc_centers <- function(spe) {
  ## get relevant column data - x,y coords and spot height/width
  spatial_coords <- spatialCoords(spe)

  if (ncol(spatial_coords) == 0) {
    stop("No spatial coordinates were set in image object")
  }
  spatial_sizes <- colData(spe)[, c("spot_width", "spot_height")]

  new_coords <- data.frame(spatial_coords[, 1] + spatial_sizes[, 1] / 2, spatial_coords[, 2] + spatial_sizes[, 2] / 2)
  colnames(new_coords) <- c("x_center", "y_center")
  return(new_coords)
}
