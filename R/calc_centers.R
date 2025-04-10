
#' Calculate centers of spots/samples for distance based analysis
#' @description calc_centers calculates the centers of spots for numerical analysis
#' @import SpatialExperiment
#' @export
#' @param spe SpatialExperiment object containing omics data
calcCenters<-function(spe){

  ##get relevant column data - x,y coords and spot height/width
  spatial_coords <- spatialCoords(spe)
  spatial_sizes <-colData(spe)[,c('spot_width','spot_height')]

  new_coords<-data.frame(spatial_coords[,1]+spatial_sizes[,1]/2,spatial_coords[,2]+spatial_sizes[,2]/2)
  colnames(new_coords)<-c('x_center','y_center')

}
