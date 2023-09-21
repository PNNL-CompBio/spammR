#' Creates a spatial heatmap for a given feature in a SpatialExperiment object (spe)
#' and provides the option to label each sample in the x-y plane
#' @export
# Assumes that every sample (column) in assay(spe) has a corresponding x,y coordinate given in spatialCoords(spe)
# Assumes that colData(spe) and spatialCoords(spe) are provided in the same order of samples (rows).
# Default values are specified for function paramters 'metric_display', 'label_column' and 'interactive,' if the user doesn't specify those.
#' @param spe: SpatialExperiment object
#' @param spatial_coord_type Position type for the given spatial coordinates of samples in spe. Current options are: "topleft_corner", "topright_corner"
#' Future versions of spatialPlot_feature() to include additional options for spatial_coord_type: "center", bottomleft_corner", "bottomright_corner"
#' @param feature: Name of the feature in the spe object whose values are to plotted in the spatial heat map. This should be a row name in rowData(spe)
#' @param metric_display: Legend title for spatial heatmap. If this parameter is not specified, legend title defaults to "Protein abundance measure"
#' @param label_column: Colunm in colData(spe) to be used for labeling grid squares. If not specified, default is no labels.
#' @param sample_label_color Color to be used for labels of samples/grid squares. Default is white.
#' @param interactive: Boolean value (TRUE/FALSE) indicating whether the plot should have interactive mouse hovering. If not specified, this defaults to TRUE.
#' @return spatial_plot: Spatial heatmap of the chosen feature
#  Note: grid squares can only be labeled when interactive = FALSE due to current ggplotly limitations.
spatialPlot_feature<-function(spe,spatial_coord_type,feature,metric_display = "Protein abundance measure",label_column=NA,sample_label_color="white",interactive=TRUE){
  library(ggplot2)
  library(SpatialExperiment)
  library(dplyr)
  library(plotly)
  spatial<-spatialCoords(spe)%>%
    as.data.frame()
  spatial$x = as.numeric(spatial$x)
  spatial$y = as.numeric(spatial$y)
  x = spatial$x
  y = spatial$y
  feature_values_toplot = assay(spe)[feature,rownames(spatial)]
  spatial_meta = colData(spe)
  title = paste("Spatial signature for ", feature, sep="")
  lab = ""
  if (is.na(label_column)){
    lab = NA
  }else{
    lab = spatial_meta[,label_column]
  }
  # Switched to using geom_rect instead of goem_raster because geom_raster positioning gets distorted when the plot is made interactive.
  # If not doing interactive hovering, then geom_raster may be more favorable to use.
  x_left = c()
  x_right = c()
  y_botton = c()
  y_top = c()
  if (spatial_coord_type == "topright_corner"){
    x_left = x-1
    x_right = x
    y_bottom = y-1
    y_top = y
  }else if (spatial_coord_type == "topleft_corner"){
    x_left = x
    x_right = x + 1
    y_bottom = y-1
    y_top = y
  }
  midpoint_x = (x_left + x_right)/2
  midpoint_y = (y_bottom + y_top)/2
  p<- ggplot(spatial, aes(xmin = x_left, xmax = x_right, ymin = y_bottom, ymax = y_top, fill=feature_values_toplot, label = lab))+
    geom_rect()+
    scale_fill_viridis_c()+
    geom_label(aes(x=midpoint_x,y=midpoint_y),label.size = NA, fill=NA, colour = sample_label_color)+
    labs(fill = metric_display)+
    theme_bw()+
    xlim(0,NA)+
    ylim(0,NA)+
    xlab("x")+
    ylab("y")+
    ggtitle(title)
  # geom_raster version
  # p<-ggplot(spatial,aes(x=x,y=y, fill=feature_values_toplot, label = lab))+
  #   geom_raster(hjust = 0, vjust = 0) +
  #   scale_fill_viridis_c()+
  #   geom_label(label.size = NA, fill=NA)+
  #   labs(fill = metric)+
  #   theme_bw()+
  #   xlim(0,NA)+
  #   ylim(0,NA)+
  #   ggtitle(title)
  if (interactive){
    p <- ggplotly(p)
  }
  return(p)
}
