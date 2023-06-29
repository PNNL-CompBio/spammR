# spatialPlot_feature: Creates a spatial heatmap for a given feature and labels each sample in the x-y plane
# Assumes that every sample (column) in assay(spe) has a corresponding x,y coordinate given in spatialCoords(spe)
# Assumes that colData(spe) and spatialCoords(spe) are provided in the same order of samples (rows).
# Default values are specified for metric_display, label_column and interactive, if the user doesn't specify those.
#' @param spe: SpatialExperiment object
#' @param feature: Name of the feature in the spe object whose values are to plotted in the spatial heat map. This should be a row name in rowData(spe)
#' @param metric_display: Legend title for spatial heatmap. If this parameter is not specified, legend title defaults to "Protein abundance measure"
#' @param label_column: Colunm in colData(spe) to be used for labeling grid squares. If not specified, default is no labels.
#' @param interactive: Boolean value (TRUE/FALSE) indicating whether the plot should have interactive mouse hovering. If not specified, this defaults to TRUE.
#  Note: grid squares can only be labeled when interactive = FALSE
spatialPlot_feature<-function(spe,feature,metric_display = "Protein abundance measure",label_column=NA,interactive=TRUE){
  spatial<-spatialCoords(spe)%>%
    as.data.frame()
  spatial$x = as.numeric(spatial$x)
  spatial$y = as.numeric(spatial$y)
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
  p<- ggplot(spatial, aes(xmin = x-1, xmax = x, ymin = y-1, ymax = y, fill=feature_values_toplot, label = lab))+
    geom_rect()+
    scale_fill_viridis_c()+
    geom_label(aes(x=(x+x-1)/2,y=(y+y-1)/2),label.size = NA, fill=NA)+
    labs(fill = metric)+
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
    ggplotly(p)
  }else{
    p
  }
}
