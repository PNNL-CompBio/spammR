#' spatial_heatmap: Creates a spatial heatmap for a given feature in a SpatialExperiment object (spe) and provides the option to label each sample in the x-y plane
#' @details Assumes that every sample (column) in assay(spe) has a corresponding x,y coordinate given in spatialCoords(spe)
#' @details Assumes that colData(spe) and spatialCoords(spe) are provided in the same order of samples (rows).
#' @details Assumes that the spe object contains background image data, stored in imgData(spe) if plotBackground_img is specified to be TRUE below.
#' @details Default values are specified for function parameters 'metric_display', 'label_column' and 'interactive,' if the user doesn't specify those.
#' @import ggplot2
#' @import ggnewscale
#' @import SpatialExperiment
#' @import dplyr
#' @import plotly
#' @import scales
#' @import ggpubr
#' @import png
#' @export
#' @param spe SpatialExperiment (SPE) object
#' @param assay_name Name of assay in the spe object that contains data to be plotted
#' @param plotBackground_img Boolean (TRUE or FALSE) to indicate whether a background image should be plotted. Default is FALSE. If TRUE, the parameters image_boundaries and image_sample_id must be specified, and the image data must be present in the spe object, under imgData(spe)
#' @param image_sample_ids c(sample_id, image_id) The names of the background image's sample_id and image_id fields in the spe object; this provides a unique identifier for the background image to be used for plotting (if there are multiple images under imgData(spe)) and only plots samples associated with the specified image sample_id. Example: c("Image0","Raw_noMarkings"). Image data stored under the spe object can be viewed by imgData(spe)
#' @param image_boundaries Background image's corners'coordinates. These are need to make sure that the background image lines up withe samples' coordinates correctly. Must specify in the following format: c(xmin_image, ymin_image, xmax_image, ymax_image). For example: c(0,0,21,25). These must be in the same coordinate system as the spatial coordinates for the samples in the SPE object (spatialCoords(spe)).
#' @param spatial_coord_type Position type for the given spatial coordinates of samples in spe. Current options are: "topleft_corner", "topright_corner"
#' @param spatial_coord_names Names of x and y spatial coordinates respectively in the spe. Example: c("Xcoord,Ycoord) or c(X,Y).
#' @details Future versions of spatialPlot_feature() to include additional options for spatial_coord_type: "center", bottomleft_corner", "bottomright_corner"
#' @param feature_type Example: "GeneName". Default is whatever identifier is used in rownames. The name of feature_type must be present as a column in rowData(spe)
#' @param feature Name of the feature in the spe object whose values are to plotted in the spatial heat map. This should be a row name in rowData(spe)
#' @param metric_display Legend title for spatial heatmap. If this parameter is not specified, legend title defaults to "Protein abundance measure"
#' @param label_column Column in colData(spe) to be used for labeling grid squares. If not specified, default is no labels.
#' @param sample_label_color Color to be used for labels of samples/grid squares. Default is white.
#' @param sample_label_size Font size for sample labels. Default is 1.75.
#' @param plot_title Title to be given to the spatial heatmap. Default is "Spatial signature for XYZ" where XYZ is the name of the specified feature
#' @param interactive Boolean value (TRUE/FALSE) indicating whether the plot should have interactive mouse hovering. If not specified, this defaults to TRUE. Note: grid squares can only be labeled when interactive = FALSE due to current ggplotly limitations.
#' @returns spatial_plot: Spatial heatmap of the chosen feature
#' @examples
#' data(pancMeta)
#' data(pancDataList)
#'
spatial_heatmap<-function(spe,
                          assay_name,
                          plotBackground_img=TRUE,
                          image_sample_ids,
                          image_boundaries,
                          spatial_coord_type,
                          spatial_coord_names,
                          feature_type=NA,
                          feature,
                          metric_display = "Protein abundance measure",
                          label_column=NA,
                          sample_label_color="white",
                          sample_label_size=1.75,
                          plot_title=NULL,
                          interactive=TRUE){
 # library(ggplot2)
  #library(ggnewscale)
  #library(SpatialExperiment)
  #library(dplyr)
  #library(plotly)
  #library(ggpubr)
  #library(png)
  #library(scales)
  spatial = as.data.frame(spatialCoords(spe))
  xcoord_name = spatial_coord_names[1]
  ycoord_name = spatial_coord_names[2]
  spatial[,xcoord_name] = as.numeric(spatial[,xcoord_name])
  spatial[,ycoord_name] = as.numeric(spatial[,ycoord_name])
  x = spatial[,xcoord_name]
  y = spatial[,ycoord_name]
  f = SummarizedExperiment::assays(spe,withDimnames=FALSE)[[ assay_name ]]
  feature_values_toplot = c()
  if (is.na(feature_type)){ # default is whatever is used for rownames of f
    feature_values_toplot = as.numeric(f[feature,rownames(spatial)])
  }else{
    rowNum_toplot = grep(feature,SummarizedExperiment::rowData(spe)[,feature_type])
    feature_values_toplot = as.numeric(f[rowNum_toplot,rownames(spatial)])
  }
  spatial_meta = SummarizedExperiment::colData(spe)
  if (is.null(plot_title)){
    title = paste("Spatial signature for ", feature, sep="")
  }else{
    title = plot_title
  }
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
  # Background image
  img_sample_id = image_sample_ids[1]
  img_image_id = image_sample_ids[2]
  # Row corresponding to the background image of interest, to be used for plotting
  imgData_rowNum = which(SpatialExperiment::imgData(spe)$sample_id==img_sample_id & SpatialExperiment::imgData(spe)$image_id==img_image_id)
  background_img = SpatialExperiment::imgData(spe)$data[[imgData_rowNum]]
  # Background image boundaries
  xmin_image = image_boundaries[1]
  ymin_image = image_boundaries[2]
  xmax_image = image_boundaries[3]
  ymax_image = image_boundaries[4]
  #img_png = ping::readPNG(background_img)
  # Scale feature_values_toplot to show relative values. Scale from 0 to 1. This scale is useful when comparing spatial plots of different proteins
  rescaled_feature_values = scales::rescale(feature_values_toplot,to=c(0,1))
  #spatial = cbind(spatial,rescaled_feature_values)
  p<- ggplot2::ggplot(spatial, aes(xmin = x_left, xmax = x_right, ymin = y_bottom, ymax = y_top, fill=feature_values_toplot, label = lab))+
    background_image(background_img)+
    ggplot2::geom_rect()+
    ggplot2::scale_fill_viridis_c()+
    ggplot2::geom_label(aes(x=midpoint_x,y=midpoint_y),label.size = NA, fill=NA, colour = sample_label_color, size=sample_label_size)+
    ggplot2:: labs(fill = metric_display)+
    ggnewscale::new_scale_fill()+
    ggplot2::geom_rect(aes(xmin = x_left, xmax = x_right, ymin = y_bottom, ymax = y_top, fill=rescaled_feature_values))+
    ggplot2:: geom_label(aes(x=midpoint_x,y=midpoint_y),label.size = NA, fill=NA, colour = sample_label_color, size=sample_label_size)+
    ggplot2::scale_fill_viridis_c(limits=c(0,1))+
    ggplot2::labs(fill = "Scaled values (min=0, max=1)")+
    ggplot2::theme_bw()+
    ggplot2::xlim(xmin_image,xmax_image)+
    ggplot2::ylim(ymin_image,ymax_image)+
    ggplot2::coord_fixed(ratio=1,expand=FALSE)+  # expand=FALSE to make sure the origin for the image is where it should be (without padding), to make sure the image lines up correctly with samples' coordinates
    ggplot2::xlab("x")+
    ggplot2::ylab("y")+
    ggplot2::ggtitle(title)
  if (interactive){
    p <- plotly::ggplotly(p)
  }
  return(p)
}
