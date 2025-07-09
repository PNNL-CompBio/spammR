#' Plot data in heatmap together with image
#' @description `spatial_heatmap()` Creates a spatial heatmap for a given feature in a SpatialExperiment object (spe) and provides the option to label each sample in the x-y plane
#' @details Assumes that every sample (column) in assay(spe) has a corresponding x,y coordinate given in spatialCoords(spe)
#' @details Assumes that colData(spe) and spatialCoords(spe) are provided in the same order of samples (rows).
#' @details Assumes that the spe object contains background image data, stored in imgData(spe) if plotBackground_img is specified to be TRUE below.
#' @details Default values are specified for function parameters 'metric_display', 'label_column' and 'interactive,' if the user doesn't specify those.
#' @import ggplot2
#' @import ggnewscale
#' @import SpatialExperiment
#' @import dplyr
#' @importFrom plotly ggplotly
#' @import scales
#' @import ggpubr
#' @import png
#' @export
#' @param spe SpatialExperiment (SPE) object
#' @param assay_name Name of assay in the spe object that contains data to be plotted
#' @param plotBackground_img Boolean (TRUE or FALSE) to indicate whether a background image should be plotted. Default is FALSE. If TRUE, the parameters image_boundaries and image_sample_id must be specified, and the image data must be present in the spe object, under imgData(spe)
#' @param sample_id The names of the background image's sample_id fields in the spe object; this provides a unique identifier for the background image to be used for plotting (if there are multiple images under imgData(spe)) and only plots samples associated with the specified image sample_id. Example: c("Image0","Raw_noMarkings"). Image data stored under the spe object can be viewed by imgData(spe)
#' @param image_id The name of the background image image_id. Together with the `sample_id` this provides a unique imge to plot.
#' @param image_boundaries Background image's corners'coordinates. These are need to make sure that the background image lines up withe samples' coordinates correctly. Must specify in the following format: c(xmin_image, ymin_image, xmax_image, ymax_image). For example: c(0,0,21,25). These must be in the same coordinate system as the spatial coordinates for the samples in the SPE object (spatialCoords(spe)).
#' @param spatial_coord_type Position type for the given spatial coordinates of samples in spe. Current options are: "topleft_corner", "topright_corner". Default is blank, which assumes coordinates are bottom right.
#' @param spatial_coord_names Names of x and y spatial coordinates respectively in the spe. Example: c("Xcoord","Ycoord") or c("X","Y").
#' @param spot_size Is a vector of length 2 describing the width and height of each spot.
#' @details Future versions of spatialPlot_feature() to include additional options for spatial_coord_type: "center", bottomleft_corner", "bottomright_corner"
#' @param feature_type Example: "GeneName". Default is whatever identifier is used in rownames. The name of feature_type must be present as a column in rowData(spe)
#' @param feature Name of the feature in the spe object whose values are to plotted in the spatial heat map. This should be a row name in rowData(spe)
#' @param metric_display Legend title for spatial heatmap. If this parameter is not specified, legend title defaults to "Protein abundance measure"
#' @param label_column Column in colData(spe) to be used for labeling grid squares. If not specified, default is no labels.
#' @param sample_label_color Color to be used for labels of samples/grid squares. Default is white.
#' @param sample_label_size Font size for sample labels. Default is 1.75.
#' @param plot_title Title to be given to the spatial heatmap. Default is "Spatial signature for XYZ" where XYZ is the name of the specified feature
#' @param interactive Boolean value (TRUE/FALSE) indicating whether the plot should have interactive mouse hovering. If not specified, this defaults to TRUE. Note: grid squares can only be labeled when interactive = FALSE due to current ggplotly limitations.
#' @returns A spatial heatmap of the chosen feature
#'
#' @examples
#' data(pancMeta)
#' data(smallPancData)
#' data(protMeta)
#' img0.spe <- convert_to_spe(smallPancData$Image_0,
#'                         pancMeta,
#'                         protMeta,
#'                         feature_meta_colname='pancProts',
#'                         image_files=system.file("extdata",'Image_0.png',package='spammR'),
#'                         image_samples_common_identifier='Image0',
#'                         spatialCoords_colnames=c('x_pixels','y_pixels'),
#'                         samples_common_identifier = 'Image0',
#'                         image_id='with_grid')
#' res = spatial_heatmap(img0.spe, 
#'         feature='INS', 
#'         sample_id='Image0', 
#'         image_id='with_grid', 
#'         feature_type='PrimaryGeneName',
#'         spatial_coord_names=c('x_pixels','y_pixels'), 
#'         spot_size=unlist(colData(img0.spe)[1,c('spot_width','spot_height')]), 
#'         image_boundaries=unlist(colData(img0.spe)[1,c('x_origin','y_origin','x_max','y_max')]),
#'         label_column='IsletOrNot', interactive=FALSE)
#' 
spatial_heatmap <- function(spe,
                          feature, ##feature to plot!
                          feature_type=NA, #element of rowdata to use
                          assay_name = 'proteomics',
                          plotBackground_img=TRUE,
                          sample_id,
                          image_id,
                          image_boundaries,
                          spatial_coord_type='',
                          spatial_coord_names = c('Xcoord','Ycoord'),
                          spot_size = c(1,1), ##change this to default to 1 pixel
                          metric_display = "Protein abundance measure",
                          label_column=NA,
                          sample_label_color="white",
                          sample_label_size=1.75,
                          plot_title=NULL,
                          interactive=FALSE){
  ##first get the spatial coordinates from the metadata
  xcoord_name = spatial_coord_names[1]
  ycoord_name = spatial_coord_names[2]
  spatial = SummarizedExperiment::colData(spe)[,c(xcoord_name,ycoord_name)]###this isnt getting rownames: as.data.frame(spatialCoords(spe))

  spatial[,xcoord_name] = as.numeric(spatial[,xcoord_name])
  spatial[,ycoord_name] = as.numeric(spatial[,ycoord_name])
  x = spatial[,xcoord_name]
  y = spatial[,ycoord_name]
  rownames(spatial) <- rownames(SummarizedExperiment::colData(spe))

  ##now we can get the feature data
  f = SummarizedExperiment::assays(spe,withDimnames = FALSE)[[ assay_name ]]
  feature_values_toplot = c()

  if (is.na(feature_type)) { # default is whatever is used for rownames of f
    rowNum_toplot  = which(rownames(f) %in% feature)
    } else{
      rowNum_toplot = which(SummarizedExperiment::rowData(spe)[,feature_type] %in% feature)
    }

  if (length(feature) == 1) { ##we are just plotting a single row
      feature_values_toplot  = f[rowNum_toplot,]
  }else { ##we are plotting more than one value and need to average
      feature_values_toplot = colMeans(f[rowNum_toplot,],na.rm=TRUE)
      feature_values_toplot = feature_values_toplot[is.finite(feature_values_toplot)]
  }
  if(length(feature_values_toplot)==0)
      stop(paste("Cannot plot heatmap:",paste(feature,collapse=','),"not in dataset"))
  #create a title of the plot
  spatial_meta = SummarizedExperiment::colData(spe)
  if (is.null(plot_title)) {
    if(length(feature) == 1)
      title = paste("Expression of", feature, "in",sample_id)
    else
      title = paste("Expression of",length(feature), "features in",sample_id)

  }else {
    title = plot_title
  }
  lab = ""
  if (is.na(label_column)) {
    lab = NA
  }else {
    lab = spatial_meta[,label_column]
  }
  # Switched to using geom_rect instead of goem_raster because geom_raster positioning gets distorted when the plot is made interactive.
  # If not doing interactive hovering, then geom_raster may be more favorable to use.
  x_left = c()
  x_right = c()
  y_botton = c()
  y_top = c()
  if (spatial_coord_type == "topright_corner") {
    x_left = x - spot_size[1]
    x_right = x
    y_bottom = y - spot_size[2]
    y_top = y
  }else if (spatial_coord_type == "topleft_corner") {
    x_left = x
    x_right = x + spot_size[1]
    y_bottom = y - spot_size[2]
    y_top = y
  }else{ ##default to bottom left corner
    x_left = x
    x_right = x + spot_size[1]
    y_bottom = y
    y_top = y + spot_size[2]
  }
  midpoint_x = (x_left + x_right)/2
  midpoint_y = (y_bottom + y_top)/2
  # Background image has to match sample and image identifier
  img_sample_id = sample_id
  img_image_id = image_id

  # Row corresponding to the background image of interest, to be used for plotting
  imgData_rowNum = which(SpatialExperiment::imgData(spe)$sample_id == img_sample_id & 
                             SpatialExperiment::imgData(spe)$image_id == img_image_id)
  background_img = SpatialExperiment::imgData(spe)$data[[imgData_rowNum]]
  # Background image boundaries
  xmin_image = image_boundaries[1]
  ymin_image = image_boundaries[2]
  xmax_image = image_boundaries[3]
  ymax_image = image_boundaries[4]
  #img_png = ping::readPNG(background_img)
  # Scale feature_values_toplot to show relative values. Scale from 0 to 1. This scale is useful when comparing spatial plots of different proteins
  rescaled_feature_values = scales::rescale(feature_values_toplot,to = c(0,1))

  spatial = cbind(spatial,feature_values_toplot,rescaled_feature_values,
                  x_left,x_right,y_bottom,y_top,midpoint_x,midpoint_y) |>
    as.data.frame()
 # print(spatial)
  p <- ggplot2::ggplot(spatial, ggplot2::aes(xmin = x_left, xmax = x_right, 
                                             ymin = y_bottom, ymax = y_top, 
                                             alpha = 0.8, 
                                             fill = feature_values_toplot, label = lab)) +
    ggpubr::background_image(background_img) +
    ggplot2::geom_rect() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::geom_label(ggplot2::aes(x = midpoint_x,y = midpoint_y),
                        label.size = NA, 
                        fill = NA, colour = sample_label_color, 
                        size = sample_label_size) +
    ggplot2::labs(fill = metric_display) +
  #  ggnewscale::new_scale_fill()+
  #  ggplot2::geom_rect(ggplot2::aes(xmin = x_left, xmax = x_right, ymin = y_bottom, ymax = y_top, fill=rescaled_feature_values))+
  #  ggplot2:: geom_label(ggplot2::aes(x=midpoint_x,y=midpoint_y),label.size = NA, fill=NA, colour = sample_label_color, size=sample_label_size)+
  #  ggplot2::scale_fill_viridis_c(limits=c(0,1))+
   # ggplot2::labs(fill = "Scaled values (min=0, max=1)")+
    ggplot2::theme_bw() +
    ggplot2::xlim(xmin_image,xmax_image) +
    ggplot2::ylim(ymin_image,ymax_image) +
    ggplot2::coord_fixed(ratio = 1,expand = FALSE) +  # expand=FALSE to make sure the origin for the image is where it should be (without padding), to make sure the image lines up correctly with samples' coordinates
    ggplot2::xlab("x") +
    ggplot2::ylab("y") +
    ggplot2::ggtitle(title)
  if (interactive) {
    p <- plotly::ggplotly(p)
  }
  return(p)
}
