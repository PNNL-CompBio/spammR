#' Reduce a high resolution image measurement (origin) to a lower resultion
#' measurement (target)
#' @description `spat_reduce()` Creates a single SpatialExperiment object
#' with all measurements mapped to the same coordinate space for plotting
#' and analysis
#' @details add detauls here 
#' @import SpatialExperiment
#' @import IRanges
#' @param spe_target Low resolution image measurements
#' @param spe_origin High-resolution image measreuements to be mapped to target
#' @param target_dims Names of width and height of target spots, respectively
#' @param stat Statistic used to summarize regions, choices are currently
#' 'mean' and 'median'
#' @return A copy of the spe_target object with the altExperiment slot containing
#' a reduced version of the origin 
#' @example Add example here
#' 
spat_reduce <- function(spe_target,
                        spe_origin,
                        origin_assay = 'metabolomics',
                        target_dims = c('spot_width','spot_height'),
                        stat = 'mean'
                        ){
  
  #first let's get overlapping regions
  #hoping this will be easier when we move to polygon geometry
  target_coords <- spatialCoords(spe_target)
  origin_coords <- spatialCoords(spe_origin)
  
  #get data itself
  odat <- SummarizedExperiment::assay(spe_origin, origin_assay)
  
  #get mean of overlapping regions
  new_dat <- vapply(seq_len(nrow(target_coords)),function(i){
      x_min = target_coords[i,1]
      x_max = x_min + unlist(colData(spe_target)[i,target_dims[1]])
      y_min = target_coords[i,2]
      y_max = y_min + unlist(colData(spe_target)[i,target_dims[2]])
      
      xvals <- intersect(which(origin_coords[,1] < x_max),
                which(origin_coords[,1] > x_min))
      
      yvals <- intersect(which(origin_coords[,2] < y_max),
                          which(origin_coords[,2] > y_min))
      
      ovals <- odat[,intersect(xvals,yvals)]
      return(rowMeans(ovals))

  }, numeric(nrow(odat)))
  
  colnames(new_dat) <- rownames(colData(spe_target))
  
  res_spe <- convert_to_spe(dat = new_dat,
                            feature_meta = rowData(spe_origin),
                            sample_meta = colData(spe_target),
                            assay_name = origin_assay)
  #then add them to the target
  altExp(spe_target, origin_assay) <- res_spe
  
  #return updated target
  return(spe_target)
}