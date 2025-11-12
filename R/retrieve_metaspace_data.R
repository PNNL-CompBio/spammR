#' Retreive data objects from metspace
#' @description `retrieve_metaspace_data()`Collects data from 
#' http://metaspace2020.org by project id and returns it in 
#' SpatialExperiment object
#' @export
#' @import reticulate
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @param project_id Identifier for project
#' @param fdr FDR value above which to collect data
#' @param assay_name Name of assay to include in object
#' @param sample_id Name of sample
#' @param rotate Set to TRUE if x and y coordinates need to be swapped
#' @param y_offset Number of pixels to adjust y based on image
#' @param x_offset Number of pixels to adjust x based on image
#' @export
#' @examples
#' mspe <- retrieve_metaspace_data()
#' 
#' name
#' 
#' 

.onLoad <- function(libname, pkgname){
  library(reticulate)
  if (!py_available('pandas'))
    py_install('pandas')
  if (!py_available('numpy'))
    py_install('numpy')
  if (!py_available('metaspace'))
    py_install('metaspace')
  py_require(python_version = ">=3.11")
  
}

retrieve_metaspace_data <- function(project_id = "2024-02-15_20h37m13s", 
                                    fdr = 0.2, 
                                    assay_name = 'metabolites', 
                                    sample_id = 'sample', 
                                    rotate = FALSE,
                                    y_offset = 0,
                                    x_offset=0){

    #load function
    path <- system.file('spammR','inst/python')
    ms <- reticulate::import_from_path(path = path,
                                      module = 'convert_metaspace_entry')
    #convert database to tuple
    db <- reticulate::tuple('SwissLipids','2018-02-02')
    
    message("Downloading ion data from metaspace...")
    ##will download two matrices into one
    datas <- ms$download_ion_data(project_id = project_id,
                                  fdr_val = fdr,
                                  database = db)   
    
    #get feature data to data frame
    fdata <- datas[[2]] |>
      tibble::column_to_rownames('ion')
    
    #get get ion data to matrix
    idata <- datas[[1]]
    dat <- idata[,c('ion','sample_id','intensity')] |>
      tidyr::pivot_wider(names_from = 'sample_id',values_from = 'intensity') |>
      tibble::column_to_rownames('ion')
      
    #get sample information to data frame
    sdata <- idata[,c('sample_id','x_coord','y_coord')] |>
      dplyr::distinct() |>
      tibble::column_to_rownames('sample_id') |>
      dplyr::mutate(spot_height = 1) |>
      dplyr::mutate(spot_width = 1)

    if (rotate) {
      sdata <- sdata |>
        dplyr::rename(x_pixels = 'y_coord', y_pixels = 'x_coord')
    }else{
      sdata <- sdata |>
        dplyr::rename(x_pixels = 'x_coord', y_pixels = 'y_coord')
    }    
    
    #can't always get the images to align...
    sdata$x_pixels <- unlist(sdata$x_pixels) + x_offset
    sdata$y_pixels <- unlist(sdata$y_pixels) + y_offset
    
    #convert to spatial object
    mspe <- spammR::convert_to_spe(dat = dat, 
                           sample_meta = sdata,
                           feature_meta = fdata,
                           spatial_coords_colnames = c('x_pixels','y_pixels'),
                           assay_name = assay_name, 
                           sample_id = sample_id)
    return(mspe)
    
    
}