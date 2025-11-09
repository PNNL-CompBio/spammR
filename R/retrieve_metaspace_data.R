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
#' @export
#' @examples
#' mspe <- retrieve_metaspace_data()
#' 
#' name
retrieve_metaspace_data <- function(project_id = "2024-02-15_20h37m13s", 
                                    fdr = 0.2, 
                                    assay_name = 'metabolites', 
                                    sample_id = 'sample', 
                                    rotate = FALSE){
    #load function
    ms <- reticulate::import_from_path(path = '../inst/python/',
                                       module = 'convert_metaspace_entry')    
  
    #convert database to tuple
    db <- reticulate::tuple('SwissLipids','2018-02-02')
    
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
    #convert to spatial object
    mspe <- spammR::convert_to_spe(dat = dat, 
                           sample_meta = sdata,
                           feature_meta = fdata,
                           spatial_coords_colnames = c('x_pixels','y_pixels'),
                           assay_name = assay_name, 
                           sample_id = sample_id)
    return(mspe)
    
    
}