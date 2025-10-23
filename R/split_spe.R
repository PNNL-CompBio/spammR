#' Divide SPE object into smaller objects by column features
#'
#' @description `split_spe()` splits an SPE object into smaller objects 
#' by a column feature. Ideally this should be in the SpatialExperiment
#' class or some utility package.
#' @export
#' @param spe SpatialExperiment object containing spatial omics data 
#' and spatial diffex results
#' @param split_colname Column of rowData that maps to gene set
#' @param assay_name Optional name of assay to use in split data
#'
#' @returns A list of SpatialExperiment objects containing a subset of the data
#'
#' @examples
#' data(smallPancData)
#' data(pancMeta)
#' data(protMeta)
#' pooledPanc <- dplyr::bind_cols(smallPancData)
#' panc.spe <- convert_to_spe(pooledPanc, pancMeta, protMeta, 
#'           feature_meta_colname = "pancProts")
#' split_list <- split_spe(panc.spe, split_colname = "Image")
#'
split_spe <- function(spe, split_colname, assay_name = NULL) {
    if (missing(split_colname) | !split_colname %in% names(colData(spe))) {
        stop("Need a column to split on that is in the spe object")
    }
  
    vals <- SummarizedExperiment::colData(spe)[, split_colname, 
                                               drop = TRUE] |> unique()
  
    spe.list <- vapply(vals, function(v) {
        dvals <- which(colData(spe)[, split_colname] == v)
        mdat <- colData(spe)[dvals, ]
        if (!is.null(assay_name)) {
            dat <- assay(spe, assay_name)[, dvals]
        } else {
            dat <- assay(spe)[, dvals]
        }
    
        res <- convert_to_spe(
          dat = dat,
          sample_meta = mdat,
          feature_meta = rowData(spe),
          sample_id = as.character(v)
        )
        imgData(res) <- imgData(spe) ## if the image is there, keep it the same!
        return(list(res))
    }, list(1))
    names(spe.list) <- as.character(vals)
    return(spe.list)
}
