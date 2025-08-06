#' Divide SPE object into smaller objects by column features
#' 
#' @description `split_spe()` splits an SPE object into smaller objects by a column feature. Ideally this shoudl be in the SpatialExperiment
#' class
#' @export
#' @param spe SpatialExperiment object containing spatial omics data and spatial diffex results
#' @param split_colname Column of rowData that maps to gene set
#' @param assay_name Optional name of assay to use in split data
#' @returns A list of SpatialExperiment objects containing a subset of the data
#' 
#' @examples
#' data(smallPancData)
#' data(pancMeta)
#' data(protMeta)
#' pooledPanc <- dplyr::bind_cols(smallPancData)
#' panc.spe <- convert_to_spe(pooledPanc,pancMeta,protMeta,feature_meta_colname='pancProts')
#' split_list <- split_spe(panc.spe,split_colname='Image')
#' 
#' 
split_spe <- function(spe, split_colname,assay_name){
    
    if(missing(split_colname) | !split_colname%in%colData(spe))
        stop("Need a column to split on that is in the spe object")
    
    vals <- SummarizedExperiment::colData(spe)[,split_colname] |> unique() 
    
    spe.list <- lapply(vals, function(v){
        dvals <- which(colData(spe)[,split_colname]==v)
        mdat <- colData(spe[dvals,])
        if(!missing(assay_name))
            dat <- assay(spe,assay_name)[,mdat]
        else
            dat <- assay(spe)[,dat]
    
       res = convert_to_spe(dat = dat,
                    sample_meta = mdat,
                    feature_meta = rowData(spe),  
                     sample_id=as.character(v),
                    feature_meta_colname = 'pancProts')
       imgData(res) <- imgData(spe) ## if the image is there, keep it the same!
    })
    return(spe.list)
}