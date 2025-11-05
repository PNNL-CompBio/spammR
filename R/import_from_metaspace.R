
#' Create a `SpatialFeatureExperiment` object from metaspace
#' @description `convert_to_sfe()` Puts omics data (omics measurements,
#'  metadata and image(s) corresponding to samples' tissue) into a 
#'  SpatialFeatureExperiment (SFE) object. Most spammR functions require the input 
#'  data to be an SPE or SFE objects.
#' @export
#' @import SummarizedExperiment
#' @import sf
#' @import SpatialExperiment
#' @import SpatialFeatureExperiment
#' @import zellkonverter
#' @import reticulate
#' @param project_id metaspace project id
#' @returns spe.out a SpatialFeatureExperiment (SFE) object that contains all 
#' data and image(s). Ready to be used as input in spammR functions that 
#' require SPE object as input.
import_from_metaspace <- function(project_id = "2024-02-15_20h37m13s", 
                                  fdr_val = 0.2) {
  
  #get anndata object from metaspace
  mc <- reticulate::import('metaspace_converter')

  
  ad = mc$metaspace_to_anndata(dataset_id = project_id, 
                                                fdr = fdr_val,
                                                add_optical_image = TRUE)
  #convert to scexperiment
#  sce <- zellkonverter::AnnData2SCE(ad)
  
  coords <- as.matrix(reducedDim(ad, "spatial"))
  colnames(coords) = c("x","y")
  spe <- SpatialExperiment(
    assay = assay(ad,"X"), 
    colData = ad@colData, 
    spatialCoords = coords,
  )
  spe[["sample_id"]] <- project_id
  addImg(spe,)
  img <- ad$uns$spatial$image$images$hires #4d matrix
  save_mat(img,filename='test.png',dev='png')
  #reformat to sfe
  
  
}