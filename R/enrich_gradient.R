#' Calculate functional or pathway enrichment from a gradient
#' @description `enrich_gradient()` calculates over-representation statistics 
#' (ORA) using a ranking of genes from `distance_based_analysis()`
#' For ORA using an external or already defined interest list of genes and 
#' gene sets, use leapR functions directly
#' @export
#' @import leapR
#' @import SummarizedExperiment
#' @param spe SpatialExperiment object containing spatial omics data and 
#' spatial diffex results
#' @param geneset in GMT format
#' @param assay_name name of assay
#' @param feature_column Column of rowData that maps to gene set
#' @param ranking_column Column of rowData that maps to ranks
#' @param method Rank enrichmemnt method used by leapR. Default is `ztest`
#' @returns A dataframe containing results from over-representation analysis 
#' of members of gene sets in the interest list of genes based on filtering 
#' criteria above.
#'
#' @examples
#' data(smallPancData)
#' data(pancMeta)
#' data(protMeta)
#' img0.spe <- convert_to_spe(smallPancData$Image_0,
#'   pancMeta,
#'   protMeta,
#'   feature_meta_colname = "pancProts",
#'   image_files = system.file("extdata", "Image_0.png", package = "spammR"),
#'   image_sample_ids = "Image0",
#'   spatial_coords_colnames = c("x_pixels", "y_pixels"),
#'   sample_id = "Image0",
#'   image_ids = "Image0"
#' )
#' img0.spe <- distance_based_analysis(img0.spe,
#'   "proteomics",
#'   sampleCategoryCol = "IsletOrNot",
#'   sampleCategoryValue = "Islet"
#' )
#' library(leapR)
#' data("krbpaths")
#' rank.res <- enrich_gradient(img0.spe,
#'   geneset = krbpaths,
#'   feature_column = "PrimaryGeneName",
#'   ranking_column = "IsletDistancespearmanCor"
#' )
#' head(rank.res)
enrich_gradient <- function(spe,
                            assay_name,
                            geneset,
                            feature_column, # name to be mapped to gene set
                            ranking_column, 
                            method = 'ztest', 
                            ...) { #other methods to pass to leapR
  

    # here we use the leapR package for the rank enrichment
    es <- leapR::leapR(
      eset = spe, assay_name = assay_name, geneset = geneset,
      enrichment_method = "enrichment_in_order",
      method = method,
      id_column = feature_column, primary_columns = ranking_column, ...
    ) |>
      subset(!is.na(pvalue)) |>
      dplyr::arrange(pvalue)
  
    return(es)
}
