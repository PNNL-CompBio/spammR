#' Calculate functional or pathway enrichment from a gradient
#' @description `enrich_gradient()` calculates over-representation statistics (ORA) using a ranking of genes from `distance_based_analysis()`
#' For ORA using an external or already defined interest list of genes and gene sets, use leapR functions directly
#' @export
#' @import leapR
#' @import SummarizedExperiment
#' @param spe SpatialExperiment object containing spatial omics data and spatial diffex results
#' @param geneset in GMT format
#' @param feature_column Column of rowData that maps to gene set
#' @param ranking_column Column of rowData that maps to ranks
#' @returns a dataframe containing results from over-representation analysis of members of gene sets in the interest list of genes based on filtering criteria above.
#' @examples
#'
#' data(pancDataList)
#' data(pancMeta)
#' data(protMeta)
#' img0.spe<-convert_to_spe(pancDataList$Image_0,pancMeta,protMeta,feature_meta_colname='pancProts',image_files=system.file("extdata",'Image_0.png',package='spammR'),image_samples_common_identifier='Image0',samples_common_identifier = 'Image0',image_ids='Image0')
#' img0.spe<-distance_based_analysis(img0.spe,'proteomics',sampleCategoryCol='IsletOrNot',sampleCategoryValue='Islet')
#' library(leapR)
#' data('krbpaths')
#' rank.res <- enrich_gradient(img0.spe, geneset=krbpaths,feature_column='PrimaryGeneName',ranking_column='IsletDistancespearmanCor')
#'

enrich_gradient <-function(spe,
                      geneset,
                      feature_column, #primary gene name to be mapped to enrichment data
                      ranking_column){


  rvals <- SummarizedExperiment::rowData(spe)|>
    as.data.frame()|>
    dplyr::rename(feature=feature_column,rank=ranking_column)|>
    dplyr::select(feature,rank)

  es <- leapR::leapR(rvals, geneset=geneset,
            enrichment_method='enrichment_in_order',
            id_column='feature',primary_columns='rank')|>
    subset(!is.na(pvalue))|>
    dplyr::arrange(pvalue)

  return(es)
}
