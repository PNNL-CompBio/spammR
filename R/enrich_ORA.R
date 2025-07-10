#' Calculate functional or pathway enrichment
#' @description `enrich_ora()` calculates over-representation statistics (ORA) using an interest list of genes from differential 
#' expression results in spammR and gene sets (either the ones provided in spammR or user supplied)
#' This function uses results from `calc_spatial_diff_ex`. Interest list of genes for ORA is obtained from spatialDiffEx results
#' based on the criteria specified in this function.It is a wrapper for the `leapR` package which is required
#' For ORA using an external or already defined interest list of genes and gene sets, use leapR functions directly
#' @export
#' @import leapR
#' @importFrom Biobase fData
#' @import SummarizedExperiment
#' @param spe SpatialExperiment object containing spatial omics data and spatial diffex results
#' @param geneset in GMT format
#' @param feature_column Column of rowData that maps to gene set
#' @param pval_type_forThresh Choose from "adjusted_pval" or "pval". Type of p-value that should be used for
#' filtering statistically significant results. Default is adjusted p-value for multiple hypotheses correction.
#' @param pval_thresh value to use for filtering based on pval_type_forThreshold. Default is 0.05. 
#' Values less than pval_thresh will be kept.
#' @param logFC_lowerThresh Lower threshold for log Fold Change, to be used for filtering spatialDiffEx results. Default is NA
#' @param logFC_upperThresh Upper threshold for log Fold Change, to be used for filtering spatialDiffEx results. Default is NA
#' @param geneset_name Name of geneset provided
#' @param sortResultsBy For sorting ORA results, choose from the following column names: "BH_pvalue" (default)
#' @param comparison_name Example: "RSPv_vs_others" Text to indicate in results data frame, which spatial groups were 
#' compared for the interest list of genes
#' @returns A dataframe containing results from over-representation analysis of members of gene sets in the 
#' interest list of genes based on filtering criteria above.
#' 
#' @examples
#' data(smallPancData)
#' data(pancMeta)
#' data(protMeta)
#' pooledPanc <- dplyr::bind_cols(smallPancData)
#' panc.spe <- convert_to_spe(pooledPanc,pancMeta,protMeta,feature_meta_colname='pancProts')
#' diffex.spe <- calc_spatial_diff_ex(panc.spe,category_col='IsletOrNot')
#' library(leapR)
#' data('krbpaths')
#' ora.res <- enrich_ora(diffex.spe,geneset=krbpaths,geneset_name='krbpaths', feature_column='PrimaryGeneName')
#' 

enrich_ora <- function(spe,
                      geneset,
                      feature_column, #primary gene name to be mapped to enrichment data
                      pval_type_forThresh='adjusted_pval',
                      pval_thresh = 0.05,
                      logFC_lowerThresh=NA,
                      logFC_upperThresh=NA,
                      geneset_name = 'msigdb',
                      sortResultsBy,
                      comparison_name=''){


  # Filter spatialDiffEx results to create interest list of genes for ORA based on user-specified criteria
  # Current assumption is that spatialDiffEx results file has columns for protein names and corresponding gene names. "PG.genes"
  # If gene names are not present in spatialDiffEx results, a helper function (which I will add later), can be run to obtain gene names from the UniProt db files
  # and have a "PG.genes" column added to spatialDiffEx results excel file.

  res <- as(spe,'ExpressionSet') ##convert to expressionSet   

  if (pval_type_forThresh == "adjusted_pval") {
      pval_col_text = "adj.P.Val"
  }else if (pval_type_forThresh == "pval") {
      pval_col_text = "P.Value"
  }else {
      # Throw error "Invalid value for pval_type_forThresh"
  }
  
  fvals <- colnames(Biobase::fData(res))
  
  pval_column = fvals[grep(pval_col_text,fvals)]
  logfc_column = fvals[grep("logFC",fvals)]
  
  # Filter based on p-value criteria
#  int_list = sp_diffEx[!is.na(sp_diffEx[,colnum_pvalue]),]
#  int_list = int_list[int_list[,colnum_pvalue] < pval_thresh,]

  # Filter based on Log fold change criteria
  if (!is.na(logFC_lowerThresh)) {
    low_vals = which(Biobase::fData(res)[,logfc_column] < logFC_lowerThresh)
    if(length(low_vals)==0)
        stop("No values below `logFC_lowerThresh")
    res = res[low_vals,] 
    #int_list = int_list[int_list[,colnum_logfc] > logFC_lowerThresh,]
  }
  if (!is.na(logFC_upperThresh)) {
    high_vals = which(Biobase::fData(res)[,logfc_column] > logFC_upperThresh)
    if(length(high_vals)==0)
        stop("No values above `logFC_upperThresh")
    res = res[high_vals,]      
    #int_list = int_list[int_list[,colnum_logfc] < logFC_upperThresh,]
  }
  
  ora.res <- leapR(eset=res, geneset = geneset,enrichment_method='enrichment_in_sets',
                   id_column=feature_column,
                   primary_column=pval_column,threshold=pval_thresh,
                   greaterthan=FALSE)
  
  return(ora.res)
  # int_list_geneNames = int_list[,feature_column]
  # 
  # es = c()
  # es_sorted = c()
  # # Over-representation analysis: Run Enrichment in Sets (es) in leapR
  # background_genes = SummarizedExperiment::rowData(spe)[,feature_column]
  # es = leapR::leapR(geneset = geneset,
  #           enrichment_method = 'enrichment_in_sets',
  #           targets = int_list_geneNames, background = background_genes)
  # # Sort results by column name specified by user
  # es_sorted =  es[order( es[,"BH_pvalue"]),]
  # num_sets = dim(es_sorted)[1]
  # es_sorted = cbind(rep(comparison_name,num_sets), rownames(es_sorted), es_sorted)
  # colnames(es_sorted)[1:2] = c("Comparison", geneset_name)

  #return(es_sorted)
}
