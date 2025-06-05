#' Identify differentially abundant features across image
#' @description `calc_spatial_diff_ex()` Calculates differential expression analysis using annotations in a SpatialExperiment object
#' @import limma
#' @import SpatialExperiment
#' @export
#' @param spe Spatial Experiment object containing data to be used for differential expression analysis
#' @param assay_name Name of the dataset stored in the spe object, that is to be used for the differential expression analysis. Example: znormalized_log2
#' @param log_transformed Is the data given in spe log2 transformed TRUE or FALSE
#' @param category_col Name of the column that specifies category of each sample. Example: "IsletOrNot"
#' #Categories from `category_col` will be compared in the differential expression analysis
#' @param compare_vals A vector containing names of categories from category_col to be compared. Only required if there are more than two values in `category_col`
#' @returns A Spatial Experiment object containing differential expression results, stored in rowData(diffEx.spe)
#  and assays(diffEx.spe) which contains the dataset on which differential expresssion analysis was carried out
#'
#' @examples
#' data(smallPancData)
#' data(pancMeta)
#' data(protMeta)
#' pooledData<-dplyr::bind_cols(smallPancData)
#' pooled.panc.spe <- convert_to_spe(pooledData,
#'                 pancMeta,
#'                 protMeta,
#'                 feature_meta_colname = 'pancProts',
#'                 samples_common_identifier='')
#' diffex.spe <- calc_spatial_diff_ex(pooled.panc.spe,
#'                 category_col='IsletOrNot')
#' 
calc_spatial_diff_ex <- function(spe,
                               assay_name='proteomics',
                               log_transformed=FALSE,
                               category_col,
                               compare_vals){
  #collect samples by factor
  factors <- unique(SummarizedExperiment::colData(spe)[[category_col]])
  if(length(factors)<1){
    ##throw error we need at least two categories
  }else if(length(factors)>2){
    if(missing(compare_vals) || length(setdiff(compare_vals,factors))>0){
      ##throw error - need exactly 2 values in category_col or 2 values in category_col to compare
    }
    factors=compare_vals
  }

  ##now select the samples for each category
  samp1<-which(SummarizedExperiment::colData(spe)[[category_col]]==factors[1])
  samp2<-which(SummarizedExperiment::colData(spe)[[category_col]]==factors[2]) #Later, limma call does samp2 vs. samp1 analysis
  comparison_name = paste(factors[1],"_vs_",factors[2],sep="")

  ##create design matrix with two factors
  fac <- factor(rep(c(2,1),c(length(samp2), length(samp1))))
  #print(fac)
  design <- stats::model.matrix(~fac)
  #print(design)
  dat = SummarizedExperiment::assays(spe)[[ assay_name ]]
#  rownames(dat) = SummarizedExperiment::rowData(spe)[rownames(dat),feature_colname] # Rownames for dat, so that results from limma later will also have corresponding rownames
  if(!log_transformed){
    dat = log2(dat)
  }
  fit <- limma::lmFit(dat[,c(samp2,samp1)], design)
  fit <- limma::eBayes(fit)

  diffex.spe = spe # initialize
  SummarizedExperiment::assays(diffex.spe) = list()
  SummarizedExperiment::assays(diffex.spe,withDimnames=FALSE) [[ assay_name ]] <- dat # Output spe will only have one entry in assays, which will be the dataset used for the differential expression analysis
  # withDimnames=FALSE drops rownames from dat while saving into the assay. We want this because rownames(spe.out) may not necessarily have rownames, depending on how the spe was defined.
  # Differential expression results will be stored in rowData(diffex.spe)
  res <- limma::topTable(fit, coef=2, number=Inf) # Sorting by P-value not needed here because later we are returning all results in the order of the genes in the input SPE, to be consistent.
  colnames_res<-paste(paste(comparison_name, colnames(res), "limma",sep="."))
  #res = data.frame(cbind(rownames(res),res))

  colnames(res) = colnames_res #c(feature_colname,colnames_res)
  # Make sure the results are in the same order of the features in the input SPE object
  diffEx <- cbind(SummarizedExperiment::rowData(spe),res)#|>
    #dplyr::left_join(res)
  SummarizedExperiment::rowData(diffex.spe) <- diffEx

  return(diffex.spe)
}
