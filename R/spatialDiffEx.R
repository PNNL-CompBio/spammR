#' spatialDiffEx: does differential expression using annotations in object
#' BiocManager::install("limma")
#' @export
#' @param spe Spatial Experiment object
#' @param category_col Name of the column that specifies category of each sample. Example: "IsletStatus"
#' #Categories from category_col will be compared in the differential expression analysis
#' @param compare_vals A vector containing names of categories from category_col to be compared. example: c('Proximal','Distal')
#' @returns Spatial Experiment object containing results from differential expression analysis, in addition to what was already present in the input spe
spatialDiffEx<-function(spe,category_col, compare_vals){
  library(limma)

  #collect samples by factor
  samp1<-which(colData(spe)[[category_col]]==compare_vals[1])
  if(length(compare_vals)>1){
    samp2<-which(colData(spe)[[category_col]]==compare_vals[2])
  }else{
    samp2=setdiff(1:ncol(colData(spe)),samp1)
  }

  fac <- factor(rep(c(2,1),c(length(samp2), length(samp1))))
  #print(fac)
  design <- model.matrix(~fac)
  #print(design)
  fit <- lmFit(assay(spe)[,c(samp2,samp1)], design)
  fit <- eBayes(fit)

  # print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="P") #Do we want to sort by "P" or "P.value" or "adj.P.val"
  # res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  colnames(res)<-paste(paste(category_col,'limma'),colnames(res))
  res<-res%>%
    tibble::rownames_to_column('X')
  rd<-rowData(spe)%>%
    as.data.frame()%>%
    full_join(res)
  rowData(spe)<-rd
  return(spe)
}
