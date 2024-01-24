#' spatialDiffEx: does differential expression using annotations in object
#' @import limma
#' @export
#' @param spe Spatial Experiment object
#' @param log_transformed Is the data given in spe log2 transformed ("Y") or not ("N")
#' @param category_col Name of the column that specifies category of each sample. Example: "IsletStatus"
#' #Categories from category_col will be compared in the differential expression analysis
#' @param compare_vals A vector containing names of categories from category_col to be compared. example: c('Proximal','Distal')
#' If length(compare_vals) = 1, i.e. only one category is specified (example: 'Proximal'), then that category will be compared against all others. Example: Proximal vs. Not proximal
#' @returns Spatial Experiment object containing results from differential expression analysis, in addition to what was already present in the input spe
spatialDiffEx<-function(spe,log_transformed,category_col, compare_vals){
  library(limma)

  #collect samples by factor
  samp1<-which(colData(spe)[[category_col]]==compare_vals[1])
  if(length(compare_vals)>1){
    samp2<-which(colData(spe)[[category_col]]==compare_vals[2])
  }else{
    samp2=setdiff(1:nrow(colData(spe)),samp1)
  }

  fac <- factor(rep(c(2,1),c(length(samp2), length(samp1))))
  #print(fac)
  design <- model.matrix(~fac)
  #print(design)
  dat = matrix()
  if (log_transformed=="Y"){
    dat = assay(spe)
  }else if(log_transformed=="N"){
    dat = log2(assay(spe))
  }
  fit <- lmFit(dat[,c(samp2,samp1)], design)
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
