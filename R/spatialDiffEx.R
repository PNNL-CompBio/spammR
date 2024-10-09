#' spatialDiffEx: does differential expression analysis using annotations in a SpatialExperiment object
#' @import limma
#' @export
#' @param spe Spatial Experiment object containing data to be used for differential expression analysis
#' @param assay_name Name of the dataset stored in the spe object, that is to be used for the differential expression analysis. Example: znormalized_log2
#' @param log_transformed Is the data given in spe log2 transformed ("Y") or not ("N")
#' @param category_col Name of the column that specifies category of each sample. Example: "IsletStatus"
#' #Categories from category_col will be compared in the differential expression analysis
#' @param compare_vals A vector containing names of categories from category_col to be compared. example: c('Proximal','Distal')
#' If length(compare_vals) = 1, i.e. only one category is specified (example: 'Proximal'), then that category will be compared against all others. Example: Proximal vs. Not proximal
#' @param feature_colname Name of column in rowData(spe), that is to be used for identifying features. Example: "PG.ProteinNames"
#' @returns a list with diffEx_df, a data frame containing the differential expression results and spe_out: Spatial Experiment object containing diffEx, stored in rowData(spe_out)
#  and assays(spe_out) which contains the dataset on which differential expresssion analysis was carried out

spatialDiffEx<-function(spe,assay_name,log_transformed,category_col, compare_vals){
  library(limma)

  #collect samples by factor
  samp2<-which(colData(spe)[[category_col]]==compare_vals[1]) #Later, limma call does samp2 vs. samp1 analysis
  comparison_name = c()
  if(length(compare_vals)>1){
    samp1<-which(colData(spe)[[category_col]]==compare_vals[2])
    comparison_name = paste(compare_vals[1],"_vs_",compare_vals[2],sep="")
  }else{
    samp1=setdiff(1:nrow(colData(spe)),samp2)
    comparison_name = paste(compare_vals[1],"_vs_others",sep="")
  }

  fac <- factor(rep(c(2,1),c(length(samp2), length(samp1))))
  #print(fac)
  design <- model.matrix(~fac)
  #print(design)
  dat = assays(spe)[[ assay_name ]]
  if(log_transformed=="N"){
    dat = log2(dat)
  }
  fit <- lmFit(dat[,c(samp2,samp1)], design)
  fit <- eBayes(fit)

  spe_out = spe # initialize
  assays(spe_out) = list()
  assays(spe_out) [[ assay_name ]] <- dat # Output spe will only have one entry in assays, which will be the dataset used for the differential expression analysis
  # Differential expression results will be stored in rowData(spe_out)
  res <- topTable(fit, coef=2, number=Inf) # Sorting by P-value not needed here because later we are returning all results in the order of the genes in the input SPE, to be consistent.
  colnames(res)<-paste(paste(comparison_name, colnames(res), "limma",sep="."))
  res<-res%>%
    tibble::rownames_to_column(feature_colname) # Turn rownames into a column whose name is given by feature_colname, so it can be matched with the column in the spe
  diffEx<-rowData(spe)%>%
    as.data.frame()%>%
    full_join(res)
  rowData(spe_out) <- as.matrix(diffEx)
  #return(spe)
  return(list("diffEx_df" = diffEx, "spe_out"= spe_out))
}
