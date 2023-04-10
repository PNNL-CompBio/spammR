#' calcSigScore: calculates a signature score on a particular SingleCellExperiment object
#' This is a numeric score based on the mean rank of the proteins in teh signature
#' @param sce SingleCellExperiment Object
#' @param sigProts list of gene names
#' @param sigNamme name of signature
#' @import SingleCellExperiment
calcSigScore<-function(sce, sigProts,sigName){

  library(SingleCellExperiment)
  allranks<-apply(exprs(sce),2,function(x) rev(rank(x)))

  allpercs<-apply(exprs(sce),2,function(x) percent_rank(x))
  rownames(allpercs)<-rownames(exprs(sce))
  siggenes<-allpercs[intersect(rownames(allpercs),sigProts),]
  sigScore=apply(siggenes,2,mean,na.rm=TRUE)

  colData(sce)[[sigName]]<-sigScore
  sce
}
