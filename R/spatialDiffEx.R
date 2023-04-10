#'spatialDiffEx: does differential expression using annotations in object
#' @export
## here we do differential expression again
spatialDiffEx<-function(sce,column='pulpAnnotation', vals=c('red','white')){
  library(limma)

  #collect samples by factor
  samp1<-which(colData(sce)[[column]]==vals[1])
  if(length(vals)>1){
    samp2<-which(colData(sce)[[column]]==vals[2])
  }else{
    samp2=setdiff(1:ncol(colData(sce)),samp1)
  }

  fac <- factor(rep(c(2,1),c(length(samp2), length(samp1))))
  #print(fac)
  design <- model.matrix(~fac)
  #print(design)
  fit <- lmFit(exprs(sce)[,c(samp2,samp1)], design)
  fit <- eBayes(fit)

  # print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="P")
  # res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  colnames(res)<-paste(paste(column,'limma'),colnames(res))
  res<-res%>%
    tibble::rownames_to_column('X')
  rd<-rowData(sce)%>%
    as.data.frame()%>%
    full_join(res)
  rowData(sce)<-rd
  return(sce)
}
