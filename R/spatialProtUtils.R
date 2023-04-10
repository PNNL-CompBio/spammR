########
#spatialUtils
#
#script containing helper files for the spatialProtExperiment object
#
#


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


#' expToLongForm: helper function that moves expression
#' matrix to long form
expToLongForm<-function(sce,rowname='prot'){
  
  exprs(sce)%>%
    as.matrix()%>%
    as.data.frame()%>%
    tibble::rownames_to_column(rowname)%>%
    tidyr::pivot_longer(c(2:(1+ncol(exprs(spat.phos)))),names_to='Voxel',values_to='LogRatio')
}

#' calcCorrelationWithScore: calculates the correlation of each
#' element with a numeric vector of the score, such as distance or 
#' an immune score, puts value in RowData
#' @param sce
#' @param scoreName: name of score
#' @param protVal: name of feature, e.g. prot or substrate
#' @param method: spearman or pearson
#' @return SingleCellExperiment objecti
calcCorrelationWithScore<-function(sce, 
                                   scoreName,
                                   protVal='prot',
                                   method='spearman'){
    
    #join value with expression
    score.dat<-colData(sce)%>%
      as.data.frame()%>%
      dplyr::select(all_of(scoreName))%>%
      tibble::rownames_to_column('Voxel')%>%
      left_join(expToLongForm(sce,protVal))
    
    #calculate correlation
    prot.cor<-score.dat%>%
      dplyr::rename(pname=protVal,score=scoreName)%>%
      group_by(pname)%>%
      summarize(corVal=cor(score,LogRatio,method=method))
   
    newDataFrame<-rowData(sce)
    newDataFrame[[paste(scoreName,'correlation')]]<-unlist(tibble::column_to_rownames(prot.cor,'pname'))
    rowData(sce)<-newDataFrame
    
    sce
}


#' plotFeatureGrid: plots a numeric value of a single feature or set of features
#' @param sce
#' @param features 
#' 
plotFeatureGrid<-function(sce,feats,featname){
  require(ggplot2)
  require(SingleCellExperiment)
  
  p<-calcSigScore(sce,feats,featname)%>%
    plotSigGrid(featname)

    return(p)
}


#'plotSigGrid: plots a numeric value using the Xcoord and Ycoord columns
#'takes a numeric score from the colDAta
#'@param sce SingleCellExperiment
#'
plotSigGrid<-function(sce,sigName){
  require(ggplot2)
  require(SingleCellExperiment)

  vars<-colData(sce)%>%
    as.data.frame()%>%
    dplyr::rename(signature=sigName)
  
  p<-ggplot(vars,aes(x=Xcoord,y=Ycoord,fill=signature))+
    geom_raster()+scale_fill_viridis_c()+
    theme_bw()+
    ggtitle(paste0(sigName,' signature'))
  
  return(p)
}


#'spatialDiffEx: does differential expression using annotations in object
#'
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



#' buildnetwork
#' Builds network using SingleCellExpression object and values
#' that were applied to the row data
#' @import PCSF
buildNetwork<-function(sce,featName,beta=.5,nrand=100){
   require(dplyr)
  
  ##emtpy weight vectors
  phos.vals<-c()
  prot.vals<-c()
  
  if('Phosphosite'%in%names(rowData(sce))){
    phos.vals<-rowData(sce)[,featName]
    names(phos.vals)<-rowData(sce)[,'Phosphosite']
    phos.vals<-phos.vals[which(phos.vals>0)]
  }

  if('Protein'%in%names(rowData(sce))){
    prot.vals<-rowData(sce)[,featName]
    names(prot.vals)<-rowData(sce)[,'Protein']
    prot.vals<-prot.vals[which(prot.vals>0)]
  }
  
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
  data("STRING")
  
  if(length(phos.vals)>0){
    #read kinase substrate database stored in data folder
    KSDB <- read.csv('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                                  stringsAsFactors = FALSE)
  
    kdat<-KSDB%>%group_by(GENE)%>%
      dplyr::select(SUB_GENE,SUB_MOD_RSD)%>%
      rowwise()%>%
      dplyr::mutate(subval=paste(SUB_GENE,SUB_MOD_RSD,sep='-'))
    
    allvals<-unique(kdat$subval)
    
    mval<-mean(STRING$cost)
    adf<-apply(kdat,1,function(x)
      #for each substrate interaction, add a link from the kinase gene -> substreate -> substrate gene
      data.frame(from=c(x[['GENE']],x[['subval']]),to=c(x[['subval']],x[['SUB_GENE']]),
                 cost=c(mval/2,mval/3)))%>%  ##arbitrary costs based on mean cost of edges around network
      do.call(rbind,.)
  }else{
    adf<-data.frame()
  } 
  
  ##first get diffex proteins
  ppi <- construct_interactome(rbind(STRING,adf))
  
  ##now run the code
  terms=c(phos.vals,prot.vals)
  #print(terms)
  subnet<-NULL
  try(
    subnet <- PCSF_rand(ppi,abs(terms), n=nrand, r=0.2,w = 4, b = beta, mu = 0.0005)
  )
  
  if(is.null(subnet))
    return("")
  
  lfcs<-terms[match(names(V(subnet)),names(terms))]
  lfcs[is.na(lfcs)]<-0.0
  
  ##assign proteins first
  types<-rep('proteins',length(names(V(subnet))))
  
  names(types)<-names(V(subnet))
  
  ##then assign phosphosites
  types[intersect(names(V(subnet)),allvals)]<-'phosphosite'
  
 
  subnet<-igraph::set.vertex.attribute(subnet,'logFoldChange',value=lfcs)
  subnet<-igraph::set.vertex.attribute(subnet,'nodeType',value=types)
  subnet<-igraph::set.edge.attribute(subnet,'interactionType',value='protein-protein interaction')
  
  
  write_graph(subnet,format='gml',file=paste0(featName,'network.gml'))
  return(list(graph=subnet,fname=paste0(featName,'.networkgml')))
}



#' adjustPhosphoWithGlobal
#' For some phospho analysis, we might want to investigate the phosphosites that are changing
#' indpendently of the proteomiocs data. Therefore we consume two independent data objects and subtract one from the other
#' @param phos.obj phosphoproteomic data
#' @param prot.obj proteomic data
#' @return phos.obj
adjustPhophoWithGlobal<-function(phos.obj,prot.obj){
  
  #first match column headeers
  samps <- intersect(colnames(exprs(phos.obj)),colnames(exprs(prot.obj)))
  
  #change new expresion values to subtract protein values when available
  newExpr<-expToLongForm(phos.obj,'Phosphosite')%>%
    tidyr::separate(Phosphosite,into=c('Protein','Site'),sep='_')%>%
    dplyr::rename(oldLogRatio='LogRatio')%>%
    left_join(expToLongForm(prot.obj,'Protein'))%>%
    mutate(diff=(oldLogRatio-LogRatio))%>%
    subset(!is.na(diff))%>%
    tidyr::unite(c(Protein,Site),col='Phosphosite',sep='_')%>%
    dplyr::select(-c(oldLogRatio,LogRatio))%>%
    tidyr::pivot_wider(names_from=Voxel,values_from=diff)%>%
    tibble::column_to_rownames('Phosphosite')%>%
    as.matrix()
    
  #create new singleCellExperiment
  new.phos<-SingleCellExperiment(assays=list(logcounts=as(newExpr,'dgCMatrix')),
                                         colData=colData(phos.obj),
                                         rowData=rownames(newExpr))
  
  new.phos
  
  
}

