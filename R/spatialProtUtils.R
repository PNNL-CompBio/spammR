########
#spatialUtils
#
#script containing helper files for the spatialProtExperiment object
#
#



#' calcCorrelationWithScore: calculates the correlation of each
#' element with a numeric vector of the score, such as distance or
#' an immune score, puts value in RowData
#' @export
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
#' @export
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
#' @export
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





#' buildnetwork
#' Builds network using SingleCellExpression object and values
#' @export
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
#' @export
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

