
# spammR

<!-- badges: start -->
<!-- badges: end -->

The goal of spammer is to provide tools for generic analysis of omics data that is measured in spatial context. 

## Installation

You can install the development version of spammer like so:

``` r
library(devtools)
devtools::github('pnnl-compbio/spammer')
```

## Example

This is a basic example which demonstrates how to organize spatial omics data into a SpatialExperiment object, which is the required format for input data in SpammR. After constructing a SpatialExperiment object, we then demonstrate the use of different function in SpammR. 

``` r
library(spammR)
## load test data
data('pancData')
data('pancMeta')

## build into SpatialExperiment object


rowData <- data.frame(protein=rownames(pancData))
rownames(rowData)<-rownames(pancData)
prowD <- rowData[rownames(pancData),]
pcolD <- pancMeta[colnames(pancData),]


##initialize hte objects
spe<-SpatialExperiment(assays=list(logcounts=as(as.matrix(pancData),'dgCMatrix')),
                                  colData=pcolD,
                                  rowData=prowD)

## run differential expression

diffex<-spatialDiffEx(spe,column='IsletOrNot',vals=c('Islet','NonIslet'))

## basic example code
```

