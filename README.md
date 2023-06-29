
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

This is a basic example which demonstrates how to organize spatial omics data into a SpatialExperiment object, which is the required format for input data in SpammR. After constructing a SpatialExperiment object, we then demonstrate the use of different functions in SpammR. 

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
pancMeta$Xcoord = as.numeric(pancMeta$Xcoord)
pancMeta$Ycoord = as.numeric(pancMeta$Ycoord)
x <- pancMeta$Xcoord
y <- pancMeta$Ycoord
spat_xy = as.matrix(data.frame(x,y))
rownames(spat_xy) = rownames(pcolD)
##initialize the SpatialExperiment object
spe_multiImages <-SpatialExperiment(assays=list(logcounts=as(as.matrix(pancData),'dgCMatrix')),
                       colData=pcolD,
                       rowData=prowD,
                       spatialCoords = spat_xy,
                       sample_id = rownames(pcolD))
spe_multiImage

# Make another SpatialExperiment object, which only has data for one image (I picked image 0)
# This is to ensure that we don't have multiple samples with the same x,y coordinates
# Simulating a non-grid experiment dataset
pick_rows = which(pcolD$Image==0)
pcolD2 = pcolD[pick_rows,]
x = pancMeta[pick_rows,"Xcoord"]
y = pancMeta[pick_rows,"Ycoord"]
spat_xy_nongrid = as.matrix(data.frame(x,y))
rownames(spat_xy_nongrid) = rownames(pcolD2)
pancData_sub = pancData[,colnames(pancData) %in% rownames(pcolD2)]
spe_singleImg = SpatialExperiment(assays=list(logcounts=as(as.matrix(pancData_sub),'dgCMatrix')),
                                colData=pcolD2,
                                rowData=prowD,
                                spatialCoords = spat_xy_nongrid,
                                sample_id = rownames(pcolD2))
spe_singleImg

## run differential expression

diffex<-spatialDiffEx(spe,column='IsletOrNot',vals=c('Islet','NonIslet'))

## basic example code
```

