
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

### 0. Organize data into a SpatialExperiment object

``` r
library(spammR)
library(SpatialExperiment)
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
# The example datat that we're using represents data from multiple images in an experiment.
# In the intiial version of SpammR, we will not worry about providing tools for analyzing data from multiple images.
# In the later versions of SpammR, we will add functionality to deal with data from multiple images.
spe_multiImages <-SpatialExperiment(assays=list(logcounts=as(as.matrix(pancData),'dgCMatrix')),
                       colData=pcolD,
                       rowData=prowD,
                       spatialCoords = spat_xy,
                       sample_id = rownames(pcolD))
spe_multiImages

# Make another SpatialExperiment object, which only has data for one image (I picked image 0)
# This is to ensure that we don't have multiple samples with the same x,y coordinates
# We will use this object for demonstration of the spatial plotting functionality in SpammR for data from a single image.
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
```


### 1. Produce a spatial heatmap for a chosen feature from the data
SpammR function to be used: spatialPlot_feature(...)

Plot feature values for a single feature on an x-y coordinate system

Define parameters needed for this plotting. Example:
``` r
feat = "sp|A0A024RBG1|NUD4B_HUMAN"
metric_lab = "Protein abundance measure" # Metric represented by color scale; this will be used as the legend label
label_col = "sample_id" # name of the column to be used for labeling sample locations
```
1a. Basic spatial heatmap where all parameters are specified by the user; non-interactive
``` r
spatialPlot_feature(spe_singleImg, feat, metric_lab, label_col, interactive = FALSE) # Squares can only be labeled when the plot is not interactive.
```
1b. Same thing but now an interactive plot; hovering over a square gives it's coordinates, label and colored value
``` r
spatialPlot_feature(spe_nongrid, feat, metric_lab, label_col)
# or
spatialPlot_feature(spe_nongrid, feat, metric_lab, label_col, interactive = TRUE)
``` 
# 1c. A spatial heatmap without labels for grid squares
spatialPlot_feature(spe_nongrid, feat, metric_lab, NA)
# or
spatialPlot_feature(spe_nongrid, feat, metric_lab)
# 1d. If the user wants to use the default legend label defined in the function "Protein abundance measure" and no labels
spatialPlot_feature(spe_nongrid, feat)
# 1e. Same thing as 1d but now specify which column to use for labeling grid squares
spatialPlot_feature(spe_nongrid, feat, label_column = label_col) 


## run differential expression

diffex<-spatialDiffEx(spe,column='IsletOrNot',vals=c('Islet','NonIslet'))

## basic example code
```

