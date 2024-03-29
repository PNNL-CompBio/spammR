
# spammR

<!-- badges: start -->
<!-- badges: end -->

The goal of spammR is to provide tools for generic analysis of omics data that is measured in spatial context. 

## Installation

You can install the development version of spammR like so:

``` r
library(devtools)
devtools::install_github('pnnl-compbio/spammR')
```

## Example

This is a basic example which demonstrates how to organize spatial omics data into a SpatialExperiment object, which is the required format for input data in SpammR. After constructing a SpatialExperiment (SPE or spe) object, we then demonstrate the use of different functions in SpammR. 

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
# The example data that we're using here represents data from multiple images in an experiment.
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
spat_xy_singleImg = as.matrix(data.frame(x,y))
rownames(spat_xy_singleImg) = rownames(pcolD2)
pancData_sub = pancData[,colnames(pancData) %in% rownames(pcolD2)]
spe_singleImg = SpatialExperiment(assays=list(logcounts=as(as.matrix(pancData_sub),'dgCMatrix')),
                                colData=pcolD2,
                                rowData=prowD,
                                spatialCoords = spat_xy_singleImg,
                                sample_id = rownames(pcolD2))
spe_singleImg
```


### 1. Produce a spatial heatmap for a chosen feature from the data
SpammR function to be used: spatialPlot_feature(...)

This function plots feature values for a single feature as a heatmap on an x-y coordinate system

First, define parameters needed for this plotting. Let's look at the spatial expression of the insulin protein.
``` r
"sp|P01308|INS_HUMAN" #insulin protein
metric_lab = "Protein abundance measure" # Metric represented by color scale; this will be used as the legend label
label_col = "sample_id" # name of the column in colData(spe) to be used for labeling sample locations
```
spatialPlot_feature() accepts the following input parameters:
- spe: A SpatialExperiment object containing data to be plotted. Must have spatial coordinates for every sample stored in spatialCoords(spe)
- feature: Feature (example: protein) whose values are to be plotted. Should be a row name in rowData(spe)
- (optional) metric_display: legend label for the color legend. If not specified, defaults to "Protein abundance measure"
- (optional) label_column: name of the column in colData(spe) to be used for labeling sample locations. If not specified, no labels are shown for samples.
- (optional) interactive: Boolean value (TRUE/FALSE) to make the plot have interactive mouse hovering spatially. Default is an interactive plot.
  
1a. Basic spatial heatmap where all parameters are specified by the user; non-interactive
``` r
spatialPlot_feature(spe_singleImg, feat, metric_lab, label_col, interactive = FALSE)
# Note: Grid squares can only be labeled when the plot is not interactive.
```
1b. Same thing as 1a but now an interactive plot; hovering over a square gives it's coordinates, label and colored value
``` r
spatialPlot_feature(spe_singleImg, feat, metric_lab, label_col)
# or
spatialPlot_feature(spe_singleImg, feat, metric_lab, label_col, interactive = TRUE)
``` 
1c. A spatial heatmap without labels for grid squares
``` r
spatialPlot_feature(spe_singleImg, feat, metric_lab, NA)
# or
spatialPlot_feature(spe_singleImg, feat, metric_lab)
```
1d. If the user wants to use the default legend label defined in the function "Protein abundance measure" and no labels
``` r
spatialPlot_feature(spe_singleImg, feat)
``` 
1e. Same thing as 1d but now specify which column to use for labeling grid squares
``` r
spatialPlot_feature(spe_singleImg, feat, label_column = label_col)
``` 


### 2. Differential expression analysis comparing two types of samples in the data given in the SPE object
2a. Run differential expression using spatialDiffEx() function, to compare "Proximal" vs. "Distal" samples in the Pancreas dataset stored in spe_singleImg

The spatialDiffEx() function requries the following input parameters:
- spe:  Spatial Experiment object
- logTransformed: Boolean indicating whether the data given in spe is log2 transformed (TRUE) or not(FALSE)
- category_col:  Name of the column that specifies category of each sample. Example: "IsletStatus"
Categories from category_col will be compared in the differential expression analysis
- compare_vals: A vector containing names of categories from category_col to be compared. example: c('Proximal','Distal')

``` r
diffex_spe <-spatialDiffEx(spe_singleImg,TRUE,category_col ='IsletStatus',compare_vals=c('Proximal','Distal'))
```
The output diffex_spe is a SpatialExperiment object containing results from differential expression analysis, in addition to what was already present in the input spe

2b. How to access the differential expression results in the returned SPE object?
``` r
rowData(diffex_spe)
```
2c. Produce a volcano plot to visualize the results from differential experession analysis

In spammR, we provide functions to create a typical volcano plot i.e. a plot shoing -log10(pvalue) vs. log2(fold change), along with annotations of points based on criteria of significance provided by the user.

There are two functions in spammR that can be used for this depending on the data structure in which the differential expression results are stored. The functions are volcanoPlot_DiffExSpe() and volcanoPlot_DiffExResults(). 

The input parameters needed for both functions are the same except for the first parameter (the object containing the differential expression results).

First, define the needed input parameters for both functions
``` r
logFC_colname = "IsletStatus.limma.logFC"
pval_colname = "IsletStatus.limma.P.Value"
pval_corrected = TRUE
title = "Differential abundance results for Distal vs. Proximal samples; Pancreas dataset spe_SingleImage"
thresh = 0.05
sigLabel_colname = "X"
```
If the differential expression results are stored in a data frame, then use volcanoPlot_DiffExResults().
``` r
diffEx_df = rowData(diffex_spe)
volcanoPlot_DiffExResults(diffEx_df, logFC_colname, pval_colname, pval_corrected, title, thresh, sigLabel_colname)
```
If the differential expression results are stored in an SPE object, then use volcanoPlot_DiffExSpe(). volcanoPlot_DiffExSpe() calls on volcanoPlot_DiffExResults():
``` r
volcanoPlot_DiffExSpe(diffex_spe, logFC_colname, pval_colname, pval_corrected, title, thresh, sigLabel_colname)
```
## basic example code
```

