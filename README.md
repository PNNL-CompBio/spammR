
# spammR

<!-- badges: start -->
<!-- badges: end -->

The goal of spammR is to provide tools for generic analysis of omics data that is measured in spatial context. For more details and examples please explore our [recent publication on BioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.26.672472v1). 

## Installation

You can install the development version of spammR like so:

``` r
library(devtools)
devtools::install_github('pnnl-compbio/spammR',build_vignettes = TRUE,force=TRUE)
```

## Example

This is a basic example which demonstrates how to organize spatial omics data into a SpatialExperiment object, which is the required format for input data in SpammR. After constructing a SpatialExperiment (SPE or spe) object, we then demonstrate the use of different functions in SpammR. 

### Organize data into a SpatialExperiment object

``` r
library(spammR)

## load test data
data('smallPancData')
data('pancMeta')
data('protMeta')
```

## build into SpatialExperiment object

``` r
img0.spe<-convert_to_spe(smallPancData$Image_0,pancMeta,protMeta,feature_meta_colname='pancProts',image_files=system.file("extdata",'Image_0.png',package='spammR'),image_samples_common_identifier='Image0',spatialCoords_colnames=c('x_pixels','y_pixels'),samples_common_identifier = 'Image0',image_ids='with_grid')

```

## Produce a spatial heatmap for a chosen feature from the data

This function plots feature values for a single feature as a heatmap on an x-y coordinate system
First, define parameters needed for this plotting. Let's look at the spatial expression of the insulin protein.

 
1a. Basic spatial heatmap where all parameters are specified by the user; non-interactive
``` r
res = spatial_heatmap(img0.spe, feature='INS', sample_id='Image0', image_id='with_grid', spatial_coord_names=c('x_pixels','y_pixels'), spot_size=unlist(colData(img0.spe)[1,c('spot_width','spot_height')]), image_boundaries=unlist(colData(img0.spe)[1,c('x_origin','y_origin','x_max','y_max')]),label_column='IsletOrNot', interactive=FALSE)
res
```
1b. Same thing as 1a but now an interactive plot; hovering over a square gives it's coordinates, label and colored value
``` r
res = spatial_heatmap(img0.spe, feature='INS', sample_id='Image0', image_id='with_grid', spatial_coord_names=c('x_pixels','y_pixels'), spot_size=unlist(colData(img0.spe)[1,c('spot_width','spot_height')]), image_boundaries=unlist(colData(img0.spe)[1,c('x_origin','y_origin','x_max','y_max')]),label_column='IsletOrNot', interactive=TRUE)
res
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

