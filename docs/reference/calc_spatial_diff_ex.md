# Identify differentially abundant features across image

`calc_spatial_diff_ex()` Calculates differential expression analysis
using annotations in a SpatialExperiment object

## Usage

``` r
calc_spatial_diff_ex(
  spe,
  assay_name = "proteomics",
  count_based = FALSE,
  log_transformed = FALSE,
  category_col,
  compare_vals
)
```

## Arguments

- spe:

  Spatial Experiment object containing data to be used for differential
  expression analysis

- assay_name:

  Name of the dataset stored in the spe object, that is to be used for
  the differential expression analysis. Example: znormalized_log2

- count_based:

  Set to TRUE of the data are count based, e.g. RNA-Seq

- log_transformed:

  Is the data given in spe log2 transformed TRUE or FALSE

- category_col:

  Name of the column that specifies category of each sample. Example:
  "IsletOrNot" \#Categories from `category_col` will be compared in the
  differential expression analysis

- compare_vals:

  A vector containing names of categories from category_col to be
  compared. Only required if there are more than two values in
  `category_col`

## Value

A Spatial Experiment object containing differential expression results,
stored in `rowData(diffEx.spe)` and `assays(diffEx.spe)` which contains
the dataset on which differential expresssion analysis was carried out

## Examples

``` r
data(smallPancData)
data(pancMeta)
data(protMeta)
pooledData <- dplyr::bind_cols(smallPancData)
pooled.panc.spe <- convert_to_spe(pooledData,
  pancMeta,
  protMeta,
  feature_meta_colname = "pancProts"
)
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
#> Note: Only mapping metadata for 2986 features out of 3000 data points
diffex.spe <- calc_spatial_diff_ex(pooled.panc.spe,
  category_col = "IsletOrNot"
)
#> Warning: Partial NA coefficients for 2 probe(s)
#> We found 0 features with a logFC greater than 1 and 
#>                  an ajusted p-value less than 0.05
```
