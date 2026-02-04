# Calculate degree of missingness for each featurea and sample

calc_missingness() calculates the fraction of samples that are missing
data for each feature, and fraction of each feature that are missing for
each sample.

## Usage

``` r
calc_missingness(spe, group_name = NULL)
```

## Arguments

- spe:

  SpatialExperiment object

- group_name:

  Name of group to use as denominator

## Value

SpatialExperiment object with missingFeature column added to colData and
missingSample added to rowData

## Examples

``` r
library(spammR)
data(smallPancData)
data(pancMeta)
data(protMeta)
pooledPanc <- dplyr::bind_cols(smallPancData)
panc.spe <- convert_to_spe(pooledPanc, pancMeta, protMeta, 
          feature_meta_colname = "pancProts")
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
#> Note: Only mapping metadata for 2986 features out of 3000 data points
calc_missingness(panc.spe)
#> class: SpatialExperiment 
#> dim: 2986 27 
#> metadata(0):
#> assays(1): proteomics
#> rownames(2986): sp|A0A024RBG1|NUD4B_HUMAN sp|A0A096LP55|QCR6L_HUMAN ...
#>   sp|Q7Z2T5|TRM1L_HUMAN sp|Q7Z2W4|ZCCHV_HUMAN
#> rowData names(7): pancProts Entry ... PrimaryGeneName missingSample
#> colnames(27): 0_S_1_1 0_S_1_2 ... 7_S_3_2 7_S_3_3
#> colData names(17): Image x_coord ... sample_id missingFeature
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(0) :
#> imgData names(1): sample_id
```
