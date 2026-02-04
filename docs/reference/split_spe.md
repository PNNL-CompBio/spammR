# Divide SPE object into smaller objects by column features

`split_spe()` splits an SPE object into smaller objects by a column
feature. Ideally this should be in the SpatialExperiment class or some
utility package.

## Usage

``` r
split_spe(spe, split_colname, assay_name = NULL)
```

## Arguments

- spe:

  SpatialExperiment object containing spatial omics data and spatial
  diffex results

- split_colname:

  Column of rowData that maps to gene set

- assay_name:

  Optional name of assay to use in split data

## Value

A list of SpatialExperiment objects containing a subset of the data

## Examples

``` r
data(smallPancData)
data(pancMeta)
data(protMeta)
pooledPanc <- dplyr::bind_cols(smallPancData)
panc.spe <- convert_to_spe(pooledPanc, pancMeta, protMeta, 
          feature_meta_colname = "pancProts")
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
#> Note: Only mapping metadata for 2986 features out of 3000 data points
split_list <- split_spe(panc.spe, split_colname = "Image")
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
```
