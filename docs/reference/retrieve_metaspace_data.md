# Retrieve data objects from metaspace

`retrieve_metaspace_data`Collects data from http://metaspace2020.org by
project id and returns it in SpatialExperiment object

## Usage

``` r
retrieve_metaspace_data(
  project_id = "2024-02-15_20h37m13s",
  fdr = 0.2,
  assay_name = "lipids",
  sample_id = "sample",
  rotate = FALSE,
  drop_zeroes = TRUE,
  y_offset = 0,
  x_offset = 0
)
```

## Arguments

- project_id:

  Identifier for project

- fdr:

  FDR value above which to collect data

- assay_name:

  Name of assay to include in object

- sample_id:

  Name of sample

- rotate:

  Set to TRUE if x and y coordinates need to be swapped

- drop_zeroes:

  Set to TRUE to drop zeroes and set the values to NAs. If False the
  zeroes will be included in the image data.

- y_offset:

  Number of pixels to adjust y based on image

- x_offset:

  Number of pixels to adjust x based on image

## Value

a `SpatialExperiment` object with the metaspace data and coordinates.

## Details

This relies on underlying python code and environment to run, including
the metaspace2020 package, the numpy package, and the pandas package.

## Examples

``` r
# \donttest{
rdat <- retrieve_metaspace_data(project_id = "2024-02-15_20h37m13s",
                                fdr = 0.2,
                                assay_name = 'lipids',
                               sample_id = 'sample')
#> Downloading ion data from metaspace...
#> Warning: index contains duplicated values: row names not set
head(rowMeans(assay(rdat,'lipids')))
#> C46H78NO10P-H-  C47H83O13P-H-    C20H32O2-H-    C22H32O2-H-    C16H30O2-H- 
#>             NA             NA             NA             NA             NA 
#> C42H80NO10P-H- 
#>             NA 
# }
```
