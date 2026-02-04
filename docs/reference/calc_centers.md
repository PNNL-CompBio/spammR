# Calculate centers of spots/samples for distance based analysis

calc_centers calculates the centers of spots for numerical
analysis.Currently not used.

## Usage

``` r
calc_centers(spe)
```

## Arguments

- spe:

  SpatialExperiment object containing omics data

## Value

Data frame of new coordinates

## Examples

``` r
data(pancMeta)
data(protMeta)
data(smallPancData)
img0.spe <- convert_to_spe(smallPancData$Image_0,
  pancMeta,
  protMeta,
  feature_meta_colname = "pancProts",
  image_files = system.file("extdata", "Image_0.png", package = "spammR"),
  sample_id = "Image0",
  spatial_coords_colnames = c("x_pixels", "y_pixels"),
  image_sample_ids = "Image0", image_ids = "Image0"
)
#> Note: Only mapping metadata for 2986 features out of 3000 data points
## can't actually call internal functuon
```
