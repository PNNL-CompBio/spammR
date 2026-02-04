# Reduce a high resolution image measurement (origin) to a lower resolution measurement (target)

`spat_reduce()` Creates a single SpatialExperiment object with all
measurements mapped to the same coordinate space for plotting and
analysis

## Usage

``` r
spat_reduce(
  spe_target,
  spe_origin,
  origin_assay = "lipids",
  target_dims = c("spot_width", "spot_height"),
  replace_nas = FALSE,
  stat = "mean"
)
```

## Arguments

- spe_target:

  Low resolution image measurements

- spe_origin:

  High-resolution image measurements to be mapped to target

- origin_assay:

  Name of assay in origin assay to be reduced

- target_dims:

  Names of width and height of target spots, respectively

- replace_nas:

  Set to TRUE if we want to replace NA values with zeroes. This will
  increase the number of data points in common for correlation
  calculations.

- stat:

  Statistic used to summarize regions, choices are currently 'mean' and
  'median'

## Value

A copy of the spe_target object with the `altExp` slot containing a
reduced version of the origin

## Details

This function allows us to map two distinct omics measurements to the
same coordinate system. It takes on 'origin' measurement and summarizes
the values (typically via mean) in the ROI in a 'target' measurement so
that the source and target values can be compared.
