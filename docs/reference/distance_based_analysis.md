# Identify multiomic features correlated with an distance to an image feature

distance_based_analysis: Identifies proteins/features that show a strong
correlation between distance from a specified ROI's samples and
protein/feature abundance differences between samples.

## Usage

``` r
distance_based_analysis(
  spe,
  assay_name,
  spotHeightCol = "spot_height",
  spotWidthCol = "spot_width",
  sampleCategoryCol,
  sampleCategoryValue,
  featuresNameCol,
  corr_type = "spearman",
  corr_thresh = 0.5,
  min_samples = 5,
  allowOverlaps = TRUE
)
```

## Arguments

- spe:

  SpatialExperiment object containing spatial omics data

- assay_name:

  Name of the assay stored in spe that is to be used for distance based
  analysis. Example: "znormalized_log2"

- spotHeightCol:

  Column containing height of spot

- spotWidthCol:

  Column containing width of spot

- sampleCategoryCol:

  Column name in metadata (colData(spe)) that should be used for
  selecting samples of certain type, to define the "origin" region for
  distance based analysis

- sampleCategoryValue:

  Sample category to be used for defining the "origin" region for
  distance based analysis

- featuresNameCol:

  Name of column containing features (example: proteins) in
  rowData(spe). It is assumed that the data provided in assay(spe) is in
  the same order as the order in which the features are listed under the
  featuresNameCol

- corr_type:

  Choose from "pearson" (default), "spearman." Correlation method to be
  used for calculating correlation between distance between samples and
  protein abundance differences. Both types of correlation provide a
  measure of monotonic association between two variables. Pearson is
  better suited for linear relationships while Spearman is better for
  non-linear monotonic relationships.

- corr_thresh:

  Minimum correlation value to be used for identifying proteins that
  have a correlation between protein abundance differences and distance
  between samples. Values greater than or equal to this theshold will be
  used.

- min_samples:

  of a minimum number of sample points for calculating correlation. For
  proteins with less than this number of sample points, correlation
  value is reported as NA.

- allowOverlaps:

  allow overlaps of regions

## Value

a `SpatialExperiment` object with the distance based data in the colData

## Examples

``` r
data(pancMeta)
data(protMeta)
data(smallPancData)
img0.spe <- convert_to_spe(smallPancData$Image_0,
  pancMeta,
  protMeta,
  sample_id = "Image0",
  spatial_coords_colnames = c("x_pixels", "y_pixels"),
  feature_meta_colname = "pancProts",
  image_files = system.file("extdata", "Image_0.png", package = "spammR"),
  image_ids = "Image0"
)
#> Note: Only mapping metadata for 2986 features out of 3000 data points

img0.spe <- distance_based_analysis(img0.spe,
  "proteomics",
  sampleCategoryCol = "IsletOrNot",
  sampleCategoryValue = "Islet"
)
```
