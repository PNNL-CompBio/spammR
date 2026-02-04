# Create a `SpatialExperiment` object from data

`convert_to_spe()` Puts omics data (omics measurements, metadata and
image(s) corresponding to samples' tissue) into a SpatialExperiment
(SPE) object. Most spammR functions require the input data to be an SPE
object.

## Usage

``` r
convert_to_spe(
  dat,
  sample_meta,
  feature_meta,
  assay_name = "proteomics",
  sample_id = "sample",
  sample_colname = NULL,
  feature_meta_colname = NULL,
  spatial_coords_colnames = NULL,
  rescale_image = FALSE,
  image_files = NULL,
  image_ids = NULL,
  image_sample_ids = rep(sample_id, length(image_files))
)
```

## Arguments

- dat:

  Matrix or data frame with omics measurements. Rows are features,
  columns are samples

- sample_meta:

  Data frame of metadata with either `sample_colname` or rownames as
  samples (if empty)

- feature_meta:

  Data frame of metadata with either `sample_colname` or rownames as
  samples (if empty)

- assay_name:

  Name to be given to the data in omics_measurements_file. Example:
  "abundance", "log2", "znormalized_log2" or any other descriptive name

- sample_id:

  Name of sample, defaults to "sample01"

- sample_colname:

  Column name in `sample_meta` table whose entries contain sample
  identifiers provided as column names in `dat`.

- feature_meta_colname:

  Name of column in `feature_meta`, that is to be used for identifying
  features. If missing defaults to rownames, which should match rownames
  of `dat`

- spatial_coords_colnames:

  A list containing names of columns in `meta_dat` that are spatial
  coordinates. Default value is NULL in which case no spatial
  coordinates are entered

- rescale_image:

  A boolean set to true if you need to rescale image. Recommended when
  coordinates are pulled from metaspace and are known to fully cover the
  image

- image_files:

  A list containing paths of image files to be stored in the
  SpatialExperiment object. More images can be added later, without
  using this function.

- image_ids:

  A list containing image names/identifiers for image paths provided in
  image_files

- image_sample_ids:

  A list of sample identifiers for each of the images provided in
  image_files, defaults to the value from sample_id

## Value

spe.out a SpatialExperiment (SPE) object that contains all data and
image(s). Ready to be used as input in spammR functions that require SPE
object as input.

## Examples

``` r
data(pancMeta)
data(protMeta)
data(smallPancData)
# We can put all samples into the same object (for statistical power)
pooledData <- dplyr::bind_cols(smallPancData)
pooled.panc.spe <- convert_to_spe(pooledData,
  pancMeta,
  protMeta,
  feature_meta_colname = "pancProts"
)
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
#> Note: Only mapping metadata for 2986 features out of 3000 data points

# or we can add the inmage to a single data capture
img0.spe <- convert_to_spe(smallPancData$Image_0,
  pancMeta,
  protMeta,
  sample_id = "Image0",
  feature_meta_colname = "pancProts",
  image_files = system.file("extdata", "Image_0.png", package = "spammR"),
  spatial_coords_colnames = c("x_pixels", "y_pixels"),
  image_ids = "Image0"
)
#> Note: Only mapping metadata for 2986 features out of 3000 data points
```
