# Impute missing values based on spatial coordinates

`impute_spe()` carries out imputation for missing data in the primary
assay using a specified method from a range of methods. Returns
additional assay in the same SPE

## Usage

``` r
impute_spe(
  spe,
  assay_name = NULL,
  method = NULL,
  group_colname,
  k = NULL,
  protein_missingness = NULL
)
```

## Arguments

- spe:

  SPE containing data to be imputed.

- assay_name:

  name of assay with data to be imputed

- method:

  Method of imputation to be used. See details.

- group_colname:

  Column name in metadata that specifies the group information to use
  for group_mean or knn_group. Example: ROI_abbreviation.

- k:

  K value to be used for k-nearest neighbor imputation

- protein_missingness:

  Proportion of samples allowed to have missing data for a protein in
  the given imputation method. Example: When method="global_mean," an
  protein_missingness of 0.5 indicates that any protein missing data for
  more than 50% of samples across the entire spatial tissue covered by
  all samples will be excluded from the imputation method algorithm, and
  that protein's missing values will not be imputed. When
  method="group_mean," then protein_missingness of 0.5 indicates that a
  protein must have data for at least 50% of samples in the specified
  group to be used in the imputation algorithm and to be imputed.

## Value

An SPE with an 'imputed' assay with the appropriate imputation called

## Details

Methods options and descriptions:

- zero: replace missing values with 0

- median : replace missing values with global median per protein

- median_half : replace missing values with 1/2 global median per
  protein

- mean: replace missing values with global mean per protein

- group_mean: replace missing values with mean per group, e.g. group
  (example: ROI) for each protein

- knn: imputation based on k-nearest neighbors, with proteins as
  neighbors, based on data from all samples across all groups. NOTE:
  There will still be NA values if the protein is not expressed in this
  group.

- group_knn: imputation based on k-nearest neighbors, with proteins as
  neighbors, based on data from specified group (e.g. ROI, tissue).
  NOTE: There will still be NA values if the protein is not expressed in
  this group.

- spatial_knn: imputation based on k-nearest neighors in space

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
  feature_meta_colname = "pancProts",
  sample_id = ""
)
#> Spatial object created without spatial coordinate 
#>          column names provided. Distance based analysis will not be enabled.
#> Note: Only mapping metadata for 2986 features out of 3000 data points
# we can try two imputation methods and compare the difference
res <- impute_spe(pooled.panc.spe, method = "mean")
res2 <- impute_spe(pooled.panc.spe, method = "group_mean", 
                  group_colname = "Image")
mean(assay(res, "imputed") - assay(res2, "imputed"), na.rm = TRUE)
#> [1] 0.003118943
```
