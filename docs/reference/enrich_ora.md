# Calculate functional or pathway enrichment

`enrich_ora()` calculates over-representation statistics (ORA) using an
interest list of genes from differential expression results in spammR
and gene sets (either the ones provided in spammR or user supplied) This
function uses results from `calc_spatial_diff_ex`. Interest list of
genes for ORA is obtained from spatialDiffEx results based on the
criteria specified in this function.It is a wrapper for the `leapR`
package which is required For ORA using an external or already defined
interest list of genes and gene sets, use leapR functions directly

## Usage

``` r
enrich_ora(
  spe,
  geneset,
  feature_column,
  pval_type_forThresh = "adjusted_pval",
  pval_thresh = 0.05,
  logFC_lowerThresh = NA,
  logFC_upperThresh = NA,
  geneset_name = "msigdb",
  sortResultsBy,
  comparison_name = ""
)
```

## Arguments

- spe:

  SpatialExperiment object containing spatial omics data and spatial
  diffex results

- geneset:

  in GMT format

- feature_column:

  Column of rowData that maps to gene set

- pval_type_forThresh:

  Choose from "adjusted_pval" or "pval". Type of p-value that should be
  used for filtering statistically significant results. Default is
  adjusted p-value for multiple hypotheses correction.

- pval_thresh:

  value to use for filtering based on pval_type_forThreshold. Default is
  0.05. Values less than pval_thresh will be kept.

- logFC_lowerThresh:

  Lower threshold for log Fold Change, to be used for filtering
  spatialDiffEx results. Default is NA

- logFC_upperThresh:

  Upper threshold for log Fold Change, to be used for filtering
  spatialDiffEx results. Default is NA

- geneset_name:

  Name of geneset provided

- sortResultsBy:

  For sorting ORA results, choose from the following column names:
  "BH_pvalue" (default)

- comparison_name:

  Example: "RSPv_vs_others" Text to indicate in results data frame,
  which spatial groups were compared for the interest list of genes

## Value

A dataframe containing results from over-representation analysis of
members of gene sets in the interest list of genes based on filtering
criteria above.

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
diffex.spe <- calc_spatial_diff_ex(panc.spe, category_col = "IsletOrNot")
#> Warning: Partial NA coefficients for 2 probe(s)
#> We found 0 features with a logFC greater than 1 and 
#>                  an ajusted p-value less than 0.05
library(leapR)
data("krbpaths")
ora.res <- enrich_ora(diffex.spe, geneset = krbpaths, 
         geneset_name = "krbpaths", feature_column = "PrimaryGeneName")
```
