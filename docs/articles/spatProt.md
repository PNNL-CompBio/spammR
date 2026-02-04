# spammR Spatial Proteomics Example

## Getting started

To install the package currently you must install directly from `GitHub`
along with the leapR dependency as shown below. Before release we hope
to move to Bioconductor.

The `leapR` package is designed for flexible pathway enrichment and
currently must be installed before spammR.

      ##install if not already installed
      library(devtools)
      devtools::install_github('PNNL-CompBio/leapR')
      devtools::install_github('PNNL-CompBio/spammR')

Once the package is installed you load the library, including the test
data.

``` r
## load spammR
library(spammR)
```

## Collecting data to analyze

`spammR` enables the analysis of disparate sets of multiomic data:
image-based data and numerical measurements of omics data. It is
incredibly flexible as to the *type* of multiomic data. We assume each
omics measurement is collected in a single sample, and that there are
specific spatial coordinates for that sample in the image. We leverage
the `SpatialExperiment` object to store the data for each
image/measurement pair.

### Data overview and examples

The `spammR` package requires omics data with spatial coordinates for
the functions to run successfully. Here we describe the data required
and show examples.

#### Omics Measurement Data

`SpatialExperiment` can hold multiple omics measurements mapping to the
same sample identifier in different ‘slots’. This data can be a tabular
data frame or matrix with rownames referencing measurements in a
particular sample (e.g. gene, species) and column names representing
sample identifiers. An example of this can be found by loading
`pancDataList.rda` file from Figshare.

To evaluate the features of this package we are using pancreatic data
from [Gosline et al.]() that is captured using mass spectrometry
measured from 7 independent regions of a single human pancreas. Each
image is segmented into nine ‘voxels’, with one voxel per image
representing a cluster of islet cells.

``` r
download.file("https://api.figshare.com/v2/file/download/55158821",
              mode = "wb", quiet = TRUE, dest = "pdl.rda")
load("pdl.rda")
utils::head(pancDataList$Image_0[, 1:8])
```

    ##                            0_S_1_1  0_S_1_2  0_S_1_3  0_S_2_1  0_S_2_2  0_S_2_3
    ## sp|A0A024RBG1|NUD4B_HUMAN 13.06042 13.42317 12.42396 13.02470 12.56442 12.69023
    ## sp|A0A096LP55|QCR6L_HUMAN 15.10920 15.27460 15.16780 15.01030 15.46639 14.73712
    ## sp|A0AV96|RBM47_HUMAN     17.40246 17.29727 17.25559 17.34851 17.12866 17.17658
    ## sp|A0AVT1|UBA6_HUMAN      18.00653 18.43015 18.24663 18.17563 18.38961 18.29268
    ## sp|A0FGR8|ESYT2_HUMAN     16.59018 16.48890 16.50134 16.55334 16.32316 16.43397
    ## sp|A0MZ66|SHOT1_HUMAN     18.19277 18.73633 18.54485 18.20005 18.74041 18.70762
    ##                            0_S_3_1  0_S_3_2
    ## sp|A0A024RBG1|NUD4B_HUMAN       NA 13.40670
    ## sp|A0A096LP55|QCR6L_HUMAN 14.81792 15.63741
    ## sp|A0AV96|RBM47_HUMAN     17.27792 17.11678
    ## sp|A0AVT1|UBA6_HUMAN      18.12080 18.10220
    ## sp|A0FGR8|ESYT2_HUMAN     16.43682 16.20515
    ## sp|A0MZ66|SHOT1_HUMAN     18.60198 18.62946

``` r
file.remove("pdl.rda")
```

    ## [1] TRUE

Here the rownames represent protein identifiers and the column names
represent individual samples. Each element of the list contains the
measurements from a different sample:

``` r
print(length(pancDataList))
```

    ## [1] 7

``` r
head(pancDataList[[2]][, 1:8])
```

    ##                            1_S_1_1  1_S_1_2  1_S_1_3  1_S_2_1  1_S_2_2  1_S_2_3
    ## sp|A0A024RBG1|NUD4B_HUMAN       NA       NA       NA       NA       NA       NA
    ## sp|A0A096LP55|QCR6L_HUMAN       NA       NA       NA       NA       NA       NA
    ## sp|A0AV96|RBM47_HUMAN     17.68866 17.58076 17.51335 17.65900 17.52926 17.44999
    ## sp|A0AVT1|UBA6_HUMAN      18.02235 18.35493 18.03651 18.09462 18.08796 18.04204
    ## sp|A0FGR8|ESYT2_HUMAN     17.50177 17.51525 17.34338 17.41622 17.34254 17.46062
    ## sp|A0MZ66|SHOT1_HUMAN     18.57445 18.67262 18.83672 18.69094 18.48785 18.75190
    ##                            1_S_3_1  1_S_3_2
    ## sp|A0A024RBG1|NUD4B_HUMAN       NA       NA
    ## sp|A0A096LP55|QCR6L_HUMAN       NA       NA
    ## sp|A0AV96|RBM47_HUMAN     17.37103 17.50367
    ## sp|A0AVT1|UBA6_HUMAN      18.21052 17.94707
    ## sp|A0FGR8|ESYT2_HUMAN     17.31495 17.53217
    ## sp|A0MZ66|SHOT1_HUMAN     18.82805 18.48856

This list is used below in our analysis examples.

#### Sample metadata

The samples metadata table contains mappings between samples and
metadata. An example can be found in `data(pancMeta)`. Most importantly
we require the image mapping information, which includes: - *Image
coordinates:* to map the image to a coordinate space we need to know the
`x_origin`, and `y_origin` (assumed to be zero) as well as `x_max` and
`y_max`, which is the top right of the image. The package plots the
*entire* image so specifying these coordinates ensures that all other
points are properly mapped. - *Sample coordinates:* Each sample has its
own `x_coord` and `y_coord`. - *Spot size:* `spot_height` and
`spot_width`.

``` r
data(pancMeta)
head(pancMeta)
```

    ##         Image x_coord y_coord IsletStatus IsletOrNot Plex Grid.Number x_pixels
    ## 0_S_3_1     0       3       1    Proximal   NonIslet 127N           1      475
    ## 0_S_2_1     0       2       1       Islet      Islet 128N           2      380
    ## 0_S_1_1     0       1       1    Proximal   NonIslet 127C           3      285
    ## 0_S_3_2     0       3       2    Proximal   NonIslet 128C           4      475
    ## 0_S_2_2     0       2       2    Proximal   NonIslet 129N           5      380
    ## 0_S_1_2     0       1       2    Proximal   NonIslet 129C           6      285
    ##         y_pixels x_origin y_origin x_max y_max spot_width spot_height
    ## 0_S_3_1      170        0        0   860   725         90         140
    ## 0_S_2_1      170        0        0   860   725         90         140
    ## 0_S_1_1      170        0        0   860   725         90         140
    ## 0_S_3_2      315        0        0   860   725         90         140
    ## 0_S_2_2      315        0        0   860   725         90         140
    ## 0_S_1_2      315        0        0   860   725         90         140

This metadata contains information for all 7 images, so we do not need a
separate metadata file for each image, the `convert_to_spe` function
will simply take the metadata relevant to the data file.

#### Image files

There can be multiple image files associated with a single set of omics
measurements. Currently we have tested working with files in `png`
format. Each image we have is stained so that we can identify the Islet
cells. Each image also has a grid superimposed to show where the sample
measurements came from. The grid is not necessary, of course, but can
help alibrate the coordinates.

``` r
library(cowplot)

cowplot::ggdraw() + cowplot::draw_image(system.file("extdata",
                                                    "Image_1.png",
                                                    package = "spammR"))
```

Now we can use this image and others to visualize omics data.

#### Omics metadata

The last set of metadata relates to the `rows` of the omics measurement
data. When using gene-based data, this will be the genes or proteins in
the dataset. When using metagenomics, this will refer to the species.
One column of this table must uniquely map to the rownames of the omics
data.

``` r
data(protMeta)
head(protMeta[, c("pancProts", "EntryName", "PrimaryGeneName")])
```

    ##                   pancProts   EntryName PrimaryGeneName
    ## 1 sp|A0A024RBG1|NUD4B_HUMAN NUD4B_HUMAN          NUDT4B
    ## 2 sp|A0A096LP55|QCR6L_HUMAN QCR6L_HUMAN          UQCRHL
    ## 3     sp|A0AV96|RBM47_HUMAN RBM47_HUMAN           RBM47
    ## 4      sp|A0AVT1|UBA6_HUMAN  UBA6_HUMAN            UBA6
    ## 5     sp|A0FGR8|ESYT2_HUMAN ESYT2_HUMAN           ESYT2
    ## 6     sp|A0MZ66|SHOT1_HUMAN SHOT1_HUMAN           SHTN1

This data helps us find better gene identifiers.

### Loading data into spatial experiment object.

Now that we have all the data loaded we can build a `SpatialExperiment`
object either using ALL samples or just the samples in a single image.
We can pool all the data for more statistical power.

``` r
pooledData <- dplyr::bind_cols(pancDataList)
pooled.panc.spe <- convert_to_spe(pooledData, ## pooled data table
  pancMeta, ## pooled metadata
  protMeta, ## protein identifiers
  feature_meta_colname = "pancProts", # column name
)
```

    ## Spatial object created without spatial coordinate 
    ##          column names provided. Distance based analysis will not be enabled.

    ## Note: Only mapping metadata for 6662 features out of 6693 data points

``` r
print(pooled.panc.spe)
```

    ## class: SpatialExperiment 
    ## dim: 6662 63 
    ## metadata(0):
    ## assays(1): proteomics
    ## rownames(6662): sp|A0A024RBG1|NUD4B_HUMAN sp|A0A096LP55|QCR6L_HUMAN ...
    ##   sp|Q9Y3M8|STA13_HUMAN sp|Q9Y6X3|SCC4_HUMAN
    ## rowData names(6): pancProts Entry ... GeneNames PrimaryGeneName
    ## colnames(63): 0_S_1_1 0_S_1_2 ... 3_S_3_2 3_S_3_3
    ## colData names(16): Image x_coord ... spot_height sample_id
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):
    ## spatialCoords names(0) :
    ## imgData names(0):

We can also create a list of `SpatialExperiment` objects, one for each
of the 3 images we have.

Now we can use these individual image objects or the combined ‘pooled’
object for analysis.

### Spatial data with image

Here we loop over all of the images in `imglist` to plot the expression
of the insulin protein in each image. We expect insulin (or INS) to be
highest in voxels containing islet cells, which we label using the
`label_column` ‘IsletOrNot’ which was loaded into the metadata for us.

``` r
allimgs <- lapply(imglist, function(x) {
    spe <- img.spes[[x]]
    res <- spatial_heatmap(spe,
      feature = "INS",
      feature_type = "PrimaryGeneName",
      sample_id = x,
      image_id = "with_grid",
      label_column = "IsletOrNot",
      interactive = FALSE
  )
  return(res)
})

allimgs[[2]]
```

![](spatProt_files/figure-html/spatial%20data-1.png)

To go further and visualize entire pathways we need to first identify
which groups of proteins are of interest using a more unsupervised
approach.

## Expression and pathway analysis

Now that we have the ability to overlay omic measurements with image
ones, we can identify new features to plot and visualize them. First we
can employ standard differential expression approaches using the voxel
labels and the `limma` pathway.

### Differential expression

First we want to identify specific proteins that are up-regulated in the
islet cells (or regions labeled ‘islet’) compared to other regions. We
can then plot the set of proteins.

``` r
islet_res <- calc_spatial_diff_ex(pooled.panc.spe,
    assay_name = "proteomics",
    log_transformed = FALSE,
    category_col = "IsletOrNot"
)

# we filter the significant proteins first
sig_prots <- subset(rowData(islet_res), 
                    NonIslet_vs_Islet.adj.P.Val.limma < 0.01)
# then separate into up-regulated and down-regulated based on fold chnage
ups <- subset(sig_prots, NonIslet_vs_Islet.logFC.limma > 0)
downs <- subset(sig_prots, NonIslet_vs_Islet.logFC.limma < 0)

print(paste(
  "We found", nrow(sig_prots), "significantly differentally \
  expressed proteins including",
  nrow(ups), "upregulated proteins and", nrow(downs), "downregulated"
))
```

    ## [1] "We found 241 significantly differentally \n  expressed proteins including 168 upregulated proteins and 73 downregulated"

Now we can plot those differentially expressed proteins across images.

### Plot differentially expressed proteins

If we are interested in the combined expression of proteins we can also
visualize those.

``` r
spe.plot <- img.spes[[2]]

hup <- spatial_heatmap(spe.plot,
    feature = rownames(downs),
    sample_id = "Image_1",
    image_id = "with_grid",
    label_column = "IsletOrNot",
    interactive = FALSE
)

hup
```

![](spatProt_files/figure-html/plot%20pathway-1.png)

### Pathway enrichment measurements

Now we can calculate the enriched pathways in the islets.

``` r
library(leapR)
data("krbpaths")
ora.res <- enrich_ora(islet_res, geneset = krbpaths, 
                      geneset_name = "krbpaths", 
                      feature_column = "PrimaryGeneName")
print(ora.res[grep("INSULIN", rownames(ora.res)), 
              c("ingroup_n", "pvalue", "BH_pvalue")])
```

    ##                                                                                                           ingroup_n
    ## KEGG_INSULIN_SIGNALING_PATHWAY                                                                                    5
    ## BIOCARTA_INSULIN_PATHWAY                                                                                          2
    ## REACTOME_GLUCOSE_REGULATION_OF_INSULIN_SECRETION                                                                 13
    ## REACTOME_INSULIN_SYNTHESIS_AND_SECRETION                                                                         24
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS         0
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION                                                                         16
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_ACETYLCHOLINE                                                         7
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_GLUCAGON_LIKE_PEPTIDE_1                                              11
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_FREE_FATTY_ACIDS                                                      7
    ## REACTOME_INHIBITION_OF_INSULIN_SECRETION_BY_ADRENALINE_NORADRENALINE                                              4
    ##                                                                                                                 pvalue
    ## KEGG_INSULIN_SIGNALING_PATHWAY                                                                            1.000000e+00
    ## BIOCARTA_INSULIN_PATHWAY                                                                                  2.247236e-01
    ## REACTOME_GLUCOSE_REGULATION_OF_INSULIN_SECRETION                                                          5.549704e-02
    ## REACTOME_INSULIN_SYNTHESIS_AND_SECRETION                                                                  2.522014e-07
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS 1.000000e+00
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION                                                                  3.897737e-02
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_ACETYLCHOLINE                                                 1.077909e-04
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_GLUCAGON_LIKE_PEPTIDE_1                                       9.083131e-06
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_FREE_FATTY_ACIDS                                              2.897170e-05
    ## REACTOME_INHIBITION_OF_INSULIN_SECRETION_BY_ADRENALINE_NORADRENALINE                                      1.628147e-02
    ##                                                                                                              BH_pvalue
    ## KEGG_INSULIN_SIGNALING_PATHWAY                                                                            1.0000000000
    ## BIOCARTA_INSULIN_PATHWAY                                                                                  0.9987657450
    ## REACTOME_GLUCOSE_REGULATION_OF_INSULIN_SECRETION                                                          0.7062802680
    ## REACTOME_INSULIN_SYNTHESIS_AND_SECRETION                                                                  0.0002100838
    ## REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS 1.0000000000
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION                                                                  0.5797883996
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_ACETYLCHOLINE                                                 0.0128271210
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_GLUCAGON_LIKE_PEPTIDE_1                                       0.0037831242
    ## REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_FREE_FATTY_ACIDS                                              0.0070612576
    ## REACTOME_INHIBITION_OF_INSULIN_SECRETION_BY_ADRENALINE_NORADRENALINE                                      0.3686245089

### Pathway plotting

We know that there are significantly enriched pathways in insulin
secretion, so let’s plot those.

``` r
secprots <- ora.res["REACTOME_GLUCOSE_REGULATION_OF_INSULIN_SECRETION", ] |>
    dplyr::select(ingroupnames) |>
    unlist() |>
    strsplit(split = ", ") |>
    unlist()

spe.plot <- img.spes[[2]]

hup <- spatial_heatmap(spe.plot,
    feature = secprots,
    sample_id = "Image_1",
    image_id = "with_grid",
    feature_type = "PrimaryGeneName",
    label_column = "IsletOrNot",
    plot_title = "Glucose regulation proteins",
    interactive = FALSE
)

hup
```

![](spatProt_files/figure-html/pathway%20plotting-1.png)

The average expression of the 20 proteins selected is shown to be higher
in islet cells than adjacent cells.

## Distance based measurements

We can also identify features that are correlated with distance to a
feature or a gradient in the sample. This will provide input to
rank-based statistical tools that can help identify pathways.

### Distance based measurements

First we identify a specific feature, the Islet cell, and use that to
identify proteins correlated with distance from the islet in each image.
Proteins with a negative correlation are decreasing in expression as
they are farther from the islet cells.

``` r
## for each image, let's compute the distance of each voxel to 
## the one labeled 'Islet'
rank.imgs <- lapply(
  img.spes,
  function(x) {
    distance_based_analysis(x, "proteomics",
      sampleCategoryCol = "IsletOrNot",
      sampleCategoryValue = "Islet"
    )
  }
)

## now we have the distances, let's plot some interesting proteins
negProts <- do.call(rbind, lapply(names(rank.imgs), function(x) {
    subset(as.data.frame(rowData(rank.imgs[[x]])), 
           IsletDistancespearmanPval < 0.01) |>
      subset(IsletDistancespearmanCor < (-.75)) |>
      dplyr::select(PrimaryGeneName, IsletDistancespearmanCor) |>
      dplyr::mutate(image = x)
}))

print(head(negProts))
```

    ##                       PrimaryGeneName IsletDistancespearmanCor   image
    ## sp|A6NFH5|FBP12_HUMAN          FABP12               -0.9621024 Image_0
    ## sp|O14979|HNRDL_HUMAN         HNRNPDL               -0.9114654 Image_0
    ## sp|O15145|ARPC3_HUMAN           ARPC3               -0.8861469 Image_0
    ## sp|O43292|GPAA1_HUMAN           GPAA1               -0.8355100 Image_0
    ## sp|O43829|ZBT14_HUMAN          ZBTB14               -0.8608285 Image_0
    ## sp|O60256|KPRB_HUMAN          PRPSAP2               -0.8355100 Image_0

``` r
## do any proteins show up more than once?
icounts <- negProts |>
    dplyr::group_by(PrimaryGeneName) |>
    dplyr::summarize(numImgs = dplyr::n()) |>
    dplyr::arrange(desc(numImgs))


print(icounts)
```

    ## # A tibble: 267 × 2
    ##    PrimaryGeneName numImgs
    ##    <chr>             <int>
    ##  1 CHGA                  2
    ##  2 DCAKD                 2
    ##  3 ECHS1                 2
    ##  4 PFKL                  2
    ##  5 PKM                   2
    ##  6 PSMD4                 2
    ##  7 SH3GL1                2
    ##  8 VAMP2                 2
    ##  9 AAK1                  1
    ## 10 ABCB7                 1
    ## # ℹ 257 more rows

It looks like SH3GL1 is correlated with distance to Islet in a few
images.

### Plot protein gradient

Now we can plot the expression of a protein suspected to have decreasing
expression farther from the islet cells.We start with SH3GL1 and SP3S2.

``` r
spatial_heatmap(img.spes[[3]],
    feature = "SH3GL1",
    feature_type = "PrimaryGeneName",
    sample_id = names(img.spes)[3],
    image_id = "with_grid",
    label_column = "IsletOrNot", interactive = FALSE
)
```

![](spatProt_files/figure-html/plot%20correlated%20proteins-1.png)

``` r
spatial_heatmap(img.spes[[2]],
    feature = "AP3S2",
    feature_type = "PrimaryGeneName",
    sample_id = names(img.spes)[2],
    image_id = "with_grid",
    label_column = "IsletOrNot", interactive = FALSE
)
```

![](spatProt_files/figure-html/plot%20correlated%20proteins-2.png)

The expression of this protein is lower farther from the Islet. Can we
identify trends in the proteins?

### Gradient-based enrichment

Rank-based pathway enrichment is a way to evaluate trends pathways that
are over-represented in a ranked list of genes. The `leapR` pathway has
such functionality and we can use the rankings as input.

``` r
library(leapR)
data("krbpaths")

spe <- rank.imgs[[2]]
enriched.paths <- enrich_gradient(spe,
    geneset = krbpaths,
    feature_column = "PrimaryGeneName", # mapped to enrichment data
    ranking_column = "IsletDistancespearmanCor"
)
enriched.paths[, "comp"] <- rep(names(rank.imgs)[[2]], nrow(enriched.paths))
enriched.paths[, "krbpaths"] <- rownames(enriched.paths)

enriched.paths |>
  subset(pvalue < 0.01) |>
  dplyr::group_by(krbpaths) |>
  dplyr::summarize(numImgs = dplyr::n()) |>
  dplyr::arrange(desc(numImgs))
```

    ## # A tibble: 2 × 2
    ##   krbpaths                                 numImgs
    ##   <chr>                                      <int>
    ## 1 REACTOME_DIABETES_PATHWAYS                     1
    ## 2 REACTOME_INSULIN_SYNTHESIS_AND_SECRETION       1

We can see that numerous pathways are coming up as enriched across
images, including ribosomal and translation related pathways. Now we can
select proteins from a particular pathway and visualize those as well.

### Plotting pathways from gradient

``` r
rprots <- subset(enriched.paths, krbpaths == "REACTOME_INSULIN_SYNTHESIS_AND_SECRETION") |>
    dplyr::select(comp, ingroupnames)

rprots <- unlist(strsplit(rprots[1, 2], split = ", "))

spatial_heatmap(img.spes[[2]],
    feature = rprots,
    feature_type = "PrimaryGeneName",
    sample_id = names(img.spes)[2],
    image_id = "with_grid",
    label_column = "IsletOrNot", interactive = FALSE
)
```

![](spatProt_files/figure-html/gradient%20plotting-1.png)

This shows the ribosomal protein expression across the image.

### Network plotting

We can also look at the correlation of the ribosomal proteins in a
graph. The correlation code takes a while but then we can reduce the
graph to the proteins we are most interested in, or those that are most
correlated.

``` r
library(tidygraph)
```

    ## 
    ## Attaching package: 'tidygraph'

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     active, slice

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     active, rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(ggraph)
```

    ## Loading required package: ggplot2

``` r
##correlation analysis can be slow, so let's only evaluate the top 1000 most variable proteins
varprots = apply(assay(img.spes[[2]]),1,var,na.rm = TRUE) |>
  sort(decreasing = TRUE) |>
  names()

date()
```

    ## [1] "Wed Feb  4 15:21:54 2026"

``` r
full_graph <- spatial_network(img.spes[[2]][varprots[1:1000],],
                              'proteomics','PrimaryGeneName')
```

    ## Joining with `by = join_by(rowval)`

``` r
date()
```

    ## [1] "Wed Feb  4 15:22:03 2026"

``` r
##how subset for only those 81 proteins
rgraph <- full_graph |>
  tidygraph::activate(nodes) |>
  dplyr::filter(name %in% rprots) |>#[sample(20)]) |>
  tidygraph::activate(edges) |>
  dplyr::filter(corval > 0.75)

##then we can plot
ggraph::ggraph(rgraph) + 
   geom_edge_link(aes(colour = corval)) + 
   geom_node_point() + 
   geom_node_label(aes(label = name))
```

    ## Using "stress" as default layout

![](spatProt_files/figure-html/graph%20correlations-1.png) Here are the
highly correlated edges between 20 randomly sampled ribosomal proteins.

## Summary

This vignette shows various functions to apply in managing spatial
proteomics data in spammR.

## References

1.  Gosline et al.
2.  Spatial Experiment

## Session info

    ## R version 4.5.2 (2025-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.2
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggraph_2.2.2                ggplot2_4.0.1              
    ##  [3] tidygraph_1.3.1             leapR_0.99.6               
    ##  [5] spammR_0.99.17              limma_3.66.0               
    ##  [7] SpatialExperiment_1.20.0    SingleCellExperiment_1.32.0
    ##  [9] SummarizedExperiment_1.40.0 Biobase_2.70.0             
    ## [11] GenomicRanges_1.62.1        Seqinfo_1.0.0              
    ## [13] IRanges_2.44.0              S4Vectors_0.48.0           
    ## [15] BiocGenerics_0.56.0         generics_0.1.4             
    ## [17] MatrixGenerics_1.22.0       matrixStats_1.5.0          
    ## [19] BiocStyle_2.38.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] DBI_1.2.3           deldir_2.0-4        gridExtra_2.3      
    ##   [4] s2_1.1.9            rlang_1.1.7         magrittr_2.0.4     
    ##   [7] otel_0.2.0          e1071_1.7-17        compiler_4.5.2     
    ##  [10] png_0.1-8           systemfonts_1.3.1   vctrs_0.7.1        
    ##  [13] wk_0.9.5            pkgconfig_2.0.3     fastmap_1.2.0      
    ##  [16] backports_1.5.0     magick_2.9.0        XVector_0.50.0     
    ##  [19] labeling_0.4.3      utf8_1.2.6          rmarkdown_2.30     
    ##  [22] tzdb_0.5.0          ragg_1.5.0          purrr_1.2.1        
    ##  [25] xfun_0.56           cachem_1.1.0        jsonlite_2.0.0     
    ##  [28] DelayedArray_0.36.0 tweenr_2.0.3        broom_1.0.12       
    ##  [31] R6_2.6.1            bslib_0.10.0        RColorBrewer_1.1-3 
    ##  [34] reticulate_1.44.1   boot_1.3-32         car_3.1-3          
    ##  [37] jquerylib_0.1.4     Rcpp_1.1.1          bookdown_0.46      
    ##  [40] knitr_1.51          readr_2.1.6         Matrix_1.7-4       
    ##  [43] igraph_2.2.1        tidyselect_1.2.1    rstudioapi_0.18.0  
    ##  [46] abind_1.4-8         yaml_2.3.12         viridis_0.6.5      
    ##  [49] lattice_0.22-7      tibble_3.3.1        withr_3.0.2        
    ##  [52] S7_0.2.1            evaluate_1.0.5      sf_1.0-24          
    ##  [55] desc_1.4.3          units_1.0-0         proxy_0.4-29       
    ##  [58] spData_2.3.4        polyclip_1.10-7     pillar_1.11.1      
    ##  [61] BiocManager_1.30.27 ggpubr_0.6.2        carData_3.0-5      
    ##  [64] KernSmooth_2.23-26  plotly_4.12.0       sp_2.2-0           
    ##  [67] hms_1.1.4           scales_1.4.0        class_7.3-23       
    ##  [70] glue_1.8.0          lazyeval_0.2.2      tools_4.5.2        
    ##  [73] data.table_1.18.2.1 ggnewscale_0.5.2    ggsignif_0.6.4     
    ##  [76] fs_1.6.6            graphlayouts_1.2.2  grid_4.5.2         
    ##  [79] spdep_1.4-1         impute_1.84.0       tidyr_1.3.2        
    ##  [82] ggforce_0.5.0       Formula_1.2-5       cli_3.6.5          
    ##  [85] textshaping_1.0.4   S4Arrays_1.10.1     viridisLite_0.4.2  
    ##  [88] dplyr_1.1.4         gtable_0.3.6        rstatix_0.7.3      
    ##  [91] sass_0.4.10         digest_0.6.39       classInt_0.4-11    
    ##  [94] SparseArray_1.10.8  ggrepel_0.9.6       rjson_0.2.23       
    ##  [97] htmlwidgets_1.6.4   farver_2.1.2        memoise_2.0.1      
    ## [100] htmltools_0.5.9     pkgdown_2.2.0       lifecycle_1.0.5    
    ## [103] httr_1.4.7          statmod_1.5.1       MASS_7.3-65
