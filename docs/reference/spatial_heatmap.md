# Plot data in heatmap together with image

`spatial_heatmap()` Creates a spatial heatmap for a given feature in a
SpatialExperiment object (spe) and provides the option to label each
sample in the x-y plane

## Usage

``` r
spatial_heatmap(
  spe,
  feature,
  feature_type = NA,
  assay_name = "proteomics",
  plotBackground_img = TRUE,
  sample_id,
  image_id,
  spatial_coord_type = "",
  spot_size_name = c("spot_width", "spot_height"),
  metric_display = "Protein expression",
  plot_log = FALSE,
  label_column = NA,
  fill_colors = c("slateblue", "darkgreen", "goldenrod"),
  sample_label_color = "white",
  sample_label_size = 1.75,
  plot_title = NULL,
  title_size = 10,
  show_na = FALSE,
  interactive = FALSE
)
```

## Arguments

- spe:

  SpatialExperiment (SPE) object

- feature:

  Name of the feature in the spe object whose values are to plotted in
  the spatial heat map. This should be a row name in rowData(spe)

- feature_type:

  Example: "GeneName". Default is whatever identifier is used in
  rownames. The name of feature_type must be present as a column in
  rowData(spe)

- assay_name:

  Name of assay in the spe object that contains data to be plotted

- plotBackground_img:

  Boolean (TRUE or FALSE) to indicate whether a background image should
  be plotted. Default is FALSE. If TRUE, the parameters image_boundaries
  and image_sample_id must be specified, and the image data must be
  present in the spe object, under imgData(spe)

- sample_id:

  The names of the background image's sample_id fields in the spe
  object; this provides a unique identifier for the background image to
  be used for plotting (if there are multiple images under imgData(spe))
  and only plots samples associated with the specified image sample_id.
  Example: c("Image0","Raw_noMarkings"). Image data stored under the spe
  object can be viewed by imgData(spe)

- image_id:

  The name of the background image image_id. Together with the
  `sample_id` this provides a unique imge to plot.

- spatial_coord_type:

  Position type for the given spatial coordinates of samples in spe.
  Current options are: "topleft_corner", "topright_corner". Default is
  blank, which assumes coordinates are bottom right.

- spot_size_name:

  Is a vector of length 2 describing the columns to find width and
  height of each spot.

- metric_display:

  Legend title for spatial heatmap. If this parameter is not specified,
  legend title defaults to "Protein expression"

- plot_log:

  Boolean to indicate if data should be log-transformed before plotting.
  Defaults to FALSE

- label_column:

  Column in colData(spe) to be used for labeling grid squares. If not
  specified, default is no labels.

- fill_colors:

  List of colors to use to provide gradient for heatmap plotting

- sample_label_color:

  Color to be used for labels of samples/grid squares. Default is white.

- sample_label_size:

  Font size for sample labels. Default is 1.75.

- plot_title:

  Title to be given to the spatial heatmap. Default is "Spatial
  signature for XYZ" where XYZ is the name of the specified feature

- title_size:

  Font size of title

- show_na:

  Set to TRUE if you want to plot NA values (will be grey) otherwise
  they will not be plotted

- interactive:

  Boolean value (TRUE/FALSE) indicating whether the plot should have
  interactive mouse hovering. If not specified, this defaults to TRUE.
  Note: grid squares can only be labeled when interactive = FALSE due to
  current ggplotly limitations.

## Value

A spatial heatmap of the chosen feature

## Details

Assumes that every sample (column) in assay(spe) has a corresponding x,y
coordinate given in spatialCoords(spe)

Assumes that colData(spe) and spatialCoords(spe) are provided in the
same order of samples (rows).

Assumes that the spe object contains background image data, stored in
imgData(spe) if plotBackground_img is specified to be TRUE below.

Default values are specified for function parameters 'metric_display',
'label_column' and 'interactive,' if the user doesn't specify those.

Future versions of spatialPlot_feature() to include additional options
for spatial_coord_type: "center", bottomleft_corner",
"bottomright_corner"

## Examples

``` r
data(pancMeta)
data(smallPancData)
data(protMeta)
img0.spe <- convert_to_spe(smallPancData$Image_0,
  pancMeta,
  protMeta,
  feature_meta_colname = "pancProts",
  image_files = system.file("extdata", "Image_0.png", package = "spammR"),
  spatial_coords_colnames = c("x_pixels", "y_pixels"),
  sample_id = "Image0",
  image_ids = "with_grid"
)
#> Note: Only mapping metadata for 2986 features out of 3000 data points
res = spatial_heatmap(img0.spe,
        feature='INS',
        sample_id='Image0',
        image_id='with_grid',
        feature_type='PrimaryGeneName',
        label_column='IsletOrNot', interactive=FALSE)
```
