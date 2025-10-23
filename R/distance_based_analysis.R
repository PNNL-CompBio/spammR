#' Identify multiomic features correlated with an distance to an image feature
#' @description distance_based_analysis: Identifies proteins/features that show
#' a strong correlation between distance from a specified ROI's samples and
#' protein/feature abundance differences between samples.
#' @import SpatialExperiment
#' @importFrom IRanges IRanges
#' @importFrom IRanges countOverlaps
#' @importFrom stats cor.test
#' @export
#' @param spe SpatialExperiment object containing spatial omics data
#' @param assay_name Name of the assay stored in spe that is to be used for 
#' distance based analysis. Example: "znormalized_log2"
#' @param spotHeightCol Column containing height of spot
#' @param spotWidthCol Column containing width of spot
#' @param sampleCategoryCol Column name in metadata (colData(spe)) that should 
#' be used for selecting samples of certain type, to define the "origin" region 
#' for distance based analysis
#' @param sampleCategoryValue Sample category to be used for defining the 
#' "origin" region for distance based analysis
#' @param featuresNameCol Name of column containing features (example: proteins)
#'  in rowData(spe). It is assumed that the data provided in assay(spe) is in 
#'  the same order as the order in which the features are listed under the 
#'  featuresNameCol
#' @param corr_type Choose from "pearson" (default), "spearman." Correlation 
#' method to be used for calculating correlation between distance between 
#' samples and protein abundance differences. Both types of correlation provide
#'  a measure of monotonic association between two variables. Pearson is better
#'   suited for linear relationships while Spearman is better for non-linear 
#'   monotonic relationships.
#' @param corr_thresh Minimum correlation value to be used for identifying 
#' proteins that have a correlation between protein abundance differences and 
#' distance between samples. Values greater than or equal to this theshold will
#'  be used.
#' @param min_samples of a minimum number of sample points for calculating 
#' correlation. For proteins with less than this number of sample points, 
#' correlation value is reported as NA.
#' @param allowOverlaps allow overlaps of regions
#' @returns a `SpatialExperiment` object with the distance based data in the
#' colData
#' @examples
#' data(pancMeta)
#' data(protMeta)
#' data(smallPancData)
#' img0.spe <- convert_to_spe(smallPancData$Image_0,
#'   pancMeta,
#'   protMeta,
#'   sample_id = "Image0",
#'   spatial_coords_colnames = c("x_pixels", "y_pixels"),
#'   feature_meta_colname = "pancProts",
#'   image_files = system.file("extdata", "Image_0.png", package = "spammR"),
#'   image_ids = "Image0"
#' )
#'
#' img0.spe <- distance_based_analysis(img0.spe,
#'   "proteomics",
#'   sampleCategoryCol = "IsletOrNot",
#'   sampleCategoryValue = "Islet"
#' )
#'
distance_based_analysis <- function(spe,
                                    assay_name,
                                    spotHeightCol = "spot_height",
                                    spotWidthCol = "spot_width",
                                    sampleCategoryCol,
                                    sampleCategoryValue,
                                    featuresNameCol,
                                    corr_type = "spearman",
                                    corr_thresh = 0.5,
                                    min_samples = 5,
                                    allowOverlaps = TRUE) {
  # Compute centroids for each sample based on top-left corner (Xcoord, Ycoord)
   # coordinates (SG: should be bottom!?!)
      spatial_coords <- data.frame(spatialCoords(spe))
    
      ## first we use Iranges to check to see if there are overlaps (in x-space)
      ## SG: Maybe remove this...
      if (!allowOverlaps) {
          xranges <- IRanges::IRanges(start = spatial_coords[, 1], 
                                      width = colData(spe)[[spotWidthCol]]) |>
            unique()
          yranges <- IRanges::IRanges(start = spatial_coords[, 2], 
                                      width = colData(spe)[[spotHeightCol]]) |>
            unique()
          allovers <- c(IRanges::countOverlaps(xranges), 
                        IRanges::countOverlaps(yranges))
          if (any(allovers) > 1) {
            stop("overlapping regions found, set allowOverlaps to TRUE")
          }
      }
    
      ## get centers of points
      centroid_x <- as.numeric(spatial_coords[, 1] + 
                               colData(spe)[[spotWidthCol]] / 2) #
      centroid_y <- as.numeric(spatial_coords[, 2] + 
                                 colData(spe)[[spotHeightCol]] / 2) 
      centroid_coords <- data.frame(cbind(centroid_x, centroid_y))
      rownames(centroid_coords) <- rownames(spatial_coords)
      
      # Compute distance between samples
      dist_between_samples <- as.matrix(stats::dist(centroid_coords,
        method = "euclidean",
        diag = TRUE
      ))
    
      assay_data <- SummarizedExperiment::assays(spe)[[assay_name]]
      source_samples_indices <- which(SummarizedExperiment::colData(spe)[,
                                  sampleCategoryCol] == sampleCategoryValue)
      if (length(source_samples_indices) > 1) {
          message("ERROR: More than one feature to measure distance from")
          exit()
      }
    
      ## get the distances
      samp_dists <- dist_between_samples[source_samples_indices, ]
    
      ## correlate the distances
      cor_dist <- cor(t(assay(spe)), samp_dists, 
                      method = corr_type, use = "pairwise.complete.obs")
      cor_test <- apply(assay(spe), 1, function(x) {
        pval <- 1.0
        try(pval <- stats::cor.test(x, samp_dists, 
                                    method = corr_type, 
                                    exact = FALSE)$p.val, 
            silent = TRUE)
        return(pval)
      })
    
      mdata <- data.frame(cv = cor_dist, cp = cor_test)
      colnames(mdata) <- c(
        paste0(sampleCategoryValue, "Distance", corr_type, "Cor"),
        paste0(sampleCategoryValue, "Distance", corr_type, "Pval")
      )
    
      rowData(spe) <- cbind(rowData(spe), mdata)
      return(spe)
  }
