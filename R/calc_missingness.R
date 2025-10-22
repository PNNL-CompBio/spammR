#' Calculate degree of missingness for each featurea and sample
#' @description calc_missingness() calculates the fraction of samples that are 
#' missing data for each feature, and fraction of each feature that are missing for each sample.
#' @param spe SpatialExperiment object
#' @param group_name Name of group to use as denominator
#' @return SpatialExperiment object with missingFeature column added to colData 
#' and missingSample added to rowData
#' @export
#' @examples
#' library(spammR)
#' data(smallPancData)
#' data(pancMeta)
#' data(protMeta)
#' pooledPanc <- dplyr::bind_cols(smallPancData)
#' panc.spe <- convert_to_spe(pooledPanc, pancMeta, protMeta, 
#'           feature_meta_colname = "pancProts")
#' calc_missingness(panc.spe)
#'
calc_missingness <- function(spe, group_name = NULL) {
    ## row data to use to calculate missing samples
    rd <- as.data.frame(rowData(spe))
  
    ## column data used to calculate missing features
    cd <- as.data.frame(colData(spe))
    mfeat <- apply(assay(spe), 2, function(x) length(which(is.na(x))) / length(x))
    cd <- cbind(cd, missingFeature = mfeat)
  
    if (is.null(group_name)) {
      msamp <- apply(assay(spe), 1, function(x) length(which(is.na(x))) / length(x))
  
      rd <- cbind(rd, missingSample = msamp)
    } else if (group_name %in% names(cd)) {
      groups <- unique(cd[, group_name])
      for (g in groups) {
        samps <- which(cd[, group_name] == g)
        msamp <- apply(assay(spe)[, samps], 1, 
                       function(x) length(which(is.na(x))) / length(x))
        rd[, paste0("missingSample_", g)] <- msamp
      }
    } else {
      stop(paste("group_name", group_name, "not in colData"))
    }
  
    colData(spe) <- DataFrame(cd)
    rowData(spe) <- DataFrame(rd)
    return(spe)
}
