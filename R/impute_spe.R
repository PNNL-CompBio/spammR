#' Impute missing values based on spatial coordinates
#' @description `impute_spe()` carries out imputation for missing data in the 
#' primary assay using a specified method from a range of methods. 
#' Returns additional assay in the same SPE
#' @importFrom matrixStats rowMedians
#' @importFrom impute impute.knn
#' @importFrom spdep knearneigh
#'
#' @param spe SPE containing data to be imputed.
#' @param assay_name name of assay with data to be imputed
#' @param method Method of imputation to be used. See details.
#' @param group_colname Column name in metadata that specifies the group 
#' information to use for group_mean or knn_group. Example: ROI_abbreviation.
#' @param k K value to be used for k-nearest neighbor imputation
#' @param protein_missingness Proportion of samples allowed to have missing
#'  data for a protein in the given imputation method. Example: When 
#'  method="global_mean," an protein_missingness of 0.5 indicates that
#' any protein missing data for more than 50% of samples across the entire 
#' spatial tissue covered by all samples will be excluded from the imputation 
#' method algorithm, and that protein's missing values will not be imputed.
#' When method="group_mean," then protein_missingness of 0.5 indicates that 
#' a protein must have data for at least 50% of samples in the specified group 
#' to be used in the imputation algorithm and to be imputed.
#' @returns An SPE with an 'imputed' assay with the appropriate imputation 
#' called
#' @export
#'
#' @details Methods options and descriptions:
#' - zero: replace missing values with 0
#' - median : replace missing values with global median per protein
#' - median_half : replace missing values with 1/2 global median per protein
#' - mean: replace missing values with global mean per protein
#' - group_mean: replace missing values with mean per group, e.g. group 
#' (example: ROI) for each protein
#' - knn: imputation based on k-nearest neighbors, with proteins as neighbors,
#'  based on data from all samples across all groups. NOTE: There will still 
#'  be NA values if the protein is not expressed in this group.
#' - group_knn: imputation based on k-nearest neighbors, with proteins as 
#' neighbors, based on data from specified group (e.g. ROI, tissue). NOTE: 
#' There will still be NA values if the protein is not expressed in this group.
#' - spatial_knn: imputation based on k-nearest neighors in space
#' @examples
#'
#' data(pancMeta)
#' data(protMeta)
#' data(smallPancData)
#' # We can put all samples into the same object (for statistical power)
#' pooledData <- dplyr::bind_cols(smallPancData)
#' pooled.panc.spe <- convert_to_spe(pooledData,
#'   pancMeta,
#'   protMeta,
#'   feature_meta_colname = "pancProts",
#'   sample_id = ""
#' )
#' # we can try two imputation methods and compare the difference
#' res <- impute_spe(pooled.panc.spe, method = "mean")
#' res2 <- impute_spe(pooled.panc.spe, method = "group_mean", 
#'                   group_colname = "Image")
#' mean(assay(res, "imputed") - assay(res2, "imputed"), na.rm = TRUE)
#'
impute_spe <- function(spe,
                       assay_name = NULL,
                       method = NULL,
                       group_colname,
                       k = NULL,
                       protein_missingness = NULL) {
  
    .methods <- c("zero", "median", "median_half", "mean", 
                  "group_mean", "sampMin", "knn", "group_knn", "spatial_knn")
    
  
    ## check for args
    if (!method %in% .methods) {
       msg <- paste("Method", method, "must be one of", paste(.methods, 
                                                             collapse = ","))
        stop(msg)
    }
    ## check for K in knn
    if (method %in% c("knn", "group_knn", "spatial_knn") && is.null(k)) {
        msg <- paste("Method", method, "requires parameter k to be set")
        stop(msg)
    }
  
    if (is.null(assay_name)) {
        dat <- SummarizedExperiment::assay(spe)
    } else {
        dat <- SummarizedExperiment::assay(spe, assay_name)
    }
  
    metadata <- SummarizedExperiment::colData(spe)
    spcoords <- SpatialExperiment::spatialCoords(spe)
  
    imputed_data <- c()
    replace_vals <- c()
  
    #first lets see if there are proteins that we dont impute for global methods
    fix_prots <- which((rowSums(is.na(dat)) / ncol(dat)) > protein_missingness)
  
    ### first iterate through the global methods
    if (method %in% c("zero", "mean", "median", "median_half", "minSamp")) {
        ## row-wise values first
        if (method == "zero") {
            replace_vals <- rep(0, nrow(dat))
        } else if (method == "mean") {
            replace_vals <- rowMeans(dat, na.rm = TRUE)
        } else if (method == "median") {
            replace_vals <- matrixStats::rowMedians(as.matrix(dat), 
                                                    na.rm = TRUE)
        } else if (method == "median_half") {
           replace_vals <- 0.5 * matrixStats::rowMedians(as.matrix(dat), 
                                                        na.rm = TRUE)
        } else if (method == 'minSamp') {
          replace_vals <- matrixStats::colMins(as.matrix(dat),
                                               na.rm = TRUE)
        }
        
        imputed_data <- dat
        ### now do the imputation
        for (i in setdiff(seq_along(1:nrow(dat)), fix_prots)) {
            ## for all of the values we WANT to fix
            imputed_data[i, which(is.na(dat[i, ]))] <- replace_vals[i]
        }
    
      ## if we were mising things
    } else if (method == "knn") {
      kres <- impute::impute.knn(as.matrix(dat), k = k, 
                                 rowmax = protein_missingness, 
                                 colmax = 1, maxp = dim(dat)[1], 
                                 rng.seed = 12345)
      # Fix immputed values for rows (proteins) that are missing more than
      # rowmax proportion of samples.
      # By default, impute.knn uses the mean of all proteins in a sample
      # for this case, which is not that meaningful.
      # Instead, we replace the missing values (for the case when 
      # missingness in a row is greater than rowmax) with global_mean or 
      # na for that row
      imputed_data <- kres$dat
      # now we go back and uninfer values for which we have poor coverage
      for (f in fix_prots) {
        na_indices <- which(is.na(dat[f, ]))
        imputed_data[f, na_indices] <- NA # Keep those missing values as 
        #NA since there is not enough dat to impute reliably.
      }
      if (length(fix_prots) > 0) {
        message("impute.knn function uses mean imputation for rows with \
                more than the specified missingness. spammR's imputation \
                function reverts these missing values back to NA.")
      }
      #      data_imputed = knn_imputed_proteins_global$dat]
    } else if (method == "group_mean" || method == "group_knn") {
      # extract column with spatial unit specification in the metadata
      spatUnits <- unique(metadata[, group_colname]) ## whats the spatial unit?
      imputed_data <- dat
      ## iterate through each group, identify means, and replace those values
      for (s in spatUnits) {
        ROI_indices <- which(metadata[, group_colname] == s)
  
        ## first lets see if there are proteins that we dont impute
        fix_prots <- which((rowSums(is.na(dat[, ROI_indices])) / 
                              length(ROI_indices)) > protein_missingness)
  
        if (method == "group_mean") {
            replace_vals <- rowMeans(dat[, ROI_indices], na.rm = TRUE)
            ## for all of the values we WANT to fix
            for (i in setdiff(seq_along(1:nrow(imputed_data)), fix_prots)) { 
              if (any(is.na(dat[i, ROI_indices]))) {
                imputed_data[i, which(is.na(dat[i, ROI_indices]))] <- 
                  replace_vals[i]
              } ## replace the values WITHIN THIS GROUP
            }
        } else if (method == "group_knn") {
            gdat <- dat[, ROI_indices]
            imputed_gdata <- impute::impute.knn(as.matrix(gdat),
                                                k = k, 
                                                rowmax = protein_missingness, 
                                                colmax = 1, 
                                                maxp = dim(gdat)[1], 
                                                rng.seed = 12345)
            imputed_data[, ROI_indices] <- imputed_gdata$dat
            # now we go back and uninfer values for which we have poor coverage
            for (f in fix_prots) {
              na_indices <- which(is.na(gdat[f, ]))
              if (length(na_indices) > 0) {
                imputed_data[f, na_indices] <- NA
              } # Keep those missing values as NA since there is 
              #not enough dat to impute reliably.
            }
          }
          # Make sure the imputed data appears in the same order as the 
        #original data
          imputed_data <- imputed_data[, colnames(dat)]
      }
    } else if (method == "spatial_knn") {
        if (ncol(spcoords) == 0) {
            stop("Need spatial coordinates for spatial_knn method")
      }
      
        d_for_nearestNeighb <- as.matrix(spcoords)
        # Find k-nearest neighbors based on spatial distance (x,y coordinates)
        knn <- spdep::knearneigh(d_for_nearestNeighb, k = k)
        # nearest neighbor indices are in knn$nn
        nn <- knn$nn
        imputed_data <- vapply(seq_along(1:nrow(dat)), function(prot){
        #for (prot in seq_along(1:dim(dat)[1])) {
          this.prot <- dat[prot, ]
          na_indices <- which(is.na(this.prot))
          for (n in na_indices) {
            n.knn <- nn[n, ] # nearest neighbors for sample corresponding to 
            #nth index
            this.prot[n] <- mean(as.numeric(this.prot[n.knn]), na.rm = TRUE)
          }
          #imputed_data <- rbind(imputed_data, this.prot)
          return(this.prot)
        }, numeric(ncol(dat)))
        imputed_data <- t(imputed_data)
    }
  
    assays(spe, withDimnames = FALSE) <- c(assays(spe), 
                                           list(imputed = imputed_data))
    return(spe)
}
