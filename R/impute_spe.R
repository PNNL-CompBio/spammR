#' Impute missing values based on spatial coordinates
#' @description `impute_spe()` carries out imputation for missing data in a data frame (df) using a specified method from a range of methods. Accepts input dataset as a data frame (df).
#' @importFrom matrixStats rowMedians
#' @import impute
#' @import spdep
#' @export
#' @param spe SPE containing data Data frame containing data to be imputed, where rows correspond to features (which can be specified as row names of dat but it is not required that they be specified) and columns correspond to samples. Column names must correspond to names provided in the sample identifier column in the metadata parameter.
#' @param method Method of imputation to be used. See details.
#' @param spatial_unit_colname Column name in metadata that specifies spatial unit information. Example: ROI_abbreviation.
#' @param knn_k K value to be used for k-nearest neighbor imputation
#' @param spatialCoord_x_colname Column name in metadata that specifies the X coordinate for each sample.
#' @param spatialCoord_y_colname Column name in metadata that specifies the Y coordinate for each sample.
#' @param allowed_missingness_perProtein Proportion of samples allowed to have missing data for a protein in the given spatial unit specified by the imputation method. Example: When method="global_mean," an allowed_missingness_perProtein of 0.5 indicates that any protein missing data for more than 50% of samples across the entire spatial tissue covered by all samples will be excluded from the imputation method algorithm, and that protein's missing values will not be imputed. When method="mean_perSpatialUnit," then allowed_missingness_perProtein of 0.5 indicates that a protein must have data for at least 50% of samples in the specified spatial unit (example: a brain ROI) to be used in the imputation algorithm and to be imputed. This argument is only needed for methods that rely on statistics based on per protein. It is used for the following methods: median_a, median_b, global_mean, mean_perSpatialUnit, knn_proteins_global and knn_proteins_perSpatialUnit.
#' @param allowed_missingness_perSample  Proportion of proteins allowed to be missing for a sample. This argument is only needed when imputing using k-nearest neighbor methods where neighbors are samples. For all other imputation methods provided here, allowed_missingness_perSample does not apply; no samples are excluded from the imputation process for those methods. If missingness in a sample is greater than this thershold, then that sample's data are not included in the imputation for the knn samples neighbors methods.
#' @returns an SPE with an 'imputed' data frame

#' @details Methods options and descriptions
#' zero: replace missing values with 0
#' median_a : replace missing values with global median per protein
#' median_b : replace missing values with 1/2 global median per protein
#' global_mean: replace missing values with global mean per protein
#' mean_perSpatialUnit: replace msising values with mean per spatial unit (example: ROI) for each protein
#' knn_proteins_global: imputation based on k-nearest neighbors, with proteins as neighbors, based on data from all samples across all ROIs
#' knn_proteins_perSpatialUnit: imputation based on k-nearest neighbors, with proteins as neighbors, based on data from specified spatial unit (here: ROI)
#' knn_samples_global_proteinData: imputation based on k-nearest neighbors, with samples as neighbors, based on protein data from all samples.
#' knn_samples_global_spatialCoords: imputation based on knn, with samples as neighbors, based on spatial coordinates of all samples.

impute_spe <- function(spe,
                      assay_name,
                      method = c('zero','median','median_half','global_mean','spatial_unit_mean','knn_global','spatial_unit_knn','knn_global','knn_splatial'),
                      spatial_unit_colname,
                      spatialCoord_x_colname,
                      spatialCoord_y_colname,
                      knn_k=NULL,
                      allowed_missingness_perProtein=NULL,
                      allowed_missingness_perSample=NULL){

  # # Set defult values if not specified by user
  # if (is.null(allowed_missingness_perProtein)){ # Only used for knn methods
  #   allowed_missingness_perProtein = 0.75
  # }
  # if (is.null(allowed_missingness_perSample)){ # Only used for knn methods
  #   allowed_missingness_perSample = 0.80
  # }
  dat <-SummarizedExperiment::assays(spe)[[assay_name]]
  metadata<-SummarizedExperiment::colData(spe)

  data_imputed<-impute_df(dat,method,metadata,spatial_unit_colname,spatialCoord_x_colname,
                         spatialCoord_y_colname,knn_k,  allowed_missingness_perProtein,
                          allowed_missingness_perSample)

  assays(spe)<-list(assays(spe),imputed=data_imputed)
  return (spe)
}


impute_df<-function(dat,
                    method = c('zero','median','median_half','global_mean','spatial_unit_mean','knn_global','spatial_unit_knn','knn_global','knn_splatial'),
                    metadata,
                    spatial_unit_colname,
                    spatialCoord_x_colname,
                    spatialCoord_y_colname,
                    knn_k=NULL,
                    allowed_missingness_perProtein=NULL,
                    allowed_missingness_perSample=NULL){
  data_imputed = c()
  replace_vals = c()
  if (method=="zero"){
    replace_vals=0
  }else if (method=="global_mean"){
    replace_vals = rowMeans(dat,na.rm=TRUE)
  }else if (method=="median"){
    replace_vals = matrixStats::rowMedians(as.matrix(dat),na.rm=TRUE)
  }else if (method == "median_half"){
    replace_vals = 0.5* matrixStats::rowMedians(as.matrix(dat),na.rm=TRUE)
  }else if (method == "spatial_unit_mean" | method=="spatial_unit_knn"){
    # extract column with spatial unit specification in the metadata
    spatUnits = unique(SummarizedExperiment::colData(spe)[,spatial_unit_colname]) ##whats is the spatial unit?
    imputed_data = c()
    for (s in 1:length(spatUnits)){
      ROI_indices = rownames(SummarizedExperiment::colData(spe))[which(SummarizedExperiment::colData(spe)[,spatial_unit_colname]==s)]
      #        grep(spatUnits[s],colnames(dat)) ##FIX: spatial unit should not be in column name
      ROI_data = dat[,ROI_indices]
      if (method == "spatial_unit_mean"){
        imputed_ROI_data = impute_df(ROI_data,method="global_mean",allowed_missingness_perProtein=allowed_missingness_perProtein)
      }else if (method == "spatial_unit_knn"){
        imputed_ROI_data = impute_df(ROI_data,method="knn_global",knn_k=knn_k,allowed_missingness_perProtein=allowed_missingness_perProtein)
      }
      if (s==1){
        imputed_data = imputed_ROI_data
      }else{
        imputed_data = cbind(imputed_data,imputed_ROI_data)
      }
    }
    # Make sure the imputed data appears in the same order as the original data
    imputed_data = imputed_data[,colnames(dat)]
    data_imputed = imputed_data
  }else if (method == "knn_proteins_global"){
    knn_imputed_proteins_global = impute::impute.knn(as.matrix(dat),k = knn_k, rowmax = allowed_missingness_perProtein, colmax = 1, maxp = dim(dat)[1], rng.seed=12345)
    # Fix immputed values for rows (proteins) that are missing more than rowmax proportion of samples.
    # By default, impute.knn uses the mean of all proteins in a sample for this case, which is not that meaningful.
    # Instead, we replace the missing values (for the case when missingness in a row is greater than rowmax) with global_mean or na for that row
    fix_prots = which( (rowSums(is.na(dat))/dim(dat)[2]) > allowed_missingness_perProtein )
    print(fix_prots)
    for (f in fix_prots){
      na_indices = which(is.na(dat[f,]))
      #knn_imputed_proteins_global$dat[f,na_indices] = mean(dat[f,],na.rm=TRUE)
      knn_imputed_proteins_global$dat[f,na_indices] = NA #Keep those missing values as NA since there is not enough dat to impute reliably.
    }
    if (!is.null(fix_prots)){
      print("impute.knn function uses mean imputation for rows with more than the specified missingness. spammR's imputation function reverts these missing values back to NA.")
    }
    data_imputed = knn_imputed_proteins_global$dat
  }else if (method == "knn_samples_global_proteinData"){
    transposed_data = as.matrix(t(dat))
    knn_imputed_samples_global_protData = impute::impute.knn(transposed_data,k = knn_k, rowmax = allowed_missingness_perSample, colmax = 1, maxp = dim(transposed_data)[1], rng.seed=12345)
    imputed_temp =  t(knn_imputed_samples_global_protData$dat)
    fix_samples = which( (colSums(is.na(dat))/dim(dat)[1]) > allowed_missingness_perSample )
    for (f in fix_samples){
      na_indices = which(is.na(dat[,f]))
      #knn_imputed_proteins_global$data[f,na_indices] = mean(data[f,],na.rm=TRUE)
      imputed_temp[na_indices,f] = NA #Keep those missing values as NA since there is not enough data to impute reliably.
    }
    fix_prots = which( (rowSums(is.na(dat))/dim(dat)[2]) > allowed_missingness_perProtein )
    for (f in fix_prots){
      na_indices = which(is.na(dat[f,]))
      #knn_imputed_proteins_global$data[f,na_indices] = mean(data[f,],na.rm=TRUE)
      imputed_temp[f,na_indices] = NA #Keep those missing values as NA since there is not enough data to impute reliably.
    }
    # if (!is.null(fix_samples)){
    #   print("impute.knn function uses mean imputation for rows with more than the specified missingness. spammR's imputation function reverts these missing values back to NA.")
    # }
    data_imputed = imputed_temp
  }else if (method == "knn_samples_global_spatialCoords"){
    d_for_nearestNeighb = as.matrix(metadata[,c("Xcoord","Ycoord")])
    # Find k-nearest neighbors based on spatial distance (x,y coordinates)
    knn <- spdep::knearneigh(d_for_nearestNeighb, k= knn_k)
    # nearest neighbor indices are in knn$nn
    nn = knn$nn
    for (prot in 1:dim(dat)[1]){
      this.prot = dat[prot,]
      na_indices = which(is.na(this.prot))
      for (n in na_indices){
        n.knn = nn[n,] #nearest neighbors for sample corresponding to nth index
        this.prot[n] = mean(as.numeric(this.prot[n.knn]),na.rm=TRUE)
      }
      data_imputed = rbind(data_imputed,this.prot)
    }
  }
  if (is.null(replace_vals)){
    #don't change data_imputed
  }else{
    data_imputed = c() #start with this, then replace na values
    # Update replace_vals to reflect allowed_missingness_perProtein criteria.
    # replace_vals should retain NAs (missing values) for proteins above the allowed_missingness_perProtein criteria for the following methods
    dontImpute_rows = c()
    if (method=="median_a" | method=="median_b" | method=="global_mean"){
      dontImpute_rows = which ( (rowSums(is.na(dat))/dim(dat)[2]) > allowed_missingness_perProtein)
    }
    replace_vals[dontImpute_rows] = NA
    for (i in 1:dim(dat)[1]){
      this.row = dat[i,]
      this.row.na_indices = which(is.na(this.row))
      if (length(this.row.na_indices)>0){
        if (length(replace_vals)>1){
          this.row[this.row.na_indices] <- replace_vals[i]
        }else{
          this.row[this.row.na_indices] <- replace_vals
        }
      }
      data_imputed = rbind(data_imputed,this.row)
    }
  }
  data_imputed = data.frame(data_imputed, row.names=rownames(dat))
  colnames(data_imputed) = colnames(dat)
}
