#' impute_df: carries out imputation for missing data in a data frame (df) using a specified method from a range of methods. Accepts input dataset as a data frame (df).
#' @import matrixStats
#' @import impute
#' @export
#' @param dat Data frame containing data to be imputed, where rows correspond to features (which can be specified as row names of dat but it is not required that they be specified) and columns correspond to samples. Column names must correspond to names provided in the sample identifier column in the metadata parameter.
#' @param method Method of imputation to be used. See details.
#' @param metadata A dataframe containing metadata for dat, where rows correspond to samples, and columns contain meta information about samples. This is required to be specified if the imputation method is spatial specific. At the very minimum, columns corresponding to sample identifier, spatial unit (if relevant) and spatial coordinates should be specified. Default is NULL.
#' @param spatial_unit_colname Column name in metadata that specifies spatial unit information. Example: ROI_abbreviation.
#' @param spatialCoord_x_colname Column name in metadata that specifies the X coordinate for each sample.
#' @param spatialCoord_y_colname Column name in metadata that specifies the Y coordinate for each sample.
#' @returns a data frame containing the imputed dataset using the specified method.

#' @details Methods options and descriptions
#' zero: replace missing values with 0\n
#' median_a : replace missing values with global median per protein\n
#' median_b : replace missing values with 1/2 global median per protein\n
#' global_mean: replace missing values with global mean per protein\n
#' mean_perSpatialUnit: replace msising values with mean per spatial unit (example: ROI) for each protein\n
#' knn_proteins_global: imputation based on k-nearest neighbors, with proteins as neighbors, based on data from all samples across all ROIs\n
#' knn_proteins_perSpatialUnit: imputation based on k-nearest neighbors, with proteins as neighbors, based on data from specified spatial unit (here: ROI)\n
#' knn_samples_global_proteinData: imputation based on k-nearest neighbors, with samples as neighbors, based on protein data from all samples.\n
#' knn_samples_global_spatialCoords: imputation based on knn, with samples as neighbors, based on spatial coordinates of all samples.\n

impute_df <- function(dat,method,metadata, spatial_unit_colname, spatialCoord_x_colname, spatialCoord_y_colname, knn_k=NULL,allowed_missingness_perProtein=NULL, allowed_missingness_perSample=NULL){
  library(matrixStats)
  library(impute)
  # # Set defult values if not specified by user
  # if (is.null(allowed_missingness_perProtein)){ # Only used for knn methods
  #   allowed_missingness_perProtein = 0.75
  # }
  # if (is.null(allowed_missingness_perSample)){ # Only used for knn methods
  #   allowed_missingness_perSample = 0.80
  # }
  data_imputed = c()
  replace_vals = c()
  if (method=="zero"){
    replace_vals=0
  }else if (method=="global_mean"){
    replace_vals = rowMeans(dat,na.rm=TRUE)
  }else if (method=="median_a"){
    replace_vals = rowMedians(as.matrix(dat),na.rm=TRUE)
  }else if (method == "median_b"){
    replace_vals = 0.5* rowMedians(as.matrix(dat),na.rm=TRUE)
  }else if (method == "mean_perSpatialUnit" | method=="knn_proteins_perSpatialUnit"){
    # extract column with spatial unit specification in the metadata
    spatUnits = unique(metadata[,spatial_unit_colname])
    imputed_data = c()
    for (s in 1:length(spatUnits)){
      ROI_indices = grep(spatUnits[s],colnames(dat))
      ROI_data = dat[,ROI_indices]
      if (method == "mean_perSpatialUnit"){
        imputed_ROI_data = impute_df(ROI_data,method="global_mean")
      }else if (method == "knn_proteins_perSpatialUnit"){
        imputed_ROI_data = impute_df(ROI_data,method="knn_proteins_global",knn_k,allowed_missingness_perProtein,allowed_missingness_perSample)
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
    knn_imputed_proteins_global = impute.knn(as.matrix(dat),k = knn_k, rowmax = allowed_missingness_perProtein, colmax = 1, maxp = dim(dat)[1], rng.seed=12345)
    # Fix immputed values for rows (proteins) that are missing more than rowmax proportion of samples.
    # By default, impute.knn uses the mean of all proteins in a sample for this case, which is not that meaningful.
    # Instead, we replace the missing values (for the case when missingness in a row is greater than rowmax) with global_mean or na for that row
    fix_prots = which( (rowSums(is.na(dat))/dim(dat)[2]) > allowed_missingness_perProtein )
    for (f in fix_prots){
      na_indices = which(is.na(dat[f,]))
      #knn_imputed_proteins_global$dat[f,na_indices] = mean(dat[f,],na.rm=TRUE)
      knn_imputed_proteins_global$dat[f,na_indices] = NA #Keep those missing values as NA since there is not enough dat to impute reliably.
    }
    data_imputed = knn_imputed_proteins_global$dat
  }else if (method == "knn_samples_global_proteinData"){
    transposed_data = as.matrix(t(dat))
    knn_imputed_samples_global_protData = impute.knn(transposed_data,k = knn_k, rowmax = allowed_missingness_perSample, colmax = 1, maxp = dim(transposed_data)[1], rng.seed=12345)
    fix_samples = which( (rowSums(is.na(transposed_data))/dim(transposed_data)[2]) > allowed_missingness_perSample )
    for (f in fix_samples){
      na_indices = which(is.na(transposed_data[f,]))
      #knn_imputed_proteins_global$data[f,na_indices] = mean(data[f,],na.rm=TRUE)
      knn_imputed_samples_global_protData$data[f,na_indices] = NA #Keep those missing values as NA since there is not enough data to impute reliably.
    }
    data_imputed = t(knn_imputed_samples_global_protData$dat)
  }else if (method == "knn_samples_global_spatialCoords"){
    d_for_nearestNeighb = as.matrix(metadata[,c("Xcoord","Ycoord")])
    # Find k-nearest neighbors based on spatial distance (x,y coordinates)
    knn <- knearneigh(d_for_nearestNeighb, k= knn_k)
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
  return (data_imputed)
}
