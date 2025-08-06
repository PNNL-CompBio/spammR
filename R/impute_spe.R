#' Impute missing values based on spatial coordinates
#' @description `impute_spe()` carries out imputation for missing data in a data frame (df) using a specified method from a range of methods. Accepts input dataset as a data frame (df).
#' @importFrom matrixStats rowMedians
#' @importFrom impute impute.knn
#' @importFrom spdep knearneigh
#' 
#' @param spe SPE containing data Data frame containing data to be imputed, where rows correspond to features (which can be specified as row names of dat but it is not required that they be specified) and columns correspond to samples. Column names must correspond to names provided in the sample identifier column in the metadata parameter.
#' @param assay_name name of assay with data to be imputed
#' @param method Method of imputation to be used. See details.
#' @param group_colname Column name in metadata that specifies the group information to use for group_mean or knn_group. Example: ROI_abbreviation.
#' @param k K value to be used for k-nearest neighbor imputation
#' @param protein_missingness Proportion of samples allowed to have missing data for a protein in the given spatial unit specified by the imputation method. Example: When method="global_mean," an allowed_missingness_perProtein of 0.5 indicates that any protein missing data for more than 50% of samples across the entire spatial tissue covered by all samples will be excluded from the imputation method algorithm, and that protein's missing values will not be imputed. When method="mean_perSpatialUnit," then allowed_missingness_perProtein of 0.5 indicates that a protein must have data for at least 50% of samples in the specified spatial unit (example: a brain ROI) to be used in the imputation algorithm and to be imputed. This argument is only needed for methods that rely on statistics based on per protein. It is used for the following methods: median_a, median_b, global_mean, mean_perSpatialUnit, knn_proteins_global and knn_proteins_perSpatialUnit.
#' @param sample_missingness  Proportion of proteins allowed to be missing for a sample. This argument is only needed when imputing using k-nearest neighbor methods where neighbors are samples. For all other imputation methods provided here, allowed_missingness_perSample does not apply; no samples are excluded from the imputation process for those methods. If missingness in a sample is greater than this thershold, then that sample's data are not included in the imputation for the knn samples neighbors methods.
#' @returns An SPE with an 'imputed' data frame with the appropriate imputation called
#' @export 
#' 
#' @details Methods options and descriptions:
#' - zero: replace missing values with 0
#' - median : replace missing values with global median per protein
#' - median_half : replace missing values with 1/2 global median per protein
#' - mean: replace missing values with global mean per protein
#' - group_mean: replace msising values with mean per group, e.g. group (example: ROI) for each protein
#' - knn: imputation based on k-nearest neighbors, with proteins as neighbors, based on data from all samples across all groups 
#' - group_knn: imputation based on k-nearest neighbors, with proteins as neighbors, based on data from specified group (e.g. ROI, tissue)
#' - spatial_knn: imputation based on k-nearest neighors in space
#' @examples
#' 
#' data(pancMeta)
#' data(protMeta)
#' data(smallPancData)
#' #We can put all samples into the same object (for statistical power)
#' pooledData<-dplyr::bind_cols(smallPancData)
#' pooled.panc.spe <- convert_to_spe(pooledData,
#'                 pancMeta,
#'                 protMeta,
#'                 feature_meta_colname = 'pancProts',
#'                 sample_id='')
#' res <- impute_spe(pooled.panc.spe, method='mean')

#.methods <- 

impute_spe <- function(spe,
                      assay_name = 'proteomics',
                      method =  NULL ,
                      group_colname,
                      k=NULL,
                      protein_missingness=NULL,
                      sample_missingness=NULL){

    
  # # Set defult values if not specified by user
  # if (is.null(allowed_missingness_perProtein)){ # Only used for knn methods
  #   protein_missingness = 0.75
  # }
  # if (is.null(allowed_missingness_perSample)){ # Only used for knn methods
  #   protein_missingness = 0.80
  # }
    .methods = c('zero','median','median_half','mean','group_mean','knn','group_knn','spatial_knn')
    
    ##check for args
    if(!method %in% .methods)
        stop(paste('Method',method, 'must be one of',paste(.methods,collapse=',')))
    ##check for K in knn
    if(method %in% c('knn','group_knn','spatial_knn') && is.null(k))
        stop(paste('Method',method,'requires parameter k to be set'))
    
    dat <- SummarizedExperiment::assays(spe)[[assay_name]]
  metadata <- SummarizedExperiment::colData(spe)
  spcoords <- SpatialExperiment::spatialCoords(spe)

  imputed_data = c()
  replace_vals = c()
  
  ##first lets see if there are proteins that we dont impute
  fix_prots = which((rowSums(is.na(dat)) / ncol(dat)) > protein_missingness )
  
  ##and samples
  fix_samps = which((colSums(is.na(dat)) / nrow(dat)) > sample_missingness )
  
  
  
  ###first iterate through the global methods
  if(method%in%c('zero','mean','median','median_half')){
      ##row-wise values first
      if (method == "zero") {
        replace_vals = rep(0,nrow(dat))
    }else if (method == "mean") {
        replace_vals = rowMeans(dat,na.rm = TRUE)
    }else if (method == "median") {
        replace_vals = matrixStats::rowMedians(as.matrix(dat),na.rm = TRUE)
    }else if (method == "median_half") {
        replace_vals = 0.5 * matrixStats::rowMedians(as.matrix(dat),na.rm = TRUE)
    }
     imputed_data <- dat
     ###now do the imputation
     for(i in setdiff(1:nrow(dat),fix_prots)){ ##for all of the values we WANT to fix
         imputed_data[i,which(is.na(dat[i,]))] <- replace_vals[i]
     }
     
     ##if we were mising things
  }
  else if (method == "knn") {
      kres = impute::impute.knn(as.matrix(dat),k = k, rowmax = protein_missingness, colmax = 1, maxp = dim(dat)[1], rng.seed=12345)
      # Fix immputed values for rows (proteins) that are missing more than rowmax proportion of samples.
      # By default, impute.knn uses the mean of all proteins in a sample for this case, which is not that meaningful.
      # Instead, we replace the missing values (for the case when missingness in a row is greater than rowmax) with global_mean or na for that row
      imputed_data <- kres$dat
      #print(fix_prots)
      for (f in fix_prots) {
          na_indices = which(is.na(dat[f,]))
          #knn_imputed_proteins_global$dat[f,na_indices] = mean(dat[f,],na.rm=TRUE)
          imputed_data[f,na_indices] = NA #Keep those missing values as NA since there is not enough dat to impute reliably.
      }
      if (length(fix_prots)>0) {
          message("impute.knn function uses mean imputation for rows with more than the specified missingness. spammR's imputation function reverts these missing values back to NA.")
      }
#      data_imputed = knn_imputed_proteins_global$dat]
  }else if (method == "group_mean" | method == "group_knn") {
    # extract column with spatial unit specification in the metadata
    spatUnits = unique(metadata[,group_colname]) ##whats is the spatial unit?
    imputed_data = dat
    ##iterate through each group, identify means, and replace those values
    for (s in 1:length(spatUnits)) {
      ROI_indices = which(metadata[,group_colname] == s)
      #        grep(spatUnits[s],colnames(dat)) ##FIX: spatial unit should not be in column name
      if (method == "group_mean") {
        replace_vals = rowMeans(dat[,ROI_indices],na.rm = TRUE)
          
        for(i in setdiff(1:nrow(imputed_data),fix_prots)){ ##for all of the values we WANT to fix
            imputed_data[i,which(is.na(dat[i,ROI_indices]))] <- replace_vals[i] ##replace the values WITHIN THIS GROUP
        }
      } else if (method == "group_knn") {
            gdat <- dat[,ROI_indices]
            imputed_gdata = impute::impute.knn(as.matrix(gdat),k = k, rowmax = protein_missingness, colmax = 1, maxp = dim(gdat)[1], rng.seed=12345)
            imputed_data[,ROI_indices] <- imputed_gdata$dat
     }
    # Make sure the imputed data appears in the same order as the original data
    imputed_data = imputed_data[,colnames(dat)]
    } 
   }else if (method == "spatial_knn"){
       if(ncol(spcoords)==0)
           stop('Need spatial coordinates for spatial_knn method')
    d_for_nearestNeighb = as.matrix(spcoords)
    # Find k-nearest neighbors based on spatial distance (x,y coordinates)
    knn <- spdep::knearneigh(d_for_nearestNeighb, k= k)
    # nearest neighbor indices are in knn$nn
    nn = knn$nn
    for (prot in 1:dim(dat)[1]){
      this.prot = dat[prot,]
      na_indices = which(is.na(this.prot))
      for (n in na_indices){
        n.knn = nn[n,] #nearest neighbors for sample corresponding to nth index
        this.prot[n] = mean(as.numeric(this.prot[n.knn]),na.rm=TRUE)
      }
      imputed_data = rbind(imputed_data,this.prot)
    }
    #imputed_data <- data_imputed[rownames(dat),]
  }
  ##TODO: figure out missingness values for knn code

  # 
  # if (is.null(replace_vals)){
  #   #don't change data_imputed
  # }else{
  #   data_imputed = c() #start with this, then replace na values
  #   # Update replace_vals to reflect allowed_missingness_perProtein criteria.
  #   # replace_vals should retain NAs (missing values) for proteins above the allowed_missingness_perProtein criteria for the following methods
  #   dontImpute_rows = c()
  #   if (method=="median" | method=="median_half" | method=="global_mean"){
  #     dontImpute_rows = which ( (rowSums(is.na(dat))/dim(dat)[2]) > allowed_missingness_perProtein)
  #   }
  #   replace_vals[dontImpute_rows] = NA
  #   for (i in 1:nrow[1]){
  #     this.row = dat[i,]
  #     this.row.na_indices = which(is.na(this.row))
  #     if (length(this.row.na_indices)>0){
  #       if (length(replace_vals)>1){
  #         this.row[this.row.na_indices] <- replace_vals[i]
  #       }else{
  #         this.row[this.row.na_indices] <- replace_vals
  #       }
  #     }
  #     data_imputed = rbind(data_imputed,this.row)
  #   }
  # }
  # data_imputed = data.frame(data_imputed, row.names=rownames(dat))
  # colnames(data_imputed) = colnames(dat)
  # return(data_imputed)
  
  
  assays(spe, withDimnames=FALSE) <- c(assays(spe),list(imputed=imputed_data))
  return(spe)
}

