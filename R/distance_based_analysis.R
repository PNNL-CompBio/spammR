#' distance_based_analysis: Identifies proteins/features that show a strong correlation between distance from a specified ROI's samples and
#' protein/feature abundance differences between samples.

#' @import SpatialExperiment
#' @export
#' @param spe SpatialExperiment object containing spatial omics data
#' @param assayName Name of the assay stored in spe that is to be used for distance based analysis. Example: "znormalized_log2"
#' @param sample_dimensions A vector containing the x and y dimensions of samples. Example: c(1,1) Sample Dimension units should match the units of spatial coordinates specified in spatialCoords(spe)
#' @param sampleCategoryCol Column name in metadata (colData(spe)) that should be used for selecting samples of certain type, to define the "origin" region for distance based analysis
#' @param sampleCategoryValue Sample category to be used for defining the "origin" region for distance based analysis
#' @param corr_type Choose from "pearson" (default), "spearman." Correlation method to be used for calculating correlation between distance between samples and protein abundance differences. Both types of correlation provide a measure of monotonic association between two variables. Pearson is better suited for linear relationships while Spearman is better for non-linear monotonic relationships.
#' @param corr_thresh Minimum correlation value to be used for identifying proteins that have a correlation between protein abundance differences and distance between samples. Values greater than or equal to this theshold will be used.
#' @param min_samplePoints_forCorr Requirement of a minimum number of sample points for calculating correlation. For proteins with less than this number of sample points, correlation value is reported as NA.
#' @param results_dir Path of output directory where results from this function should be saved
#' @returns dist_based_results List containing a.) pairwise_calculations_betweenSamples which is a list; each entry corresponds to a feature and contains a dataframe consisting of a distance (between samples) vector and a corresponding protein abundance difference vector.
#' b.) corr_pval_all and c.) corr_pval_thresholded which are data frames containing the correlation value, p-value and number of sample points used for the correlation calculation for each feature. corr_pval_thresholded only contains results with correlation > corr_thresh.
#' dist_based_results is also saved as .RData object under results_dir.


distance_based_analysis <- function(spe,assayName,sample_dimensions,sampleCategoryCol, sampleCategoryValue,corr_type="pearson",corr_thresh,min_samplePoints_forCorr=6,results_dir){
  library(SpatialExperiment)
  # Compute centroids for each sample based on top-left corner (Xcoord, Ycoord) coordinates
  sample_dim_x = sample_dimensions[1]
  sample_dim_y = sample_dimensions[2]
  spatial_coords = data.frame(spatialCoords(spe))
  centroid_x = as.numeric(spatial_coords$x + sample_dim_x/2)
  centroid_y = as.numeric(spatial_coords$y - sample_dim_y/2)
  centroid_coords = data.frame(cbind(centroid_x,centroid_y))
  rownames(centroid_coords) = rownames(spatial_coords)
  # Compute distance between samples
  dist_between_samples = as.matrix(dist(centroid_coords, method = "euclidean", diag=T))

  assay_data = assays(spe)[[assayName]]
  source_samples_indices = which(colData(spe)[,sampleCategoryCol]==sampleCategoryValue)
  source_samples = colnames(assay_data)[source_samples_indices]
  source_samples_data = assay_data[, source_samples_indices]
  # Remove features that don't have data for any source samples
  remove_rows = which(rowSums(is.na(source_samples_data)) == length(source_samples))
  new_dat = assay_data[-remove_rows,]
  # Initialize vectors for storing results
  corr_dist_protAbund_Diffs = c()
  dist_AbundDiffs_vectors = list()
  pval_corr_dist_protAbund_Diffs = c()
  num_samplePoints_forCorrCalc = c()
  results_dir2 = paste(results_dir,"/dataType",assayName,"/",corr_type,"_correlation/",sampleCategoryValue,sep="")
  results_dir3 = paste(results_dir2,"/calculations_per_feature",sep="")
  if (!dir.exists(results_dir3)){
    dir.create(results_dir3,recursive=TRUE)
  }
  for (j in rownames(new_dat)){ # for each protein/feature
    abund_currProt = new_dat[j,] # protein abundance values for current protein/feature
    # Call helper function calcl_diff_pairwiseSamples() that we defined outside of this script, for computing pairwise differences between sample values
    diffAbund = calc_diff_pairwiseSamples(abund_currProt)
    diffAbund_sourceSamples = data.frame(diffAbund[source_samples,]) # Abundance differences between source samples and all other samples
    dist_betw_samples_sourceSamples = data.frame(dist_between_samples[source_samples,]) # Distances between source samples and other samples
    diffVec = c() # Store all protein abund difference values between different pairwise combinations of samples; Re-initialize for each protein
    distVec = c() # Store all corresponding distances between samples; Re-initialize for each protein
    sample_i = c() # for documentation
    sample_j = c() # for documentation
    for (i in 1:dim(diffAbund_sourceSamples)[1]){
      diffVec = c(diffVec, t(diffAbund_sourceSamples[i,]))
      distVec = c(distVec, t(dist_betw_samples_sourceSamples[i,]))
      sample_i = c(sample_i, rep(rownames(diffAbund_sourceSamples)[i], length(diffAbund_sourceSamples[i,])) )
      sample_j = c(sample_j, colnames(diffAbund_sourceSamples) )
    }
    num_samplePoints = sum(!is.na(diffVec))
    data_used = rep(assayName, length(diffVec))
    dist_AbundDiffs_vectors[[ j ]] = data.frame(distVec,diffVec, sample_i, sample_j, data_used)
    colnames(dist_AbundDiffs_vectors[[ j ]]) = c("Distance_between_samples","Abundance_difference_between_samples (sample_i - sample_j)","sample_i", "sample_j","Data_used")
    # if (grep(";",j)){ # Feature name has multiple names (this was a special case for the brain ROI data)
    #   j_shortened =
    # }
    write_xlsx(dist_AbundDiffs_vectors[[j]], path=paste(results_dir3,"/",sampleCategoryValue,"_",j,"_dist_calcs.xlsx",sep=""))
    if (num_samplePoints >= min_samplePoints_forCorr){ # Must have at least the specified number of values to calculate correlation
      print(j)
      print(length(diffVec))
      print(length(distVec))
      corr_results = cor.test(distVec, diffVec, method=corr_type,use="pariwise.complete.obs")
      corr_dist_protAbund_Diffs [j] = corr_results$estimate
      pval_corr_dist_protAbund_Diffs [j] = corr_results$p.value
    }else{
      corr_dist_protAbund_Diffs [j] = NA
      pval_corr_dist_protAbund_Diffs [j] = NA
    }
    num_samplePoints_forCorrCalc [j] = num_samplePoints
  }
  origin_samples = rep(sampleCategoryValue,length(corr_dist_protAbund_Diffs))
  correlation_type = rep(corr_type,length(corr_dist_protAbund_Diffs))
  corr_and_pval_all = data.frame(origin_samples,corr_dist_protAbund_Diffs,pval_corr_dist_protAbund_Diffs,correlation_type,num_samplePoints_forCorrCalc)
  corr_and_pval_all = corr_and_pval_all[order(corr_and_pval_all$corr_dist_protAbund_Diffs,decreasing=TRUE),]
  corr_and_pval_thresholded = corr_and_pval_all[which(corr_and_pval_all$corr_dist_protAbund_Diffs>corr_thresh),]
  dist_based_results = list()
  dist_based_results [["pairwise_calculations_betweenSamples"]] = dist_AbundDiffs_vectors
  dist_based_results [["corr_and_pval_thresholded"]] = corr_and_pval_thresholded
  dist_based_results [["corr_and_pval_all"]] = corr_and_pval_all
  return(dist_based_results)
  save(dist_based_results,file=paste(results_dir2,"/",sampleCategoryValue,"_dist_based_results_",corr_type,".RData",sep=""))
  #plot(distVec,diffVec,xlab=paste("Distance from samples in ",sampleCategoryValue,"\n1 distance unit = 200 um",sep=""),ylab=paste("Protein abundance gradient realtive to samples in ",sampleCategoryValue,sep=""))
}
