#' Identify multiomic features correlated with an distance to an image feature
#' @description distance_based_analysis: Identifies proteins/features that show a strong correlation between distance from a specified ROI's samples and
#' protein/feature abundance differences between samples.
#' @import SpatialExperiment
#' @importFrom IRanges IRanges
#' @importFrom IRanges countOverlaps
#' @importFrom stats cor.test
#' @export
#' @param spe SpatialExperiment object containing spatial omics data
#' @param assayName Name of the assay stored in spe that is to be used for distance based analysis. Example: "znormalized_log2"
#' @param spotHeightCol Column containing height of spot
#' @param spotWidthCol Column containing width of spot
#' @param sample_dimensions A vector containing the x and y dimensions of samples. Example: c(1,1) Sample Dimension units should match the units of spatial coordinates specified in spatialCoords(spe)
#' @param sampleCategoryCol Column name in metadata (colData(spe)) that should be used for selecting samples of certain type, to define the "origin" region for distance based analysis
#' @param sampleCategoryValue Sample category to be used for defining the "origin" region for distance based analysis
#' @param featuresNameCol Name of column containing features (example: proteins) in rowData(spe). It is assumed that the data provided in assay(spe) is in the same order as the order in which the features are listed under the featuresNameCol
#' @param corr_type Choose from "pearson" (default), "spearman." Correlation method to be used for calculating correlation between distance between samples and protein abundance differences. Both types of correlation provide a measure of monotonic association between two variables. Pearson is better suited for linear relationships while Spearman is better for non-linear monotonic relationships.
#' @param corr_thresh Minimum correlation value to be used for identifying proteins that have a correlation between protein abundance differences and distance between samples. Values greater than or equal to this theshold will be used.
#' @param min_samples of a minimum number of sample points for calculating correlation. For proteins with less than this number of sample points, correlation value is reported as NA.
#' @param allowOverlaps allow overlaps of regions
#' @returns dist_based_results List containing a.) pairwise_calculations_betweenSamples which is a list; each entry corresponds to a feature and contains a dataframe consisting of a distance (between samples) vector and a corresponding protein abundance difference vector.
#' b.) corr_pval_all and c.) corr_pval_thresholded which are data frames containing the correlation value, p-value and number of sample points used for the correlation calculation for each feature. corr_pval_thresholded only contains results with correlation > corr_thresh.
#' dist_based_results is also saved as .RData object under results_dir.
#' @examples
#' data(pancData)
#' data(pancMeta)
#' data(protMeta)
#' data(pancDataList)
#' pooled.panc.spe <- convert_to_spe(pancData,pancMeta,protMeta,feature_meta_colname='pancProts',samples_common_identifier='')
#' img0.spe<-convert_to_spe(pancDataList$Image_0,pancMeta,protMeta,feature_meta_colname='pancProts',image_files=system.file("extdata",'Image_0.png',package='spammR'),image_samples_common_identifier='Image0',samples_common_identifier = 'Image0',image_ids='Image0')

#' img0.spe<-distance_based_analysis(img0.spe,'proteomics',sampleCategoryCol='IsletOrNot',sampleCategoryValue='Islet')

distance_based_analysis <- function(spe,
                                    assayName,
                                    spotHeightCol='spot_height',
                                    spotWidthCol='spot_width',
                                    sampleCategoryCol,
                                    sampleCategoryValue,
                                    featuresNameCol,
                                    corr_type="spearman",
                                    corr_thresh = 0.5,
                                    min_samples=5,
                                    allowOverlaps=TRUE){



  # Compute centroids for each sample based on top-left corner (Xcoord, Ycoord) coordinates (SG: should be bottom!?!)
 # sample_dim_x = sample_dimensions[1]
#  sample_dim_y = sample_dimensions[2]
  spatial_coords = data.frame(spatialCoords(spe))

  ##first we use Iranges to check to see if there are overlapping areas (in x-space)
  ##SG: Maybe remove this...
  if(!allowOverlaps){

    xranges = IRanges::IRanges(start=spatial_coords[,1],width=colData(spe)[[spotWidthCol]])|>
        unique()
    yranges = IRanges::IRanges(start=spatial_coords[,2],width=colData(spe)[[spotHeightCol]])|>
      unique()
    allovers = c(IRanges::countOverlaps(xranges),IRanges::countOverlaps(yranges))
    if(any(allovers)>1){
      message(print(paste(length(allovers),"overlapping regions found, consider setting allowOverlaps to TRUE")))
      exit()
    }
  }

  ##get centers of points
  centroid_x = as.numeric(spatial_coords[,1] + colData(spe)[[spotWidthCol]]/2)#sample_dim_x/2)
  centroid_y = as.numeric(spatial_coords[,2] + colData(spe)[[spotHeightCol]]/2)#- sample_dim_y/2)
  centroid_coords = data.frame(cbind(centroid_x,centroid_y))
  rownames(centroid_coords) = rownames(spatial_coords)
  # Compute distance between samples
  dist_between_samples = as.matrix(stats::dist(centroid_coords, method = "euclidean", diag=TRUE))

  assay_data = SummarizedExperiment::assays(spe)[[assayName]]
  #rownames(assay_data) = SummarizedExperiment::rowData(spe)[,featuresNameCol]
  source_samples_indices = which(SummarizedExperiment::colData(spe)[,sampleCategoryCol]==sampleCategoryValue)
  if(length(source_samples_indices)>1){
    message('More than one feature to measure distance from, not supported yet')
    exit()
  }

  ##get the distances
  samp_dists <- dist_between_samples[source_samples_indices,]

 ##correlate the distances
  cor_dist <- cor(t(assay(spe)),samp_dists,method=corr_type,use='pairwise.complete.obs')
  cor_test <- apply(assay(spe),1,function(x){
    pval=1.0
    try(pval<-stats::cor.test(x,samp_dists,method=corr_type,exact=FALSE)$p.val,silent=TRUE)
    return(pval)})

  mdata<-data.frame(cv=cor_dist,cp=cor_test)
  colnames(mdata)<-c(paste0(sampleCategoryValue,'Distance',corr_type,'Cor'),
                     paste0(sampleCategoryValue,'Distance',corr_type,'Pval'))

  rowData(spe)<-cbind(rowData(spe),mdata)

#  source_samples = colnames(assay_data)[source_samples_indices]

 # source_samples_data = assay_data[, source_samples_indices]
  # Remove features that don't have data for any source samples
#  remove_rows = which(is.na(source_samples_data))#which(rowSums(is.na(source_samples_data)) == length(source_samples))
#  new_dat = assay_data[-remove_rows,]


  # # Initialize vectors for storing results
  # corr_dist_protAbund_Diffs = c()
  # dist_AbundDiffs_vectors = list()
  # pval_corr_dist_protAbund_Diffs = c()
  # num_samplePoints_forCorrCalc = c()
  #
  #
  # for (j in rownames(new_dat)){ # for each protein/feature
  #   abund_currProt = new_dat[j,] # protein abundance values for current protein/feature
  #   # Call helper function calcl_diff_pairwiseSamples() that we defined outside of this script, for computing pairwise differences between sample values
  #   diffAbund = calc_diff_pairwiseSamples(abund_currProt)
  #   diffAbund_sourceSamples = data.frame(diffAbund[source_samples,]) # Abundance differences between source samples and all other samples
  #   dist_betw_samples_sourceSamples = data.frame(dist_between_samples[source_samples,]) # Distances between source samples and other samples
  #   diffVec = c() # Store all protein abund difference values between different pairwise combinations of samples; Re-initialize for each protein
  #   distVec = c() # Store all corresponding distances between samples; Re-initialize for each protein
  #   sample_i = c() # for documentation
  #   sample_j = c() # for documentation
  #   for (i in 1:dim(diffAbund_sourceSamples)[1]){
  #     diffVec = c(diffVec, t(diffAbund_sourceSamples[i,]))
  #     distVec = c(distVec, t(dist_betw_samples_sourceSamples[i,]))
  #     sample_i = c(sample_i, rep(rownames(diffAbund_sourceSamples)[i], length(diffAbund_sourceSamples[i,])) )
  #     sample_j = c(sample_j, colnames(diffAbund_sourceSamples) )
  #   }
  #   num_samplePoints = sum(!is.na(diffVec))
  #   data_used = rep(assayName, length(diffVec))
  #   dist_AbundDiffs_vectors[[ j ]] = data.frame(distVec,diffVec, sample_i, sample_j, data_used)
  #   colnames(dist_AbundDiffs_vectors[[ j ]]) = c("Distance_between_samples","Abundance_difference_between_samples (sample_i - sample_j)","sample_i", "sample_j","Data_used")
  #   dist_calcs_filepath = ""
  #
  #       if (num_samplePoints >= min_samplePoints_forCorr){ # Must have at least the specified number of values to calculate correlation
  #     print(j)
  #     print(length(diffVec))
  #     print(length(distVec))
  #     corr_results = stats::cor.test(distVec, diffVec, method=corr_type,use="pariwise.complete.obs")
  #     corr_dist_protAbund_Diffs [j] = corr_results$estimate
  #     pval_corr_dist_protAbund_Diffs [j] = corr_results$p.value
  #   }else{
  #     corr_dist_protAbund_Diffs [j] = NA
  #     pval_corr_dist_protAbund_Diffs [j] = NA
  #   }
  #   num_samplePoints_forCorrCalc [j] = num_samplePoints
  # }
  # origin_samples = rep(sampleCategoryValue,length(corr_dist_protAbund_Diffs))
  # correlation_type = rep(corr_type,length(corr_dist_protAbund_Diffs))
  # corr_and_pval_all = data.frame(origin_samples,corr_dist_protAbund_Diffs,pval_corr_dist_protAbund_Diffs,correlation_type,num_samplePoints_forCorrCalc)
  # corr_and_pval_all = corr_and_pval_all[order(corr_and_pval_all$corr_dist_protAbund_Diffs,decreasing=TRUE),]
  # corr_and_pval_thresholded = corr_and_pval_all[which(corr_and_pval_all$corr_dist_protAbund_Diffs>corr_thresh),]
  # dist_based_results = list()
  # dist_based_results [["pairwise_calculations_betweenSamples"]] = dist_AbundDiffs_vectors
  # dist_based_results [["corr_and_pval_thresholded"]] = corr_and_pval_thresholded
  # dist_based_results [["corr_and_pval_all"]] = corr_and_pval_all
  #save(dist_based_results,file=paste(results_dir2,"/",sampleCategoryValue,"_dist_based_results_",corr_type,".RData",sep=""))
  return(spe)
  #plot(distVec,diffVec,xlab=paste("Distance from samples in ",sampleCategoryValue,"\n1 distance unit = 200 um",sep=""),ylab=paste("Protein abundance gradient realtive to samples in ",sampleCategoryValue,sep=""))
}
