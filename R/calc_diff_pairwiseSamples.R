#' Helper function for computing pairwise differences between sample values
#' This function is called in distance_based_analysis()
#' @param feature_values A 1 x n vector containing values for a feature for n samples. names(feature_values) should contain names of samples
calc_diff_pairwiseSamples = function(feature_values){
  pairwiseDiff = data.frame()
  for (i in 1:length(feature_values)){
    for (j in 1:length(feature_values)){
      pairwiseDiff[i,j] =  feature_values[i] - feature_values[j]
    }
  }
  rownames(pairwiseDiff) = names(feature_values)
  colnames(pairwiseDiff) = names(feature_values)
  return(pairwiseDiff)
}
