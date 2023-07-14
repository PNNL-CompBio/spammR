#' Generate a volcano plot to visualize results from differential expression analysis.
#' In spammR, these could be results obtained from the spatialDiffEx() function.
#' log10(p-value) vs. log2(fold change)
#' @export
#' @param spe SpatialExperiment object containing results from differential expression
#' @param logFC_col column name in differenital expression results that represents log10 fold change
#' @param pval_col column name in differential epxression results that represents the p-value to be plotted
#' @param ylab label for y-axis (Example: Corrected p-value)
#' @param title title for the plot
#' @param thresh Threshold for p-value to be used to annotate significant results on the plot
#' @param sigLabel_col Either a a vector with labels or a string (Example: "Gene") that is the column name in differential expression results that should be used for labeling significant results on the plot.
#' @return ggplot Volcano plot
volcanoPlot_DiffExSpe <- function(spe, logFC_col, pval_col, ylab, title, thresh=0.05, sigLabel_colname){
  df = rowData(spe)
  volcanoPlot_DiffExResults (df, logFC_col, pval_col, ylab, title, thresh, sigLabel_colname)
}