#' Generate a volcano plot to visualize results from differential expression analysis stored in a data frame.
#' Use this function directly if differential expression results are in a separate data frame and not an spe object
#' Helper function for volcanoPlot_DiffExSpe.R; use volcanoPlot_DiffExSpe.R if differential expression results are in an spe object.
#' -log10(p-value) vs. log2(fold change)
#' @import ggplot2
#' @export
#' @param diffEx_df A dataframe containing results from differential expression
#' @param logFC_colname column name in differenital expression results that represents log10 fold change
#' @param pval_colname column name in differential epxression results that represents the p-value to be plotted
#' @param pval_corrected Boolean indicating whether pval_colname represents the corrected p-value or not
#' @param title title for the plot
#' @param thresh Threshold for p-value to be used to annotate significant results on the plot
#' @param sigLabel_col Either a a vector with labels or a string (Example: "Gene") that is the column name in differential expression results that should be used for labeling significant results on the plot.
#' @return ggplot Volcano plot
volcanoPlot_DiffExResults <- function(diffEx_df, logFC_colname, pval_colname, pval_corrected, title, thresh=0.05, sigLabel_colname){
  require('ggrepel')
  diffEx_df = data.frame(diffEx_df)
  pval_type_lbl = ""
  if(pval_corrected){
    pval_type_lbl = "corrected p-value"
  }else{
    pval_type_lbl = "p-value"
  }
  yaxis_lbl = paste("-log10 (",pval_type_lbl,")",sep="")
  xdat = diffEx_df[,logFC_colname]
  ydat = -log10(diffEx_df[,pval_colname])
  # Plot -log10(p-value) vs. log2(fold change) and the significance level line
  sig_line_x = seq(min(xdat,na.rm = TRUE), max(xdat,na.rm = TRUE), by=0.1)
  sig_line_y = rep(-log10(thresh), length(sig_line_x))
  ymin = min(-log10(thresh), ydat)
  ymax = max(-log10(thresh), ydat)
  text_xcoord = (min(xdat) + max(xdat))/2
  pt_lbls = rep("",dim(diffEx_df)[1])
  sig_indices = which(diffEx_df[,pval_colname] < thresh)
  pt_lbls [sig_indices] = diffEx_df[sig_indices,sigLabel_colname]
  pt_colors = rep("black",dim(diffEx_df)[1])
  pt_colors[sig_indices] = rep("blue", length(sig_indices))
  #p_title = paste("Differential abundance analysis for ", group1_name," vs. ",group2_name,"\nP-values from t-test for each protein; Fold Change = (Mean for ",group1_name,")/(Mean for ",group2_name,")",sep="")
  p_title = title
  if (length(sig_indices > 0)){
    p_title = paste(p_title,"\nAnnotated points show ",sigLabel_colname," names",sep="")
  }
  text_xcoord = ( min(xdat,na.rm = TRUE) + max(xdat,na.rm=TRUE) ) / 2 #midpoint of x-axis
  p2 <- ggplot(diffEx_df, aes(get(logFC_colname), -log10(get(pval_colname)), label = pt_lbls)) +
    geom_point(color = pt_colors)+
    geom_hline(yintercept = -log10(thresh), color = "blue") +
    annotate("text", x=text_xcoord + 0.2*text_xcoord, y= (-log10(thresh) + 0.05), label=paste(pval_type_lbl," < ",thresh,sep=""), color = "blue")+
    geom_text_repel()+
    labs(title = p_title)+
    xlab("log2 (fold change)")+
    ylab(yaxis_lbl)
  #pdf(file = paste(plot_dir_curr,"/volcanoPlot2_DiffAbundance_ ",group1_name,"_vs_",group2_name,".pdf",sep=""),width=11, height=7)
  #print(p2)
  #dev.off()
  return(p2)
}
