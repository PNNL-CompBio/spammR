#' enrich_ORA: Does over-representation analysis (ORA) using an interest list of genes from differential expression results in spammR and gene sets (either the ones provided in spammR or user supplied)
#' This function uses results from spatialDiffEx.R and assumes input of a specific format. Interest list of genes for ORA is obtained from spatialDiffEx results based on the criteria
#' specified in this function.
#' For ORA using an external or already defined interest list of genes and gene sets, use leapR functions directly
#' @export
#' @param spe SpatialExperiment object containing spatial omics data. This will be used to obtain background genes for ORA.
#' @param spatialDiffEx_results path to excel file containing results from spatialDiffEx.R
#' @param pval_type_forThresh Choose from "adjusted_pval" or "pval". Type of p-value that should be used for filtering statistically significant results. Default is adjusted p-value for multiple hypotheses correction.
#' @param pval_thresh value to use for filtering based on pval_type_forThreshold. Default is 0.05. Values less than pval_thresh will be kept.
#' @param logFC_lowerThresh Lower threshold for log Fold Change, to be used for filtering spatialDiffEx results. Default is NA
#' @param logFC_upperThresh Upper threshold for log Fold Change, to be used for filtering spatialDiffEx results. Default is NA
#' @param geneset_spammR Geneset to be used from spammR. Choose from: NA, "M2_curated_all" (default),"M2_cp_biocarta", "M2_cp_reactome","M5_GO_all", "M5_GO_bp", "M5_GO_cc", "M5_GO_mf"
#' @param geneset_external_type If not using a spammR provided geneset, specify type of geneset_external file. Options are ".gmt" or "list." If "list," make sure it is compatible with leapR format. Default is NA.
#' @param geneset_external If not using a spammR provided geneset, provide the filepath for the gene set file or the list object containing geneset information in leapR format. Default is NA.
#' @param ora_outdir Directory path for where results from ORA should be stored. Output is stored as an excel file as well as a dataframe as an Rdat file.
#' @param sortResultsBy For sorting ORA results, choose from the following column names: "BH_pvalue" (default)
#' @param comparison_name Example: "RSPv_vs_others" Text to indicate in results data frame, which spatial groups were compared for the interest list of genes
#' @returns a dataframe containing results from over-representation analysis of members of gene sets in the interest list of genes based on filtering criteria above.

enrich_ORA <-function(spe,spatialDiffEx_results,pval_type_forThresh, pval_thresh, logFC_lowerThresh=NA, logFC_upperThresh=NA, geneset_spammR="M2_curated_all",geneset_external_type=NA,geneset_external=NA,outdir,sortResultsBy,comparison_name){
  library(leapR)
  library(readxl)
  library(writexl)
  if (!is.na(geneset_spammR)){
    # Mouse geneset databases (downloaded from MSigDB as .gmt files) have been added to the data under spammR package
    #########################
    # Load appropriate mouse geneset databases
    #########################
    geneset_name = geneset_spammR
    geneset_path = paste("Geneset_Databases/mouse/misgdb/",geneset_spammR,"/",sep="")
    symbols.gmt_filename = list.files(path=geneset_path,pattern="sybmols.gmt")
    symbols.gmt_file = data(paste(geneset_path,"/",symbols.gmt_file,sep=""))
    geneset_leapR = read_gene_sets(symbols.gmt_file)
  }else{
    # add code for handling the case of an external geneset file or list
    geneset_name = "geneset_external"
  }
  # Filter spatialDiffEx results to create interest list of genes for ORA based on user-specified criteria
  # Current assumption is that spatialDiffEx results file has columns for protein names and corresponding gene names. "PG.genes"
  # If gene names are not present in spatialDiffEx results, a helper function (which I will add later), can be run to obtain gene names from the UniProt db files
  # and have a "PG.genes" column added to spatialDiffEx results excel file.
  sp_diffEx = data.frame(read_excel(spatialDiffEx_results))
  if (pval_type_forThresh=="adjusted_pval"){
    pval_col_text = "adj.P.Val"
  }else if (pval_type_forThresh=="pval"){
    pval_col_text = "P.Value"
  }else{
    # Throw error "Invalid value for pval_type_forThresh"
  }
  colnum_pvalue = grep(pval_col_text, colnames(sp_diffEx))
  colnum_logfc = grep("logFC",colnames(sp_diffEx))
  # Filter based on p-value criteria
  int_list = sp_diffEx[!is.na(sp_diffEx[,colnum_pvalue]),]
  int_list= int_list[int_list[,colnum_pvalue]< pval_thresh,]
  # Filter based on Log fold change criteria
  if (!is.na(logFC_lowerThresh)){
    int_list =int_list[int_list[,colnum_logfc]>logFC_lowerThresh,]
  }
  if (!is.na(logFC_upperThresh)){
    int_list =int_list[int_list[,colnum_logfc]<logFC_upperThresh,]
  }
  int_list_proteinNames = int_list$PG.ProteinNames
  int_list_geneNames = int_list$PG.Genes
  es = c()
  es_sorted = c()
  # Over-representation analysis: Run Enrichment in Sets (es) in leapR
  background_genes = rowData(spe)[,"PG.Genes"]
  es= leapR(geneset=geneset_leapR,
            enrichment_method='enrichment_in_sets',
            targets=int_list_geneNames, background = background_genes)
  # Sort results by column name specified by user
  es_sorted =  es[order( es[,"BH_pvalue"]),]
  num_sets = dim(es_sorted)[1]
  es_sorted = cbind(rep(comparison_name,num_sets), rownames(es_sorted), es_sorted)
  colnames(es_sorted)[1:2] = c("Comparison", geneset_name)
  # Output enrichment results for the ROI
  if (!dir.exists(ora_dir)){
    dir.create(ora_dir,recursive=TRUE)
  }
  es_file = paste("ORA_",comparison_name,".xlsx",sep="")
  int_list_filename = paste("intList_",comparison_name,".xlsx",sep="")
  write_xlsx(es_sorted,path=paste(ora_dir,es_file_name,sep="/"))
  write_xlsx(int_list,path=paste(ora_dir,int_list_filename,sep="/"))
  return(es_sorted)
}
