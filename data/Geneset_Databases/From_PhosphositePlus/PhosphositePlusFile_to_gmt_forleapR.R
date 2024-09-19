# Load needed packages and libraries
# install.packages("devtools")
# devtools::install_github("biodataganache/leapR")
library("leapR")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("writexl")
library("writexl")
#install.packages("readxl")
library("readxl")

# Load kinase substrate dataset
indir = "/Users/sohi472/OneDrive - PNNL/Projects/Exploration_data_tools/leapR/pathwayDatabases/From_PhosphositePlus"
kinase_file = "Kinase_Substrate_Dataset.xls" #has only one sheet
kinase_info = data.frame(read_excel(paste(indir,kinase_file,sep="/")))
#Filter dataframe for mouse data
kinase_info_mouse = kinase_info[kinase_info$KIN_ORGANISM=="mouse" & kinase_info$SUB_ORGANISM == "mouse",]
kinase_to_genes = hash()
kinase_to_genes_length = hash()
for (i in 1:dim(kinase_info_mouse)[1]){
  curr_kinase = kinase_info_mouse$KINASE[i]
  substrate_gene = kinase_info_mouse$SUB_GENE[i]
  if (is.null(kinase_to_genes[[ curr_kinase ]])){
    l = list()
    l[[ 1 ]] = substrate_gene
  }else{
    l = append(kinase_to_genes[[ curr_kinase ]], substrate_gene)
  }
  kinase_to_genes[[ curr_kinase ]] = unique(l)
  kinase_to_genes_length[[ curr_kinase ]] = length(unique(l))
}
# Output a file in the gmt format to be used in leapR functions later
# Format: 
# kinase_name, a column with description of kinase (use "desc" for placeholder), followed by a list of substrate genes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cmapR")
library(cmapR)
# Can use write_gmt(lst, fname) function to create the gmt file
# "lst needs to be a nested list where each sub-list is itself a list with the following fields: 
# - head: the name of the data 
# - desc: description of the corresponding data 
# - len: the number of data items 
# - entry: a vector of the data items"
kinase_list = list()
kinases = keys(kinase_to_genes)
for (i in 1:length(kinases)){
  kinase_list[[ kinases[i] ]] = list()
  kinase_name = paste("kinase_", kinases[i],sep="")
  kinase_list[[ kinases[i] ]] [[ "head" ]] = kinase_name
  kinase_list[[ kinases[i] ]] [[ "desc" ]] = "desc"
  kinase_list[[ kinases[i] ]] [[ "len" ]] = kinase_to_genes_length[[ kinases[i] ]]
  kinase_list[[ kinases[i] ]] [[ "entry" ]] = unlist(kinase_to_genes[[ kinases[i] ]])
}
kinase_list_fname = paste(indir, "/PhosphositePlusFile_to_gmt/kinase_substrate_mouse.gmt", sep="")
write_gmt(kinase_list, kinase_list_fname)
# checked output file; looks great and as expected!
