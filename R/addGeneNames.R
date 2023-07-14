#' Obtain gene names for each feature (currently, protein) in a given dataset in an spe object, add a column with gene names to the dataset and return a new spe object
#' Run this function before running any geneset based analysis functions if the starting data doesn't have a gene name column already.
#' Gene names are needed for any pathway/geneset based analysis as the members in each geneset are given as gene names, not protein identifiers.
#' @export
#' @param spe Spatial Experiment object containing the omics dataset
#' @param featureID_colname Column name containing feature names for which gene names are needed (Example: "Protein_rowID")
#' @param species Name of the species whose data is in spe (example: "human")
#' @param db_file Path or name of database file containing information for converting from feature ID to gene name. Defaults to database files provided in spammR for the given species
#' @return spe_withGeneNames Spatial Experiment object containing everything in the input spe, along with a column for gene name.
addGeneNames <- function(spe, featureID_colname, species, db_file){
  retrun(spe)
}
