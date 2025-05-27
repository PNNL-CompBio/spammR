#' Here is the data to be included in the spammR package
#'

#'
#' Sample metadata for the `pancData` file.
#' @references \url{panc paper}
"pancMeta"

#'
#' Protein metadata for general uniprot data, including the pancProts column for the panc data
#' @source Uniprot
"protMeta"

#'
#' Pancreatic proteomics measurements from Gosline et al. proteomics study.
#' Row names indicate individual proteins, columns indicate sample identifiers described in
#' `pancMeta`. Individual values show un-transformed abundance values for each protein.
#' List of proteomic measurements in a list sorted from individual regions of the pancreas.
#' @source \url{panc paper}
"pancDataList"

#' 
#' Soil samples from the 1000 soils project
#' Row names indiciate KEGG Orthologs, Column names indicate locations
#' @source \url{}
'oneKSoilsData'

'oneKSoilsMeta'

'oneKSoilsCoords'