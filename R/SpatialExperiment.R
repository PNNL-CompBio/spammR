#' The spatial experiment class
#'
#' Based on the \linkS4class{RangedSummarizedExperiment} and \linksS4class{SingleCellExperiment} classes, this class
#' maintains omics measurements for various sptailly related regions alongside
#' their x and y-coordinates
#'
#' @details
#'
#' @return SpatialExperiment object
#'



#' @docType class
#' @aliases
#' coerce,SpatialExperiment,SpatialExperiment-method
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods is as
#' @importFrom SingleCellExperiment
#' @importClassesFrom SingleCellExperiment
SpatialExperiment <- function(..., xCoord=c(), yCoord=c()){
  se <- SingleCellExperiment(...)
 # if(!is(se, "RangedSummarizedExperiment")) {
#    se <- as(se, "RangedSummarizedExperiment")
#  }
  .rse_to_sce(se, xCoord=xCoord,yCoord=yCoord)

}



#' @importFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
#' @importFrom BiocGenerics nrow ncol
.sce_to_spe <- function(sce, reducedDims=list(), altExps=list(), rowPairs=list(), colPairs=list(), mainExpName=NULL) {
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  out <- new("SpatialExperiment", sce,
             int_elementMetadata=new("DFrame", nrows=nrow(sce)),
             int_colData=new("DFrame", nrows=ncol(sce)))

  out
}

#' @exportMethod coerce
#' @importClassesFrom SingleCellExperiment
setAs("SingleCellExperiment", "SpatialExperiment", function(from) {
  .sec_to_spe(from)
})

