#' The spatial experiment class
#'
#' Based on the \linkS4class{RangedSummarizedExperiment} and \linksS4class{SingleCellExperiment} classes, this class
#' maintains omics measurements for various sptailly related regions alongside
#' their x and y-coordinates
#'
#' @details
#' Add details here
#' @examples
#' data(pancData)
#' data(pancMeta)
#' spe<-SpatialExperiment(assays=list(logcounts=as(pancData),'dgCMatrix'),colData=pancMeta)
#' @return SpatialExperiment object
#'


#' @export
#' @rdname SpatialExperiment
#' @importFrom utils packageVersion
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("SpatialExperiment",
         slots=list(coordinates="data.frame"),#,
         #int_metadata = "list"),
         contains = "SingleCellExperiment",
         prototype = prototype(
           int_metadata=list(
             version=packageVersion("spammer"),
             mainExpName=NULL
           )
         )
)



#' @name SpatialExperiment
#' @docType class
#' @aliases
#' coerce,SingleCellExperiment,SpatialExperiment-method
#' @importFrom S4Vectors SimpleList
#' @importFrom methods is as
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
 SpatialExperiment <- function(..., coords){
   se <- SingleCellExperiment(...)
  # if(!is(se, "RangedSummarizedExperiment")) {
 #    se <- as(se, "RangedSummarizedExperiment")
 #  }
   .sce_to_spe(se, coordinates=coords)

 }



#' @name .sce_to_spe
#' @importFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
#' @importFrom BiocGenerics nrow ncol
.sce_to_spe <- function(sce, coordinates=data.frame()) {
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  out <- new("SpatialExperiment", sce)

  coordinates(out)<-coordinates

  out
}

#' @name setAs
#' @exportMethod coerce
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setAs("SingleCellExperiment", "SpatialExperiment", function(from) {
  .sce_to_spe(from)
})

