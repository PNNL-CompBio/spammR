#' Identify differentially abundant features across image
#' @description `calc_spatial_diff_ex()` Calculates differential expression
#' analysis using annotations in a SpatialExperiment object
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma voom
#' @import SpatialExperiment
#' @export
#' @param spe Spatial Experiment object containing data to be used for
#' differential expression analysis
#' @param assay_name Name of the dataset stored in the spe object, that is
#' to be used for the differential expression analysis.
#'  Example: znormalized_log2
#' @param count_based Set to TRUE of the data are count based, e.g. RNA-Seq
#' @param log_transformed Is the data given in spe log2 transformed TRUE
#' or FALSE
#' @param category_col Name of the column that specifies category of each
#' sample. Example: "IsletOrNot"
#' #Categories from `category_col` will be compared in the differential
#' expression analysis
#' @param compare_vals A vector containing names of categories from
#' category_col to be compared. Only required if there are more than two
#' values in `category_col`
#' @returns A Spatial Experiment object containing differential expression
#' results, stored in `rowData(diffEx.spe)`
#'  and `assays(diffEx.spe)` which contains the dataset on which differential
#' expresssion analysis was carried out
#'
#' @examples
#' data(smallPancData)
#' data(pancMeta)
#' data(protMeta)
#' pooledData <- dplyr::bind_cols(smallPancData)
#' pooled.panc.spe <- convert_to_spe(pooledData,
#'   pancMeta,
#'   protMeta,
#'   feature_meta_colname = "pancProts"
#' )
#' diffex.spe <- calc_spatial_diff_ex(pooled.panc.spe,
#'   category_col = "IsletOrNot"
#' )
#'
calc_spatial_diff_ex <- function(spe,
                                 assay_name = "proteomics",
                                 count_based = FALSE,
                                 log_transformed = FALSE,
                                 category_col,
                                 compare_vals) {
    # collect samples by factor
    factors <- unique(SummarizedExperiment::colData(spe)[[category_col]])
    if (length(factors) < 1) {
      ## throw error we need at least two categories
    } else if (length(factors) > 2) {
        if (missing(compare_vals) ||
              length(setdiff(compare_vals, factors)) > 0) {
        ## throw error - need exactly 2 values in category_col or 2 values in
        #category_col to compare
      }
      factors <- compare_vals
    }

    ##added to shorten line length
    f1 <- factors[1]
    f2 <- factors[2]

    ## now select the samples for each category
    samp1 <- which(SummarizedExperiment::colData(spe)[[category_col]] == f1)
    samp2 <- which(SummarizedExperiment::colData(spe)[[category_col]] == f2)
    # Later, limma call does samp2 vs. samp1 analysis
    comparison_name <- paste(factors[1], "_vs_", factors[2], sep = "")

    ## create design matrix with two factors
    fac <- factor(rep(c(2, 1), c(length(samp2), length(samp1))))
    design <- stats::model.matrix(~fac)
    # print(design)
    dat <- SummarizedExperiment::assays(spe)[[assay_name]]
    ldat <- dat
    if (!log_transformed) {
      ldat <- log2(dat)
    }

    if (count_based) {#call voom here
      ldat <- limma::voom(ldat, design)
    }

    fit <- limma::lmFit(ldat[, c(samp2, samp1)], design)
    fit <- limma::eBayes(fit)

    diffex.spe <- spe # initialize
    SummarizedExperiment::assays(diffex.spe) <- list()
    SummarizedExperiment::assays(diffex.spe,
                                 withDimnames = FALSE)[[assay_name]] <- dat
    # Output spe will only have one entry in assays, which will be the dataset
    # used for the differential expression analysis
    # withDimnames=FALSE drops rownames from dat while saving into the assay.
    # We want this because rownames(spe.out) may not necessarily have rownames,
    # depending on how the spe was defined.
    # Differential expression results will be stored in rowData(diffex.spe)
    # We deliberately do not sort the limma results (`sort.by = "none"`), such
    # that we can merge them correctly with the genes in the input SPE.
    res <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")

    ## message telling us how many genes are differentially expressedion
    nd <- subset(res, adj.P.Val < 0.05) |>
      subset(abs(logFC) > 1) |>
      nrow()
    msg <- paste("We found", nd,
                 "features with a logFC greater than 1 and \
                 an ajusted p-value less than 0.05")
    message(msg)

    colnames_res <- paste(paste(comparison_name,
                                colnames(res), "limma", sep = "."))

    colnames(res) <- colnames_res # c(feature_colname,colnames_res)
    # Make sure the results are in the same order of the features in the
    # input SPE object
    SummarizedExperiment::rowData(diffex.spe) <- cbind(rowData(diffex.spe),res)

    return(diffex.spe)
}
