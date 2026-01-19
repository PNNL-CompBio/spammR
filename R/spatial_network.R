#' Plot network of correlated features
#' @description `spatial_network()` Creates a network for a given 
#' SpatialExperiment object (spe) based on the correlation of features. Will 
#' return a fully connect graph
#' @details Using the feature values, the algorithm computes the correlation
#' between them across space, filtes by `lowest_thresh` (in the interest of
#' compute time) and then selects the most correlated edges using the
#' `top_edges` parameter. 
#' @details Only the largest connected component is plotted, but the entire 
#' graph is returned as an `tidygraph` object. 
#' @import ggplot2
#' @import ggraph
#' @import tidygraph
#' @import SpatialExperiment
#' @importFrom SingleCellExperiment altExp
#' @importFrom SingleCellExperiment altExpNames
#' @import dplyr
#' @import tidyr
#' @export
#' @param spe SpatialExperiment (SPE) object
#' @param assay_names Name of assay in the spe object that contains 
#' data to be plotted. If only one it will focus on features in the first 
#' assay. If more than one it will compute correlations across all features.
#' @param query_features List of features to query for correlation. If there 
#' are no `target_featres` will build a network between this feature all other
#' nodes. If missing, correlation will be computed between all nodes in 
#' `target_features`.
#' @param target_features List of features to compute correlation across. If
#' `query_features` is not NULL, will compute correlation only between these
#' features and the `query_features`. Otherwise, it will compute correlation
#' between all features in this list. DEFAULT: all features in spe. 
#' @param feature_names will describe which value to use for node names. IF
#' missing will use rownames
#' @param method to use for correlation, default is `spearman`
#' @return a `tidygraph` object that can be used for plotting or analysis
#' @examples
#' \dontrun{
#' data(pancMeta)
#' data(smallPancData)
#' data(protMeta)
#' img0.spe <- convert_to_spe(smallPancData$Image_0,
#'   pancMeta,
#'   protMeta,
#'   feature_meta_colname = "pancProts",
#'   image_files = system.file("extdata", "Image_0.png", package = "spammR"),
#'   spatial_coords_colnames = c("x_pixels", "y_pixels"),
#'   sample_id = "Image0",
#'   image_ids = "with_grid"
#' )
#' ##now call the network
#' res <- spatial_network(img0.spe,'proteomics')
#' }
spatial_network <- function(spe, 
                             assay_names,
                             feature_names = NULL, 
                             query_features = NULL,
                             target_features = list(), 
                             method = 'spearman'){
  
  ##TODO: check input
  stopifnot(is(spe,"SpatialExperiment"),
            !missing(assay_names))
  
  ##first let's combine any omes if they exist by iterating over the assays
  fullmat <- NULL
  full_fname <- NULL
  all_features <- target_features
  for (a in seq_along(assay_names)) {
    an <- assay_names[a] #get name of assay
    if (an %in% names(assays(spe))) { #first check if its in the primary dataset
      mval <- t(assay(spe,an))
      fname <- data.frame(rowval = colnames(mval), name = colnames(mval))
      if (!is.null(feature_names) &&
          feature_names[a] %in% colnames(rowData(spe)))
        fname$name <- rowData(spe)[colnames(mval),feature_names[a]]
      fname$class <- an
    }else if (length(altExpNames) > 0 &&
              an %in% names(assays(altExp(spe)))) { #its in the other experiment
      mval <- t(assay(altExp(spe), an))
      fname <- data.frame(rowval = colnames(mval), name = colnames(mval))
      if (!is.null(feature_names) && 
          feature_names[a] %in% colnames(rowData(altExp(spe))))
          fname$name <- rowData(altExp(spe))[colnames(mval),feature_names[a]]
      fname$class <- an
    }else{
      message(paste('Assay',an,'not in object'))
    }
    
    ##add target_features to list
    if (length(target_features) == 0) {
      all_features <- c(all_features,fname$name)
    }
    #add query_features to list
    if (!is.null(query_features) && query_features %in% fname$name) {
      all_features <- union(all_features, query_features)
    }
    
    fullmat <- cbind(fullmat,mval)
    full_fname <- rbind(full_fname, fname)
  }
  
  ##now reduce
  full_fname <- full_fname |> subset(name %in% all_features)
  fullmat <- fullmat[,intersect(colnames(fullmat), full_fname$rowval)]
  
  #default back to rowdata if we are missing feature names
  full_fname$name[which(is.na(full_fname$name))] <- 
    full_fname$rowval[which(is.na(full_fname$name))]
  
  #calculate correlations
  if (is.null(query_features)) {
    cormat <- cor(x = fullmat, method = method, use = 'p') |>
        as.data.frame(check.names = FALSE) 
  }else{
    qi <- subset(full_fname, name == query_features)
    ti <- subset(full_fname, name %in% target_features)
    cormat <- cor(x = fullmat[,qi$rowval], y = fullmat[,ti$rowval], 
                  method = method, use = 'p')
    
  }
  cormat <- data.frame(feature1 = full_fname$name, cormat, check.names = FALSE)
  #make into long format
  cordat <- cormat |>
      tidyr::pivot_longer(names_to = 'rowval', values_to = 'corval', 
                         cols = 2:ncol(cormat)) |>
    dplyr::left_join(full_fname) |>
    dplyr::rename(feature2 = 'name') |>
    dplyr::select(feature1, feature2, corval) |>
    dplyr::distinct() |>
    subset(!is.na(feature1)) |>
    subset(!is.na(feature2))
  
  #print(head(cordat))
   #get rid of self edges
  red <- cordat |>
    subset(!is.na(corval)) 
  
  #create graph
  g0 <- tidygraph::as_tbl_graph(red) 
  #add in data type
  
  g0 <- g0 %>%
    activate(edges) %>%
    filter(!edge_is_loop())
  
  g1 <- g0 |> 
    tidygraph::activate(nodes) |> 
    dplyr::left_join(dplyr::distinct(full_fname[,c('name','class')]), 
                     by = 'name')
  
  return(g1)
}
  