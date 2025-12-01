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
#' @import dplyr
#' @import tidyr
#' @export
#' @param spe SpatialExperiment (SPE) object
#' @param assay_names Name of assay in the spe object that contains 
#' data to be plotted. If only one it will focus on features in the first 
#' assay. If more than one it will compute correlations across all features.
#' @param features_of_interest List of features to focus on. If empty will 
#' compute correlation across all features 
#' @param feature_names will describe which value to use for node names. IF
#' missing will use rownames
#' @method Method to use for correlation, default is `spearman`
#' @return a `tidygraph` object that can be used for plotting or analysis

spatial_network <- function(spe, 
                             assay_names,
                             feature_names = NULL, 
                             features_of_interest = list(), 
                             method = 'spearman'){
  
  ##TODO: check input
  stopifnot(is(spe,"SpatialExperiment"),
            !missing(assay_names))
  
  ##first let's combine any omes if they exist by iterating over the assays
  fullmat <- NULL
  full_fname <- NULL
  for (a in seq_along(assay_names)) {
    an <- assay_names[a] #get name of assay
    if (an %in% names(assays(spe))) { #first check if its in the primary dataset
      mval <- t(assay(spe,an))
      fname <- data.frame(rowval = colnames(mval), name = colnames(mval))
      if (!is.null(feature_names) &&
          feature_names[a] %in% colnames(rowData(spe)))
        fname$name <- rowData(spe)[colnames(mval),feature_names[a]]
      fname$class <- an
    }else if (an %in% names(assays(altExp(spe)))) { #its in the other experiment
      mval <- t(assay(altExp(spe), an))
      fname <- data.frame(rowval = colnames(mval), name = colnames(mval))
      if (!is.null(feature_names) && 
          feature_names[a] %in% colnames(rowData(altExp(spe))))
          fname$name <- rowData(altExp(spe))[colnames(mval),feature_names[a]]
      fname$class <- an
    }else{
      message(paste('Assay',an,'not in object'))
    }
    if (length(features_of_interest) > 0) {
      mval <- mval[intersect(rownames(mval), features_of_interest),]
      fname <- fname |> subset(name %in% rownames(mval))
    }
    
    fullmat <- cbind(fullmat,mval)
    full_fname <- rbind(full_fname, fname)
  }
  #default back to rowdata if we are missing feature names
  full_fname$name[which(is.na(full_fname$name))] <- 
    full_fname$rowval[which(is.na(full_fname$name))]
  #calculate correlations
  cormat <- cor(x = fullmat, method = method, use = 'p') |>
      as.data.frame() 
  cormat <- cbind(feature1 = full_fname$name, cormat)
  #make into long format
  cordat <- cormat |>
      tidyr::pivot_longer(names_to = 'rowval', values_to = 'corval', 
                         cols = 2:ncol(cormat)) |>
    dplyr::left_join(full_fname) |>
    dplyr::rename(feature2 = 'name') |>
    dplyr::select(feature1, feature2, corval) |>
    dplyr::distinct() 
  
   #get rid of self edges
  red <- cordat |>
    subset(!is.na(corval)) 
  red$self <- vapply(seq_len(nrow(red)), 
                     function(x) {
                       red[[x,'feature1']] == red[[x,'feature2']]
                       }, logical(1))
  
  red <- red |>
    subset(!self)
  
  #create graph
  g0 <- tidygraph::as_tbl_graph(red) 
  #add in data type
  
  g1 <- g0 |> 
    tidygraph::activate(nodes) |> 
    dplyr::left_join(dplyr::distinct(full_fname[,c('name','class')]), by = 'name')
  
  return(g1)
}
  #   #count components
  # g2 <- g1 |>
  #   tidygraph::to_components()
  # 
  # ##get only the components that are greater than 2
  # csize <- vapply(g2,function(x) x %>% activate(nodes) %>% 
  #                   as_tibble() %>% 
  #                   nrow(), numeric(1))
  # 
  # print(paste('Graph has',length(g2),'components at current thresholds'))
  # 
  # cvals <- which(csize > 3)
  # g3 <- g2[[cvals[1]]]
  # for(i in cvals) {
  #   g3 <- g3 |>
  #     graph_join(g2[[i]])
  # }
  #   
  # pg <- ggraph(g2[[1]]) + 
  #   geom_edge_link(aes(colour = corval)) + 
  #   geom_node_point() + 
  #   geom_node_label(aes(label = name, colour = class))
  # 
  # print(pg)
  # 
  # return(g1)
#}