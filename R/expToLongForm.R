

#' expToLongForm: helper function that moves expression
#' matrix to long form
#' @export
expToLongForm<-function(sce,rowname='prot'){

  exprs(sce)%>%
    as.matrix()%>%
    as.data.frame()%>%
    tibble::rownames_to_column(rowname)%>%
    tidyr::pivot_longer(c(2:(1+ncol(exprs(spat.phos)))),names_to='Voxel',values_to='LogRatio')
}
