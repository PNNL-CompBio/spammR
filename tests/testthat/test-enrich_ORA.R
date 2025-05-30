test_that("ORA enrichment works", {
    data(pancDataList)
    data(pancMeta)
    data(protMeta)
    pooledPanc <- dplyr::bind_cols(pancDataList)
    panc.spe <- convert_to_spe(pooledPanc,pancMeta,protMeta,feature_meta_colname = 'pancProts',samples_common_identifier = '')
    diffex.spe <- calc_spatial_diff_ex(panc.spe,category_col = 'IsletOrNot')
    library(leapR)
    data('krbpaths')
    ora.res <- enrich_ora(diffex.spe, geneset = krbpaths,geneset_name = 'krbpaths', feature_column = 'PrimaryGeneName')
    
})
