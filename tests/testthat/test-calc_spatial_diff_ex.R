test_that("diffex works", {
 data(smallPancData)
 data(pancMeta)
 data(protMeta)
 pooledData <- dplyr::bind_cols(smallPancData)
 pooled.panc.spe <- convert_to_spe(pooledData,
                 pancMeta,
                 protMeta,
                 feature_meta_colname = 'pancProts',
                 sample_id = '')
 
 diffex.spe <- calc_spatial_diff_ex(pooled.panc.spe,
                 category_col = 'IsletOrNot')
                 
})
