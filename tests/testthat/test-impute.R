test_that("imputation works", {
 data(smallPancData)
 data(pancMeta)
 data(protMeta)
 pooledData <- dplyr::bind_cols(smallPancData)
 pooled.panc.spe <- convert_to_spe(pooledData,
                 pancMeta,
                 protMeta,
                 feature_meta_colname = 'pancProts',
                 sample_id = '')
 
 resz <- impute_spe(pooled.panc.spe, method = 'zero')
 resm <- impute_spe(pooled.panc.spe, method = 'mean')
 resd <- impute_spe(pooled.panc.spe, method = 'median')
 resk <- impute_spe(pooled.panc.spe, method = 'knn',k = 4)
 resgm <- impute_spe(pooled.panc.spe, method = 'group_mean',group_colname = 'Image')
 resgk <- impute_spe(pooled.panc.spe, method = 'group_knn',k = 4,group_colname = 'Image')
 
 ##now create a single image
 img0.spe <- convert_to_spe(smallPancData$Image_0,
                            pancMeta,
                            protMeta,
                            feature_meta_colname = 'pancProts',
                            image_files = system.file("extdata",'Image_0.png',package = 'spammR'),
                            sample_id = 'Image0',
                            spatial_coords_colnames = c('x_pixels','y_pixels'),
                            image_sample_ids = 'Image0',
                            image_ids = 'Image0')

 ress <- impute_spe(img0.spe, method = 'spatial_knn',k = 3)
 
                 
})
