test_that("distance analysis works", {
    data(pancMeta)
    data(protMeta)
    data(smallPancData)
    img0.spe<-convert_to_spe(smallPancData$Image_0,
                            pancMeta,protMeta,
                            feature_meta_colname = 'pancProts',
                            spatial_coords_colnames = c('x_pixels','y_pixels'),
                            image_files = system.file("extdata",'Image_0.png',
                                                    package = 'spammR'),
                            sample_id  =  'Image0',
                            image_ids = 'Image0')
    img0.spe<-distance_based_analysis(img0.spe,'proteomics',
                                      sampleCategoryCol = 'IsletOrNot',
                                      sampleCategoryValue = 'Islet')
})
