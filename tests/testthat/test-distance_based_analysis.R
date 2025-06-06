test_that("distance analysis works", {
    data(pancMeta)
    data(protMeta)
    data(smallPancData)
    img0.spe<-convert_to_spe(smallPancData$Image_0,
                            pancMeta,protMeta,
                            feature_meta_colname = 'pancProts',
                            image_files = system.file("extdata",'Image_0.png',
                                                    package = 'spammR'),
                            image_samples_common_identifier = 'Image0',
                            samples_common_identifier  =  'Image0',
                            image_ids = 'Image0')
    img0.spe<-distance_based_analysis(img0.spe,'proteomics',
                                      sampleCategoryCol = 'IsletOrNot',
                                      sampleCategoryValue = 'Islet')
})
