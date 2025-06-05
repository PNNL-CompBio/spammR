test_that('can create spatial experiment',{
    data(pancMeta)
    data(protMeta)
    data(smallPancData)
    #We can put all samples into the same object (for statistical power)
    pooledData <- dplyr::bind_cols(smallPancData)
    pooled.panc.spe <- convert_to_spe(pooledData,pancMeta,protMeta,feature_meta_colname = 'pancProts',samples_common_identifier = '')
    #or we can add the inmage to a single data capture
    img0.spe <- convert_to_spe(smallPancData$Image_0,
                           pancMeta,
                           protMeta,
                           feature_meta_colname = 'pancProts',
                           image_files = system.file("extdata",'Image_0.png',package = 'spammR'),
                           image_samples_common_identifier = 'Image0',
                           samples_common_identifier = 'Image0',
                           image_ids = 'Image0')
})