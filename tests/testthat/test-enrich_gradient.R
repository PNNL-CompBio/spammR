test_that("gradient enrichment works", {

  # We can put all samples into the same object (for statistical power)
  # or we can add the inmage to a single data capture
  img0.spe <- convert_to_spe(smallPancData$Image_0,
                             pancMeta, protMeta,
                             feature_meta_colname = "pancProts",
                             spatial_coords_colnames = c("x_pixels", "y_pixels"),
                             image_files = system.file("extdata", "Image_0.png",
                                                       package = "spammR"
                             ),
                             sample_id = "Image0",
                             image_ids = "Image0"
  )
  img0.spe <- distance_based_analysis(img0.spe, "proteomics",
                                      sampleCategoryCol = "IsletOrNot",
                                      sampleCategoryValue = "Islet"
  )
  
  
  rank.res <- enrich_gradient(img0.spe,
     geneset = krbpaths,
     feature_column = "PrimaryGeneName",
     ranking_column = "IsletDistancespearmanCor"
   )
  
  sigs <- subset(rank.res,BH_pvalue < 0.05)
  expect_equal(nrow(sigs), 12)
})
