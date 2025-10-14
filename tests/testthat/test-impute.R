test_that("imputation works", {

  pooled.panc.spe <- convert_to_spe(pooledData,
    pancMeta,
    protMeta,
    feature_meta_colname = "pancProts",
    sample_id = ""
  )

  resz <- impute_spe(pooled.panc.spe, method = "zero")
  resm <- impute_spe(pooled.panc.spe, method = "mean")
  resd <- impute_spe(pooled.panc.spe, method = "median")
  resk <- impute_spe(pooled.panc.spe, method = "knn", k = 4)
  resgm <- impute_spe(pooled.panc.spe, method = "group_mean", group_colname = "Image")
  resgk <- impute_spe(pooled.panc.spe, method = "group_knn", k = 4, group_colname = "Image")

  ## now create a single image
  img0.spe <- convert_to_spe(smallPancData$Image_0,
    pancMeta,
    protMeta,
    feature_meta_colname = "pancProts",
    image_files = system.file("extdata", "Image_0.png", package = "spammR"),
    sample_id = "Image0",
    spatial_coords_colnames = c("x_pixels", "y_pixels"),
    image_sample_ids = "Image0",
    image_ids = "Image0"
  )

  expect_equal(round(assay(resz,'imputed')['sp|A0PJW6|TM223_HUMAN','0_S_1_3']),0)
  expect_equal(round(assay(resm,'imputed')['sp|A0PJW6|TM223_HUMAN','0_S_1_3']),11)
  expect_equal(round(assay(resd,'imputed')['sp|A0PJW6|TM223_HUMAN','0_S_1_3']),12)
  expect_equal(round(assay(resk,'imputed')['sp|A0PJW6|TM223_HUMAN','0_S_1_3']),12)
  expect_equal(round(assay(resgm,'imputed')['sp|A0PJW6|TM223_HUMAN','0_S_1_3']),11)
  expect_equal(round(assay(resgk,'imputed')['sp|A0PJW6|TM223_HUMAN','0_S_1_3']),12)
  
  ress <- impute_spe(img0.spe, method = "spatial_knn", k = 3)
  expect_equal(round(assay(ress,'imputed')['sp|A0PJW6|TM223_HUMAN','0_S_1_3']),11)
  
})
