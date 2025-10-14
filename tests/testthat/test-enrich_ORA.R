test_that("ORA enrichment works", {

  panc.spe <- convert_to_spe(pooledData, pancMeta, protMeta, feature_meta_colname = "pancProts", sample_id = "")
  diffex.spe <- calc_spatial_diff_ex(panc.spe, category_col = "IsletOrNot")
  library(leapR)
  data("krbpaths")
  ora.res <- enrich_ora(diffex.spe, geneset = krbpaths, geneset_name = "krbpaths", feature_column = "PrimaryGeneName")
  sigs <- subset(ora.res, BH_pvalue < 0.05)
  expect_equal(nrow(sigs), 154)
  
})
