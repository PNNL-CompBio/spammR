test_that("diffex works", {
  pooled.panc.spe <- convert_to_spe(pooledData,
                                    pancMeta,
                                    protMeta,
                                    feature_meta_colname = "pancProts",
                                    sample_id = ""
  )
  diffex.spe <- calc_spatial_diff_ex(pooled.panc.spe,
    category_col = "IsletOrNot"
  )
  sigs <- subset(rowData(diffex.spe), NonIslet_vs_Islet.adj.P.Val.limma < 0.05)
  expect_equal(nrow(sigs),506)
})
