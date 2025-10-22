test_that("can create spatial experiment", {

  # or we can add the inmage to a single data capture
  img0.spe <- convert_to_spe(smallPancData$Image_0,
    pancMeta,
    protMeta,
    feature_meta_colname = "pancProts",
    image_files = system.file("extdata", "Image_0.png", package = "spammR"),
    sample_id = "Image0",
    image_sample_ids = "Image0",
    image_ids = "Image0"
  )
  expect_equal(nrow(SummarizedExperiment::assay(img0.spe)), 2986)
})
