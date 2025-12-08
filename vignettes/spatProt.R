## ----load package, message=FALSE,warning=FALSE--------------------------------
## load spammR
library(spammR)

## ----omics data---------------------------------------------------------------
download.file("https://api.figshare.com/v2/file/download/55158821",
              mode = "wb", quiet = TRUE, dest = "pdl.rda")
load("pdl.rda")
utils::head(pancDataList$Image_0[, 1:8])
file.remove("pdl.rda")

## ----data list----------------------------------------------------------------
print(length(pancDataList))
head(pancDataList[[2]][, 1:8])

## ----samp meta----------------------------------------------------------------
data(pancMeta)
head(pancMeta)

## ----img, eval=FALSE----------------------------------------------------------
# library(cowplot)
# 
# cowplot::ggdraw() + cowplot::draw_image(system.file("extdata",
#                                                     "Image_1.png",
#                                                     package = "spammR"))

## ----omics meta---------------------------------------------------------------
data(protMeta)
head(protMeta[, c("pancProts", "EntryName", "PrimaryGeneName")])

## ----spe combined-------------------------------------------------------------
pooledData <- dplyr::bind_cols(pancDataList)
pooled.panc.spe <- convert_to_spe(pooledData, ## pooled data table
  pancMeta, ## pooled metadata
  protMeta, ## protein identifiers
  feature_meta_colname = "pancProts", # column name
)
print(pooled.panc.spe)

## ----spe list, echo=FALSE, warning=FALSE, message=FALSE-----------------------
## list of image names
imglist <- c("Image_0", "Image_1", "Image_2")
           
img.spes <- lapply(
    imglist,
    function(x) {
        convert_to_spe(pancDataList[[x]],
          pancMeta,
          protMeta,
          sample_id = x,
          feature_meta_colname = "pancProts",
          spatial_coords_colnames = c("x_pixels", "y_pixels"),
          image_files = system.file("extdata",
            paste0(x, ".png"),
            package = "spammR"
          ),
          image_sample_ids = x,
          image_ids = "with_grid"
        )
    }
)
names(img.spes) <- imglist

## ----spatial data-------------------------------------------------------------
allimgs <- lapply(imglist, function(x) {
    spe <- img.spes[[x]]
    res <- spatial_heatmap(spe,
      feature = "INS",
      feature_type = "PrimaryGeneName",
      sample_id = x,
      image_id = "with_grid",
      label_column = "IsletOrNot",
      interactive = FALSE
  )
  return(res)
})

allimgs[[2]]

## ----diffex, warning=FALSE, error=FALSE, message=FALSE------------------------
islet_res <- calc_spatial_diff_ex(pooled.panc.spe,
    assay_name = "proteomics",
    log_transformed = FALSE,
    category_col = "IsletOrNot"
)

# we filter the significant proteins first
sig_prots <- subset(rowData(islet_res), 
                    NonIslet_vs_Islet.adj.P.Val.limma < 0.01)
# then separate into up-regulated and down-regulated based on fold chnage
ups <- subset(sig_prots, NonIslet_vs_Islet.logFC.limma > 0)
downs <- subset(sig_prots, NonIslet_vs_Islet.logFC.limma < 0)

print(paste(
  "We found", nrow(sig_prots), "significantly differentally \
  expressed proteins including",
  nrow(ups), "upregulated proteins and", nrow(downs), "downregulated"
))

## ----plot pathway, message=FALSE, warning=FALSE, error=FALSE------------------
spe.plot <- img.spes[[2]]

hup <- spatial_heatmap(spe.plot,
    feature = rownames(ups),
    sample_id = "Image_1",
    image_id = "with_grid",
    label_column = "IsletOrNot",
    interactive = FALSE
)

hup


## ----ora analysis-------------------------------------------------------------
library(leapR)
data("krbpaths")
ora.res <- enrich_ora(islet_res, geneset = krbpaths, 
                      geneset_name = "krbpaths", 
                      feature_column = "PrimaryGeneName")
print(ora.res[grep("INSULIN", rownames(ora.res)), 
              c("ingroup_n", "pvalue", "BH_pvalue")])

## ----pathway plotting,warning=FALSE-------------------------------------------
secprots <- ora.res["REACTOME_GLUCOSE_REGULATION_OF_INSULIN_SECRETION", ] |>
    dplyr::select(ingroupnames) |>
    unlist() |>
    strsplit(split = ", ") |>
    unlist()

spe.plot <- img.spes[[2]]

hup <- spatial_heatmap(spe.plot,
    feature = secprots,
    sample_id = "Image_1",
    image_id = "with_grid",
    feature_type = "PrimaryGeneName",
    label_column = "IsletOrNot",
    plot_title = "Glucose regulation proteins",
    interactive = FALSE
)

## ----distance, warning=FALSE, error=FALSE, message=FALSE----------------------
## for each image, let's compute the distance of each voxel to 
## the one labeled 'Islet'
rank.imgs <- lapply(
  img.spes,
  function(x) {
    distance_based_analysis(x, "proteomics",
      sampleCategoryCol = "IsletOrNot",
      sampleCategoryValue = "Islet"
    )
  }
)

## now we have the distances, let's plot some interesting proteins
negProts <- do.call(rbind, lapply(names(rank.imgs), function(x) {
    subset(as.data.frame(rowData(rank.imgs[[x]])), 
           IsletDistancespearmanPval < 0.01) |>
      subset(IsletDistancespearmanCor < (-.75)) |>
      dplyr::select(PrimaryGeneName, IsletDistancespearmanCor) |>
      dplyr::mutate(image = x)
}))

print(head(negProts))

## do any proteins show up more than once?
icounts <- negProts |>
    dplyr::group_by(PrimaryGeneName) |>
    dplyr::summarize(numImgs = dplyr::n()) |>
    dplyr::arrange(desc(numImgs))


print(icounts)

## ----plot correlated proteins, message=FALSE----------------------------------
spatial_heatmap(img.spes[[3]],
    feature = "SH3GL1",
    feature_type = "PrimaryGeneName",
    sample_id = names(img.spes)[3],
    image_id = "with_grid",
    label_column = "IsletOrNot", interactive = FALSE
)

spatial_heatmap(img.spes[[2]],
    feature = "AP3S2",
    feature_type = "PrimaryGeneName",
    sample_id = names(img.spes)[2],
    image_id = "with_grid",
    label_column = "IsletOrNot", interactive = FALSE
)

## ----rank based enrichment,warning=FALSE,error=FALSE--------------------------
library(leapR)
data("krbpaths")
enriched.paths <- do.call(rbind, lapply(names(rank.imgs), function(x) {
  spe <- rank.imgs[[x]]
  es <- enrich_gradient(spe,
    geneset = krbpaths,
    feature_column = "PrimaryGeneName", # mapped to enrichment data
    ranking_column = "IsletDistancespearmanCor"
  )
  es[, "comp"] <- rep(x, nrow(es))
  es[, "krbpaths"] <- rownames(es)
  es
}))

enriched.paths |>
  subset(BH_pvalue < 0.05) |>
  dplyr::group_by(krbpaths) |>
  dplyr::summarize(numImgs = dplyr::n()) |>
  dplyr::arrange(desc(numImgs))

## ----gradient plotting--------------------------------------------------------
rprots <- subset(enriched.paths, krbpaths == "KEGG_RIBOSOME") |>
    dplyr::select(comp, ingroupnames)

rprots <- unlist(strsplit(rprots[1, 2], split = ", "))

spatial_heatmap(img.spes[[2]],
    feature = rprots,
    feature_type = "PrimaryGeneName",
    sample_id = names(img.spes)[2],
    image_id = "with_grid",
    label_column = "IsletOrNot", interactive = FALSE
)

