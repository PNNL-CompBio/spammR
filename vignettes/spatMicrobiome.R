## ----loading-packages, warning=FALSE, message=FALSE---------------------------
library(spammR)
library(ggplot2)
## these two libraries are only used for geographic data
library(sf)
library(terra)

## ----oneksoils, echo=FALSE, warning=FALSE, message=FALSE----------------------
sd <- download.file("https://api.figshare.com/v2/file/download/55158824", 
                    dest = "sdata.rda", quiet = TRUE)
load("sdata.rda")
file.remove("sdata.rda")

## ----oneksoils meta-----------------------------------------------------------
oneKSoilsOmicsMeta <- data.frame(KO = character(nrow(oneKSoilsData)))

oneKSoilsOmicsMeta["KO"] <- rownames(oneKSoilsData)
oneKSoilsOmicsMeta["annotation"] <- seq(1, by = 1, 
                                        length.out = nrow(oneKSoilsOmicsMeta))
oneKSoilsOmicsMeta["annotation"] <- lapply(oneKSoilsOmicsMeta["annotation"], 
                                           as.character)
head(oneKSoilsOmicsMeta)

## ----coordiantes--------------------------------------------------------------
data("oneKSoilsMeta")
head(oneKSoilsMeta)

ggplot2::ggplot(oneKSoilsMeta, aes(x = Core_Layer, y = pH, fill = biome_name)) +
  geom_boxplot()

## ----import-usshape, warning=FALSE, message=FALSE, eval=FALSE-----------------
# us_map <- sf::st_read("../data/tl_2024_us_state/tl_2024_us_state.shp",
#                       quiet = TRUE)
# us_map

## ----import-usmap-from-data-file, echo=FALSE, warning=FALSE, message=FALSE----
maprda <- download.file("https://api.figshare.com/v2/file/download/55158815", 
                        quiet = TRUE, destfile = "map.rda")
load("map.rda")
file.remove("map.rda")

## ----subsetting-us-map, warning=FALSE, message=FALSE--------------------------
us_continental <- c(
  "WV", "FL", "IL", "MN", "MD", "RI", "ID", "NH", "NC", "VT",
  "CT", "DE", "NM", "CA", "NJ", "WI", "OR", "NE", "PA", "WA",
  "LA", "GA", "AL", "UT", "OH", "TX", "CO", "SC", "OK", "TN",
  "WY", "ND", "KY", "ME", "NY", "NV", "MI", "AR", "MS", "MO",
  "MT", "KS", "IN", "SD", "MA", "VA", "DC", "IA", "AZ"
)

s_us_cont <- us_map %>%
  .[.$STUSPS %in% us_continental, ] %>%
  vect()
s_us_cont

## ----raster-template, warning=FALSE, message=FALSE----------------------------
template <- terra::rast(
  s_us_cont,
  res = 0.01
)

## -----------------------------------------------------------------------------
# us_raster <- rasterize(s_us_cont, template, background=0)
# us_raster <- rasterize(as.lines(s_us_cont), template, field='STATEFP', 
# touches=TRUE)

us_raster <- terra::rasterize(s_us_cont, template, field = "STATEFP")

## -----------------------------------------------------------------------------
oneKSoilsMeta <- oneKSoilsMeta[!is.na(oneKSoilsMeta$latitude) & 
                                 !is.na(oneKSoilsMeta$longitude), ]
coords <- sf::st_as_sf(oneKSoilsMeta, coords = c("longitude", "latitude"), 
                       remove = TRUE)
sf::st_crs(coords) <- sf::st_crs(s_us_cont)
coords

## -----------------------------------------------------------------------------
coords_rast <- terra::rasterize(coords, template)
coords_rast

## -----------------------------------------------------------------------------
cells <- terra::extract(coords_rast, coords, xy = TRUE, cells = TRUE)
print(head(cells))

## -----------------------------------------------------------------------------
cells_tmp <- terra::extract(coords_rast, coords, xy = TRUE, cells = TRUE)
xy_rast <- lapply(cells["cell"], function(x) rowColFromCell(coords_rast, x))
cells <- cells_tmp %>%
  cbind(., oneKSoilsMeta, xy_rast) %>%
  dplyr::rename(., y_pixels = "cell.1", x_pixels = "cell.2") %>%
  subset(., select = -c(x, y, ID, last))
cells$y_pixels <- nrow(coords_rast) - cells$y_pixels

## ----write-raster-map, eval=TRUE----------------------------------------------
terra::writeRaster(us_raster, "us_map.png", NAflag = 0, overwrite = TRUE)

## -----------------------------------------------------------------------------
coords_test <- sf::st_as_sf(oneKSoilsMeta, 
                            coords = c("longitude", "latitude"), 
                            remove = TRUE)
sf::st_crs(coords_test) <- sf::st_crs(s_us_cont)
plot(us_raster, col = "grey")
# lines(as.polygons(us_raster), col='black')
plot(coords_test, col = "darkred", add = TRUE)

## -----------------------------------------------------------------------------
dm <- cells
rownames(dm) <- dm$Sample_ID
# dm[,1] <- NULL
dm[, "cell"] <- NULL
dm[, "longitude.1"] <- NULL
dm[, "latitude.1"] <- NULL
dm[, "Image"] <- 0
dm[, "x_max"] <- ncol(coords_rast)
dm[, "y_max"] <- nrow(coords_rast)
dm[, "x_origin"] <- 0
dm[, "y_origin"] <- 0
dm[, "spot_height"] <- 500
dm[, "spot_width"] <- 500
dm[, "above_40_deg_lat"] <- dplyr::if_else(dm[, "latitude"] >= 40, 
                                                  true = 1, false = 0)
dm[, "psychrophile"] <- dplyr::if_else(dm[, "Soil.Temperature.C"] <= 15, 
                                       true = 1, false = 0)
dm[, "Core_Layer_bin"] <- dplyr::if_else(dm[, "Core_Layer"] == "TOP", 
                                         true = "TOP", false = "BTM")
dm[, "Desert"] <- dplyr::if_else(dm[, "BIOME"] == "13", 
                                 true = "Desert", false = "NonDesert")

data_meta <- dm
head(data_meta)

## -----------------------------------------------------------------------------
data_meta <- subset(data_meta, Sample_ID %in% colnames(oneKSoilsData))
oneKSoilsData <- subset(oneKSoilsData, select = data_meta$Sample_ID)
data_meta <- data_meta[data_meta$Core_Layer == "TOP", ]
oneKSoilsData <- oneKSoilsData[, grep("TOP", x = names(oneKSoilsData))]

## -----------------------------------------------------------------------------
microbiome.spe <- spammR::convert_to_spe(
  oneKSoilsData, ## pooled data table
  data_meta, ## pooled metadata
  oneKSoilsOmicsMeta, ## protein identifiers
  feature_meta_colname = "KO", # column name
  spatial_coords_colnames = c("x_pixels", "y_pixels"),
  image_files = c("us_map.png"),
  image_sample_ids = "1000 Soils",
  sample_id = "1000 Soils",
  image_ids = c("0"),
  assay_name = "KO"
)

## -----------------------------------------------------------------------------
diff_ex <- spammR::calc_spatial_diff_ex(
  microbiome.spe,
  assay_name = "KO",
  log_transformed = FALSE,
  category_col = "Desert"
)

sig_ko <- subset(rowData(diff_ex), Desert_vs_NonDesert.adj.P.Val.limma < 0.1)
ups <- subset(sig_ko, Desert_vs_NonDesert.logFC.limma > 0)
downs <- subset(sig_ko, Desert_vs_NonDesert.logFC.limma < 0)

print(paste(
  "We found", nrow(sig_ko), "significantly differentally abundant KOs",
  nrow(ups), "upregulated KOs and", nrow(downs), "downregulated"
))

print(ups)

## -----------------------------------------------------------------------------
spatial_heatmap(
  microbiome.spe,
  feature = "K00005",
  # feature_type = 'KO',
  assay_name = "KO",
  sample_id = "1000 Soils",
  image_id = "0",
  label_column = "Desert",
  metric_display = "KO association abundance",
  plot_title = "Abundance of glycerol dehydrogenase (K00005)",
  interactive = FALSE
)

## -----------------------------------------------------------------------------
kk <- download.file("https://api.figshare.com/v2/file/download/55158818", 
                    destfile = "kegg.rda", quiet = TRUE)
load("kegg.rda")
file.remove("kegg.rda")

## -----------------------------------------------------------------------------
ora.res <- spammR::enrich_ora(diff_ex, geneset = keggPathwayToKo, 
                              geneset_name = "KEGG pathways", 
                              feature_column = "KO", pval_thresh = 0.1)
print(head(ora.res))

## -----------------------------------------------------------------------------
ora.res <- tibble::rownames_to_column(ora.res, "pathways")
pathway_prots <- subset(ora.res, pathways == "map00561") |>
  dplyr::select(ingroupnames)
pathway_prots <- unlist(strsplit(pathway_prots[1, ], split = ", "))

## -----------------------------------------------------------------------------
spatial_heatmap(
    microbiome.spe,
    feature = pathway_prots,
    assay_name = "KO",
    sample_id = "1000 Soils",
    image_id = "0",
    label_column = "Desert",
    metric_display = "Enriched KO Groups",
    plot_title = "Glycerolipid metabolism expression",
    interactive = FALSE
)

## ----session------------------------------------------------------------------
sessionInfo()

