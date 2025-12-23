library(spammR)

##load other libraries

##there are additional dependencies
if (!require('tidyverse')) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require('readxl')) {
  install.packages('readxl')
  library(readxl)
}

if (!require('geojsonR')) {
  install.packages('geojsonR')
  library(geojsonR)
}

if (!require('tidygraph')) {
  install.packages('tidygraph')
  library(tidygraph)
}

if (!require('ggraph')) {
  install.packages('ggraph')
  library(ggraph)
}

if (!require('cowplot')) {
  install.packages('cowplot')
  library(cowplot)
}
library(leapR)
library(sf)
library(terra)
library(geojsonR)

##panel 1: insulin plot

download.file("https://api.figshare.com/v2/file/download/55158821",
              mode = "wb", quiet = TRUE, dest = "pdl.rda")
load("pdl.rda")
file.remove("pdl.rda")
data(pancMeta)
data(protMeta)
imglist <- c("Image_0", "Image_1", "Image_2")

x = 'Image_1'
i0 <- convert_to_spe(pancDataList[[x]],
                   pancMeta,
                   protMeta,
                   sample_id = x,
                   feature_meta_colname = "pancProts",
                   spatial_coords_colnames = c("x_pixels", "y_pixels"),
                   image_files = system.file("extdata",
                                             paste0(x, ".png"),
                                             package = "spammR"),
                                             image_sample_ids = x,
                                             image_ids = "with_grid")
                   
res <- spatial_heatmap(i0,
                       feature = "INS",
                       feature_type = "PrimaryGeneName",
                       sample_id = x,
                       image_id = "with_grid",
                       label_column = "IsletOrNot",
                       interactive = FALSE)

ggsave('fig_panel1.pdf',res)

##panel 2: microbiome plot

#ko measurements
sd <- download.file("https://api.figshare.com/v2/file/download/55158824", 
                    dest = "sdata.rda", quiet = TRUE)
load("sdata.rda")
file.remove("sdata.rda")

#meadata
data("oneKSoilsMeta")

oneKSoilsOmicsMeta <- data.frame(KO = character(nrow(oneKSoilsData)))

oneKSoilsOmicsMeta["KO"] <- rownames(oneKSoilsData)
oneKSoilsOmicsMeta["annotation"] <- seq(1, by = 1, 
                                        length.out = nrow(oneKSoilsOmicsMeta))
oneKSoilsOmicsMeta["annotation"] <- lapply(oneKSoilsOmicsMeta["annotation"], 
                                           as.character)

maprda <- download.file("https://api.figshare.com/v2/file/download/55158815", 
                        quiet = TRUE, destfile = "map.rda")
load("map.rda")
file.remove("map.rda")

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

template <- terra::rast(
  s_us_cont,
  res = 0.01
)

us_raster <- terra::rasterize(s_us_cont, template, field = "STATEFP")

oneKSoilsMeta <- oneKSoilsMeta[!is.na(oneKSoilsMeta$latitude) & 
                                 !is.na(oneKSoilsMeta$longitude), ]
coords <- sf::st_as_sf(oneKSoilsMeta, coords = c("longitude", "latitude"), 
                       remove = TRUE)
sf::st_crs(coords) <- sf::st_crs(s_us_cont)
coords_rast <- terra::rasterize(coords, template)


cells <- terra::extract(coords_rast, coords, xy = TRUE, cells = TRUE)
print(head(cells))

cells_tmp <- terra::extract(coords_rast, coords, xy = TRUE, cells = TRUE)
xy_rast <- lapply(cells["cell"], function(x) rowColFromCell(coords_rast, x))
cells <- cells_tmp %>%
  cbind(., oneKSoilsMeta, xy_rast) %>%
  dplyr::rename(., y_pixels = "cell.1", x_pixels = "cell.2") %>%
  subset(., select = -c(x, y, ID, last))
cells$y_pixels <- nrow(coords_rast) - cells$y_pixels
terra::writeRaster(us_raster, "us_map.png", NAflag = 0, overwrite = TRUE)

coords_test <- sf::st_as_sf(oneKSoilsMeta, 
                            coords = c("longitude", "latitude"), 
                            remove = TRUE)
sf::st_crs(coords_test) <- sf::st_crs(s_us_cont)


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

spatial_heatmap(
  microbiome.spe,
  feature = "K00005",
  # feature_type = 'KO',
  assay_name = "KO",
  sample_id = "1000 Soils",
  image_id = "0",
  sample_label_size = 2.0, 
  title_size = 10, 
  label_column = "Desert",
  metric_display = "KO association abundance",
  plot_title = "Abundance of glycerol dehydrogenase (K00005)",
  interactive = FALSE
)

ggsave('fig_panel2.pdf')

###########################
##panel 3: lipid plot

##getprotein data
finame <- paste0('https://zenodo.org/records/13345212/files/DESI+spatial',
                 '%20proteomics%20rat%20brain%20proteomics%20',
                 'result.xlsx?download=1')
fi <- download.file(finame,
                    dest = 'dat.xlsx')

#read in data, convert to matrix 
dat <- readxl::read_xlsx('dat.xlsx') 

feature_data <- select(dat, c(Protein,`Protein ID`,`Entry Name`,Gene)) |>
  dplyr::distinct() |>
  tibble::column_to_rownames('Protein')

dat <- select(dat, c(Protein, starts_with('ROI'))) |>
  tibble::column_to_rownames('Protein')

##image data
img <- system.file("extdata","brain_img_0.png",package = 'spammR')


###now get coordinates
jsondat <- FROM_GeoJson(system.file('extdata','brain_roi.geojson',
                                    package = 'spammR'))

##herr er hry yhr
coords <- do.call(rbind, lapply(jsondat$features,function(x){
  roi <- x$properties$name
  xv <- x$geometry$coordinates[,1]
  yv <- x$geometry$coordinates[,2]
  if (is.null(roi))
    roi = ""
  return(list(ID = roi, x_pixels = min(xv), y_pixels = min(yv),
              spot_height = max(yv) - min(yv),
              spot_width = max(xv) - min(xv)))
})) |>
  as.data.frame() |>
  subset(ID != "") |>
  tidyr::separate(ID,into = c('ROI','Replicate'), 
                  sep = '_', 
                  remove = FALSE) |>
  tibble::remove_rownames() |>
  tibble::column_to_rownames('ID')

##now for each ROI we want an x, y, cell height and cell_width
#y-coordinates are from top, so need to udpate
coords$y_pixels = 263 - unlist(coords$y_pixels) - unlist(coords$spot_height)

spe <- spammR::convert_to_spe(dat = dat, 
                              feature_meta = feature_data,
                              sample_meta = coords,
                              spatial_coords_colnames = c('x_pixels','y_pixels'),
                              assay_name = 'proteomics',
                              sample_id = 'rat_brain', 
                              image_files = img,
                              image_id = 'rat_brain')

#lastly get metaspace data
mspe <- spammR::retrieve_metaspace_data("2024-02-15_20h37m13s", 
                                        fdr = 0.2, 
                                        assay_name = 'lipids',
                                        sample_id = 'rat_brain',
                                        rotate = TRUE, 
                                        drop_zeroes = TRUE)
## add in image to check 
mspe <- SpatialExperiment::addImg(mspe, img, scaleFactor = NA_real_,
                                  sample_id = 'rat_brain',
                                  image_id = 'rat_brain')

meta <- rowData(mspe) 
meta$Name <- vapply(meta$moleculeNames, function(x) {
  x[1]}, character(1))

rowData(mspe) <- meta

reduced <- spat_reduce(spe_target = spe, 
                       spe_origin = mspe, 
                       origin_assay = 'lipids')

## pick a variable lipid
lvs <- sort(apply(assay(mspe,'lipids'),1,var,na.rm = T),decreasing = TRUE)[1:10]

lvn <- rowData(altExp(reduced))[names(lvs),'Name']

lvn[1]

#find correlating proteins
fl <- c(lvn[1], rowData(reduced)[,'Gene'])


cg <- spatial_network(reduced, 
                      assay_names = c("proteomics","lipids"), 
                      query_features = lvn[1], #lipid
                      target_features = rowData(reduced)[,'Gene'],
                      feature_names = c('Gene','Name'))

##we can then use trick above to get the neighborhood of our lipid
gt <- cg %>%
  activate(nodes) |>
  as_tibble()

neigh_weights <-  cg |>
  activate(edges) |>
  filter(from == which(gt$name == lvn[1])) 

#then we can get the edge weights and nodes 
edges <- neigh_weights |> as_tibble()

nodes <- neigh_weights |>
  activate(nodes) |>
  as_tibble()


net <- neigh_weights |>
  activate(edges) |> filter(abs(corval) > 0.75) |>
  activate(nodes) |>
  mutate(degree = centrality_degree(mode = 'all')) |>
  filter(degree > 0) |>
  ggraph(layout = 'stress')  +
  geom_edge_link(aes(colour = corval)) + 
  geom_node_point(aes(color = class)) +
  geom_node_label(aes(label = name, color = class))

ggsave("fig_panel3_corNet.pdf", net)
##create a table with the protein correlationm
prot_cor <- data.frame(nodes[edges$to,],lipid_cor = edges$corval) |>
  dplyr::rename(Gene = 'name') |>
  subset(Gene %in% rowData(reduced)[,'Gene'])


# do enrichment

prot_cor <- prot_cor |> 
  group_by(Gene)|> 
  summarize(lipid_cor = mean(lipid_cor))

rowData(reduced) <- rowData(reduced) |>
  as.data.frame() |>
  left_join(prot_cor) |>
  mutate(upperGene = toupper(Gene))

data("krbpaths")
enrich <- leapR::leapR(geneset = krbpaths, enrichment_method = 'enrichment_in_order',
                       eset = reduced, primary_column = 'lipid_cor', 
                       id_column = 'upperGene')

sig <- subset(enrich,BH_pvalue < 0.3)


##panel 5: protein expression alongside lipid
p0 <- spatial_heatmap(mspe,feature = lvn[1], assay_name = 'lipids', 
                feature_type = 'Name',
                image_id = 'rat_brain', sample_id = 'rat_brain', 
                metric_display = 'Lipid expression',
                plot_log = TRUE)

ggsave('fig1_panel3_lipFul.pdf',p0)

p1 <- spatial_heatmap(reduced,feature = lvn[1], assay_name = 'lipids', 
                feature_type = 'Name',
                image_id = 'rat_brain', sample_id = 'rat_brain', 
                metric_display = 'Lipid expression',
                plot_log = TRUE)

ggsave('fig1_panel4_lipLim.pdf',p1)
# plot pathway
prots <- unlist(stringr::str_split(sig[1,'ingroupnames'], pattern = ', '))


p2 <- spatial_heatmap(reduced, assay_name = 'proteomics',
                      feature = prots,
                      feature_type = 'upperGene',
                      sample_id = 'rat_brain',
                      image_id = 'rat_brain', 
                      plot_log = TRUE, 
                      title_size = 10,
                      plot_title = gsub('_',' ',rownames(sig)[1]),
                      metric_display = 'Protein expression')
ggsave('fig1_panel5_corFul.pdf',p2)