## ----setup, include=FALSE, echo = FALSE---------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(spammR)

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

## ----get proteomics data, warning=FALSE, message=FALSE------------------------

#brain data
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


## ----load image, warning=FALSE, message=FALSE---------------------------------
##read in image
img <- system.file("extdata","brain_img_0.png",package = 'spammR')


## ----plot image, eval=FALSE---------------------------------------------------
# 
# cowplot::ggdraw() + cowplot::draw_image(system.file("extdata",
#                                                     "brain_img_0.png",
#                                                     package = "spammR"))

## ----format coordinates, warning=FALSE, message=FALSE-------------------------
#remove this once we can include it in spammR
 
  library(geojsonR)
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
 

## ----create proteomics spe, warning=FALSE, message=FALSE----------------------
   #create an SFE
  spe <- spammR::convert_to_spe(dat = dat, 
                      feature_meta = feature_data,
                      sample_meta = coords,
                      spatial_coords_colnames = c('x_pixels','y_pixels'),
                      assay_name = 'proteomics',
                      sample_id = 'rat_brain', 
                      image_files = img,
                      image_id = 'rat_brain')

  # myelin
  spammR::spatial_heatmap(spe, feature = 'Mbp', 
                          feature_type = 'Gene',
                sample_id = 'rat_brain', 
                image_id = 'rat_brain')


## ----correlation enrichment, warning=FALSE, message=FALSE---------------------

library(leapR)
data(krbpaths)

#first we create human gene names from rat
rowData(spe) <- rowData(spe) |> 
  as.data.frame() |>
  dplyr::mutate(upperGene = toupper(Gene))

#then calculate correlation enrichment
res <- leapR::leapR(krbpaths, 'correlation_enrichment',id_column = 'upperGene',
                    spe,'proteomics')

res[,c('ingroup_n','ingroup_mean','pvalue','BH_pvalue')] |> 
  dplyr::arrange(BH_pvalue) |>
  head()


## ----heatmaps, warning=FALSE, message=FALSE-----------------------------------

oxphos <- unlist(strsplit(res$ingroupnames[1],','))

#lets visualize only the oxphos proteins in the brain
spatial_heatmap(spe,feature = oxphos, assay_name = 'proteomics',
                sample_id = 'rat_brain',  image_id = 'rat_brain',
                plot_log = TRUE)


## ----protein variance, warning=FALSE, message=FALSE---------------------------

var_val <- apply(assay(spe,'proteomics')[oxphos,], 1, var) |>
  sort()

vv <- data.frame(Gene = rowData(spe)[names(var_val),'Gene'],
                 var = var_val)

ggplot(vv,aes(x = reorder(Gene,var), y = var)) + 
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


## ----correlation networks, warning=FALSE, message=FALSE-----------------------

#let's get gene names
gn <- rowData(spe)[oxphos,'Gene']

pg <- spatial_network(spe, features_of_interest = gn,
                      assay_names = 'proteomics',
                      feature_names = 'Gene')

##let's get a list of some interesting proteins??
pg |> activate(edges) |>
  filter(corval > 0.5) |>
  ggraph()  +
  geom_edge_link(aes(colour = corval)) + 
  geom_node_point() + 
  geom_node_label(aes(label = name))


## ----hk1 plot-----------------------------------------------------------------

spatial_heatmap(spe,feature = 'Hk1', feature_type = 'Gene',
                assay_name = 'proteomics',
                sample_id = 'rat_brain',  image_id = 'rat_brain',
                plot_log = TRUE)


spatial_heatmap(spe,feature = 'Eno1', feature_type = 'Gene',
                assay_name = 'proteomics',
                sample_id = 'rat_brain',  image_id = 'rat_brain',
                plot_log = TRUE)

