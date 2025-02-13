## ----load package, message=FALSE,warning=FALSE--------------------------------
##load spammR
library(spammR)


## ----omics data---------------------------------------------------------------
data(pancData)
utils::head(pancData[,1:8])

## ----data list----------------------------------------------------------------
data(pancDataList)
print(length(pancDataList))
head(pancDataList[[2]][,1:8])

## ----samp meta----------------------------------------------------------------
data(pancMeta)
head(pancMeta)

## ----img----------------------------------------------------------------------
library(cowplot)

cowplot::ggdraw()+cowplot::draw_image('../inst/extdata/Image_1.png')

## ----omics meta---------------------------------------------------------------
data(protMeta)
head(protMeta[,c('pancProts','EntryName','PrimaryGeneName')])

## ----spe combined-------------------------------------------------------------
pooled.panc.spe <- convert_to_spe(pancData,  ##pooled data table
                                  pancMeta,  ##pooled metadata
                                  protMeta,  ##protein identifiers
                                  feature_meta_colname='pancProts', #column name
                                  samples_common_identifier='')
print(pooled.panc.spe)


## ----spe list-----------------------------------------------------------------

##list of image names
imglist=c('Image_0','Image_1','Image_2','Image_3','Image_4','Image_7','Image_10')

img.spes=lapply(imglist,
              function(x){
                convert_to_spe(pancDataList[[x]],
                               pancMeta,
                               protMeta,
                              feature_meta_colname='pancProts',
                              spatialCoords_colnames=c('x_pixels','y_pixels'),
                              image_files=system.file('extdata',
                                                      paste0(x,'.png'),
                                                      package='spammR'),
                              image_samples_common_identifier=x,
                              samples_common_identifier = x,
                              image_ids='with_grid')
              })
names(img.spes)<- imglist

print(img.spes)

## ----spatial data-------------------------------------------------------------

allimgs = lapply(imglist,function(x){
  spe = img.spes[[x]]
  res = spatial_heatmap(spe, feature = 'INS', sample_id=x, 
                          image_id='with_grid',spatial_coord_names=c('x_pixels','y_pixels'), 
                        spot_size=unlist(colData(spe)[1,c('spot_width','spot_height')]), 
                        image_boundaries=unlist(colData(spe)[1,c('x_origin','y_origin','x_max','y_max')]),
                        label_column='IsletOrNot', interactive=FALSE)
  return(res)}
)
cowplot::plot_grid(allimgs[[2]],allimgs[[3]])


## ----do ploting---------------------------------------------------------------

imgs = lapply(imglist,function(x){
    spatial_heatmap(img.spes[[x]],
                    feature='INS',
                    sample_id=x,
                    image_id='with_grid',
                    spatial_coord_names=c('x_pixels','y_pixels'),
                    spot_size=unlist(colData(img.spes[[x]])[1,c('spot_width','spot_height')]),
                    image_boundaries=unlist(colData(img.spes[[x]])[1,c('x_origin','y_origin','x_max','y_max')]),
                    label_column='IsletOrNot',
                    interactive=FALSE)
  })
    
imgs[[5]]

## ----plot pathway-------------------------------------------------------------



## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

