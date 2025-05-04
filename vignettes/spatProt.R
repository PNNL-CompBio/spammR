## ----load package, message=FALSE,warning=FALSE--------------------------------
##load spammR
library(spammR)


## ----omics data---------------------------------------------------------------
data(pancDataList)
utils::head(pancDataList$Image_0[,1:8])

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

pooledData<-dplyr::bind_cols(pancDataList)
pooled.panc.spe <- convert_to_spe(pooledData,  ##pooled data table
                                  pancMeta,  ##pooled metadata
                                  protMeta,  ##protein identifiers
                                  feature_meta_colname='pancProts', #column name
                                  samples_common_identifier='')
print(pooled.panc.spe)


## ----spe list, echo=FALSE, warning=FALSE, message=FALSE-----------------------

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

#print(img.spes)

## ----spatial data-------------------------------------------------------------

allimgs = lapply(imglist,function(x){
  spe = img.spes[[x]]
  res = spatial_heatmap(spe, feature = 'INS',
                        feature_type='PrimaryGeneName',
                        sample_id=x, 
                        image_id='with_grid',
                        spatial_coord_names=c('x_pixels','y_pixels'), 
                        spot_size=unlist(colData(spe)[1,c('spot_width',
                                                          'spot_height')]), 
                        image_boundaries=unlist(colData(spe)[1,c('x_origin',
                                                                 'y_origin',
                                                                 'x_max',
                                                                 'y_max')]),
                        label_column='IsletOrNot', interactive=FALSE)
  return(res)}
)
cowplot::plot_grid(allimgs[[2]],allimgs[[3]],nrow=2)


## ----diffex, warning=FALSE, error=FALSE, message=FALSE------------------------


islet_res <- calc_spatial_diff_ex(pooled.panc.spe,
                               assay_name='proteomics',
                               log_transformed=FALSE,
                               category_col='IsletOrNot')

sig_prots<-subset(rowData(islet_res),NonIslet_vs_Islet.adj.P.Val.limma<0.01)
ups<-subset(sig_prots, NonIslet_vs_Islet.logFC.limma>0)
downs<-subset(sig_prots, NonIslet_vs_Islet.logFC.limma<0)

print(paste('We found',nrow(sig_prots),'significantly differentally expressed proteins including',
            nrow(ups),'upregulated proteins and',nrow(downs),'downregulated'))




