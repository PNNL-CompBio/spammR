#' Create a `SpatialExperiment` object from data
#' @description `convert_to_spe()` Puts omics data (omics measurements, metadata and image(s) corresponding to samples' tissue) into a SpatialExperiment (SPE) object. Most spammR functions require the input data to be an SPE object.
#' @export
#' @import SummarizedExperiment
#' @import SpatialExperiment
#' @param dat Matrix or data frame with omics measurements. Rows or `feature_colname` are features, columns are samples
#' @param sample_meta Data frame of metadata with either `sample_colname` or rownames as samples (if empty)
#' @param feature_meta Data frame of metadata with either `sample_colname` or rownames as samples (if empty)
#' @param assay_name Name to be given to the data in omics_measurements_file. Example: "abundance", "log2", "znormalized_log2" or any other descriptive name
#' @param sample_colname Column name in `sample_meta` table whose entries contain sample identifiers provided in `dat`.
#' @param remove_samples List of names of samples (as they occur in sample_colname) that should be removed/excluded from the data for SPE. Don't need to specify if no samples need to be removed.
#' @param feature_data_colname Name of column in `dat`, that is to be used for identifying features. If missing defaults to rownames.
#' @param feature_meta_colname Name of column in `feature_meta`, that is to be used for identifying features. If missing defaults to rownames.
#' @param spatialCoords_colnames A list containing names of columns in `meta_dat` that are spatial coordinates. Default value is NA which indicates that these columns are not present in the metadata_file, and in that case, the argument 'spatialCoords_file' must be specified.
#' @param samples_common_identifier A string (if same for all samples in omics_measurements_file) or a character vector (same length as number of samples in omics_measurements_file) corresponding to a descriptive name for samples in the current dataset. Examples: "Image0", "Experiment1", etc.
#' @param image_files A list containing paths of image files to be stored in the SpatialExperiment object. More images can be added later, without using this function.
#' @param image_ids  A list containing image names/identifiers for image paths provided in image_files
#' @param image_samples_common_identifier A list containing names of samples_common_identifier(s) corresponding to images provided in image_files. This identifier links specific samples to a experiment/condition represented by a given image.
#' @returns spe.out a SpatialExperiment (SPE) object that contains all data and image(s). Ready to be used as input in spammR functions that require SPE object as input.

#' @examples
#' data(pancData)
#' data(pancMeta)
#' data(protMeta)
#' #We can put all samples into the same object (for statistical power)
#' pooled.panc.spe <- convert_to_spe(pancData,pancMeta,protMeta,feature_meta_colname='pancProts',samples_common_identifier='')
#' #or we can add the inmage to a single data capture
#' data(pancDataList)
#' img0.spe<-convert_to_spe(pancDataList$Image_0,pancMeta,protMeta,feature_meta_colname='pancProts',image_files=system.file("extdata",'Image_0.png',package='spammR'),image_samples_common_identifier='Image0',samples_common_identifier = 'Image0',image_ids='Image0')
#'

convert_to_spe <-function(dat, ##expression data frame - rows are feature,s columns are samples
                          sample_meta, ##table of metadata for samples. rownames are columns of unless sample_colname is set
                          feature_meta, #table of metadata, row names match the rownmes of `dat` unless protein_colname is set
                          assay_name = 'proteomics',
                          sample_colname = NULL, #column with sample identifiers in `sample_meta` if not rownames
                          remove_samples = NULL, #list of samples to remove
                          feature_data_colname = NULL, ##colname of features in `dat`, if not rowname of data matrix
                          feature_meta_colname = NULL, #column with protein identifeirs in `feature_meta` if not rownames
                          spatialCoords_colnames = c('x_pixels','y_pixels'),
                          samples_common_identifier = 'sample',
                          image_files = NULL,
                          image_ids = NULL,
                          image_samples_common_identifier = NULL){

  ##first clean up samples: make sure rowanmes of metadata file are samples
  # Separate sample columns and feature meta data columns in dat
  if(!is.null(sample_colname)){
    rownames(sample_meta)<-sample_meta[[sample_colname]]
  }

  ##remove problematic samples - do we need this?
  if(!is.null(remove_samples)){
    remove_sample_colnums = which(colnames(dat) %in% remove_samples)
    dat = dat[,-remove_sample_colnums]
    sample_meta = sample_meta[!rownames(sample_meta) %in% remove_samples,]
  }

  sample_colnames <- intersect(colnames(dat),rownames(sample_meta))
  if(length(sample_colnames)==0){
    print("ERROR: rownames of sample metadata must match column names of data")
    exit()
  }

  dat_samples_only = dat[,sample_colnames]

  # Keep rows in meta data that have a corresponding sample ID in the omics measurements file
  sample_meta = sample_meta[sample_colnames,] # To be specified as colData for SPE

  ##try to collect additional data
  other_dat<-setdiff(colnames(dat),rownames(sample_meta))#sample_colnums)] # to be specified as rowData for SPE
  ##second clean up the protein metadata and make sure we have proper rowData for the spatial object
  feature_supp <- NULL
  if(length(other_dat)>0){ ##first check to see if there is other data in `dat`
    feature_supp = unique(dat[,other_dat])|>
      dplyr::mutate(proteins=rownames(dat))
  }

  if(missing(feature_meta))
    feature_meta<-feature_supp

  ##if there is a feature_meta table
  if(!is.null(feature_meta_colname) && !missing(feature_meta)){
    feature_meta<-feature_meta[!is.na(feature_meta[[feature_meta_colname]]),]
    rownames(feature_meta)<-feature_meta[[feature_meta_colname]]
  }

  if(!is.null(feature_supp)){
    feature_meta<-merge(feature_meta,feature_supp)
  }

  features <-intersect(rownames(feature_meta),rownames(dat))
  if(length(features)<nrow(dat)){
    message(paste("Mapping metadata for",length(features),'features out of',nrow(dat),'data points'))
  }

  feature_meta<-feature_meta[features,]
  dat_samples_only<-dat_samples_only[features,]

  spatialCoords_dat = as.matrix(apply(sample_meta[,spatialCoords_colnames],2,as.numeric))
  rownames(spatialCoords_dat) <- rownames(sample_meta)

  spe.out <-SpatialExperiment::SpatialExperiment(assays=list(as.matrix(dat_samples_only)),
                              colData=sample_meta,
                              rowData=feature_meta,
                              spatialCoords = spatialCoords_dat,
                              sample_id = samples_common_identifier)
  SummarizedExperiment::assayNames(spe.out) = c(assay_name)
  # Add image(s) to SPE
  if (!is.null(image_files)){
    if (is.null(image_ids)){
      # Throw error: Error: image_ids must be specified
    }
    if (!is.null(image_samples_common_identifier)){
      # Throw error: Error: image_samples_common_identifier must be specified
    }
    for (i in 1:length(image_files)){
      # scaleFactor under addImg() is a "single numeric scale factor used to rescale spatial coordinates according to the image's resolution."
      # Can turn scaleFactor into a parameter specified by user also, but not needed right now.
      spe.out = SpatialExperiment::addImg(spe.out,
                       sample_id = image_samples_common_identifier[i],
                       image_id = image_ids[i],
                       imageSource = image_files[i],
                       scaleFactor = NA_real_,
                       load = TRUE)
    }
  }
  return(spe.out)
}
