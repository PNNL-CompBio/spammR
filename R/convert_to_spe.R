#' convert_to_spe: Puts omics data (omics measurements, metadata and image(s) corresponding to samples' tissue) into a SpatialExperiment (SPE) object. Most spammR functions require the input data to be an SPE object.

#' @export
#' @import SpatialExperiment
#' @param dat Matrix or data frame with omics measurements. Rows or `feature_colname` are features, columns are samples
#' @param meta_dat Data frame of metadata with either `sample_colname` or rownames as samples (if empty)
#' @param assay_name Name to be given to the data in omics_measurements_file. Example: "abundance", "log2", "znormalized_log2" or any other descriptive name
#' @param sample_colname Column name in metadata_file whose entries contain sample identifiers provided in omics_measurements_file.
#' @param remove_samples List of names of samples (as they occur in sample_colname) that should be removed/excluded from the data for SPE. Don't need to specify if no samples need to be removed.
#' @param feature_colname Name of column in omics_measurements_file, that is to be used for identifying features.
#' @param spatialCoords_colnames A list containing names of columns in `meta_dat` that are spatial coordinates. Default value is NA which indicates that these columns are not present in the metadata_file, and in that case, the argument 'spatialCoords_file' must be specified.
#' @param samples_common_identifier A string (if same for all samples in omics_measurements_file) or a character vector (same length as number of samples in omics_measurements_file) corresponding to a descriptive name for samples in the current dataset. Examples: "Image0", "Experiment1", etc.
#' @param image_files A list containing paths of image files to be stored in the SpatialExperiment object. More images can be added later, without using this function.
#' @param image_ids  A list containing image names/identifiers for image paths provided in image_files
#' @param image_samples_common_identifier A list containing names of samples_common_identifier(s) corresponding to images provided in image_files. This identifier links specific samples to a experiment/condition represented by a given image.
#' @returns spe.out a SpatialExperiment (SPE) object that contains all data and image(s). Ready to be used as input in spammR functions that require SPE object as input.

#' @examples
#' data(pancData)
#' data(pancMeta)
#' allImages <- convert_to_spe(pancData,pancMeta,samples_common_identifier='')
#
# Example input parameters (remove this once we have it all in examples)
# samples_common_identifier1 = "Image0" # Name of a common identifier for samples
# img1 = "/Users/sohi472/Library/CloudStorage/OneDrive-PNNL/Projects/BICCN/data/brain_14ROIs_data_Aug24_2023/tissue_images_forPlotting/Image0_Raw_noMarkings_cropped_forPlotting.png"
# img2 = "/Users/sohi472/Library/CloudStorage/OneDrive-PNNL/Projects/BICCN/data/brain_14ROIs_data_Aug24_2023/tissue_images_forPlotting/Image0_ROIsMarked_cropped_forPlotting.png"
# image_files1 = c(img1,img2)
# image_ids1 = c("Raw_noMarkings","ROIsMarked") # Image names/identifiers for image paths provided in image_files
# image_samples_common_identifier1 = c("Image0","Image0") #Name of a common identifier that links specific samples to a experiment/condition represented by a given image.

convert_to_spe <-function(dat, ##expression data frame - rows are feature,s columns are samples
                          meta_dat, ##table of metadata. rownames are columns of unless sample_colname is set
                          assay_name = 'proteomics',
                          sample_colname=NULL, #column with sample identifiers if not rownames
                          remove_samples=NULL, #list of samples to remove
                          feature_colname = NULL, ##colname of features, if not rowname of data matrix
                          spatialCoords_colnames=c('Xcoord','Ycoord'),
                          samples_common_identifier='sample',
                          image_files=NULL,
                          image_ids=NULL,
                          image_samples_common_identifier=NULL){
  #ibrary("SpatialExperiment")

  # Separate sample columns and feature meta data columns in dat
  if(!is.null(sample_colname)){
    meta_dat <- meta_dat
    rownames(meta_dat)<-meta_dat[[sample_colname]]
  }

  ##remove problematic samples - do we need this?
  if(!is.null(remove_samples)){
    remove_sample_colnums = which(colnames(dat) %in% remove_samples)
    dat = dat[,-remove_sample_colnums]
    meta_dat = meta_dat[!rownames(meta_dat) %in% remove_samples,]
  }

  sample_colnames <- intersect(colnames(dat),rownames(meta_dat))
  dat_samples_only = dat[,sample_colnames]

  other_dat<-setdiff(colnames(dat),rownames(meta_dat))#sample_colnums)] # to be specified as rowData for SPE

  if(length(other_dat)>0){
    features_info = dat[,other_dat]|>
      dplyr::mutate(proteins=rownames(dat))
  }else{
    features_info = data.frame(proteins=rownames(dat))
    rownames(features_info)<-rownames(dat)
    }
  # The list of samples specified in the data and metadata to be stored in the SPE object should be exactly the same.
  # Keep rows in meta data that have a corresponding sample ID in the omics measurements file
  meta_dat_keep = meta_dat[sample_colnames,] # To be specified as colData for SPE
#    rownames(meta_dat_keep) = meta_dat[,sample_colname]

  spatialCoords_dat = as.matrix(apply(meta_dat_keep[,spatialCoords_colnames],2,as.numeric))

  spe.out <-SpatialExperiment::SpatialExperiment(assays=list(as.matrix(dat_samples_only)),
                              colData=meta_dat_keep,
                              rowData=features_info,
                              spatialCoords = spatialCoords_dat,
                              sample_id = samples_common_identifier)
  SpatialExperiment::assayNames(spe.out) = c(assay_name)
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
