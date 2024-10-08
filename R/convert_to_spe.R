#' convert_to_spe: Puts omics data (omics measurements, metadata and image(s) corresponding to samples' tissue) into a SpatialExperiment (SPE) object. Most spammR functions require the input data to be an SPE object.

#' @export
#' @param omics_measurements_file File path for omics measurements. Rows correspond to features. Samples are in columns. Can have additional columns also describing features.
#' @param assay_name Name to be given to the data in omics_measurements_file. Example: "abundance", "log2", "znormalized_log2" or any other descriptive name
#' @param metadata_file File path for metadata for samples. This can include spatial coordinates.
#' @param meta_colname_sampleIDs Column name in metadata_file whose entries contain sample identifiers provided in omics_measurements_file.
#' @param remove_samples Names of samples (as they occur in meta_colname_sampleIDs) that should be removed/excluded from the data for SPE. Don't need to specify if no samples need to be removed.
#' @param feature_colname Name of column in omics_measurements_file, that is to be used for identifying features.
#' @param spatialCoords_colnames A list containing names of columns in metadata_file that are spatial coordinates. Default value is NA which indicates that these columns are not present in the metadata_file, and in that case, the argument 'spatialCoords_file' must be specified.
#' @param spatialCoords_file If spatial coordinates are not provided in the metadata_file, then a file path for the spatial coordinates file should be specified. This file should contain only columns corresponding to spatial coordinates. Rows should represent samples.
#' @param samples_common_identifier A string (if same for all samples in omics_measurements_file) or a character vector (same length as number of samples in omics_measurements_file) corresponding to a descriptive name for samples in the current dataset. Examples: "Image0", "Experiment1", etc.
#' @param image_files A list containing paths of image files to be stored in the SpatialExperiment object. More images can be added later, without using this function.
#' @param image_ids  A list containing image names/identifiers for image paths provided in image_files
#' @param image_samples_common_identifier A list containing names of samples_common_identifier(s) corresponding to images provided in image_files. This identifier links specific samples to a experiment/condition represented by a given image.
#' @returns spe.out a SpatialExperiment (SPE) object that contains all data and image(s). Ready to be used as input in spammR functions that require SPE object as input.

# Example input parameters
# omics_measurements_file1 = "/Users/sohi472/Library/CloudStorage/OneDrive-PNNL/Projects/BICCN/data/brain_14ROIs_data_Aug24_2023/useThis/image0_data_only/image0_prot_dat_shortSampleNames.xlsx"
# assay_name1 = "abundance" # Name to be given to the current data.
# metadata_file1 = "/Users/sohi472/Library/CloudStorage/OneDrive-PNNL/Projects/BICCN/data/brain_14ROIs_data_Aug24_2023/useThis/image0_data_only/image0_samples_meta_data.xlsx"
# meta_colname_sampleIDs1 = "Sample_ID"
# remove_samples1 = c("0_TH_3","0_TEP_3","0_TH_2")
# feature_colname1 = "PG.ProteinNames"
# spatialCoords_colnames1 = c("Xcoord","Ycoord")
# samples_common_identifier1 = "Image0" # Name of a common identifier for samples
# img1 = "/Users/sohi472/Library/CloudStorage/OneDrive-PNNL/Projects/BICCN/data/brain_14ROIs_data_Aug24_2023/tissue_images_forPlotting/Image0_Raw_noMarkings_cropped_forPlotting.png"
# img2 = "/Users/sohi472/Library/CloudStorage/OneDrive-PNNL/Projects/BICCN/data/brain_14ROIs_data_Aug24_2023/tissue_images_forPlotting/Image0_ROIsMarked_cropped_forPlotting.png"
# image_files1 = c(img1,img2)
# image_ids1 = c("Raw_noMarkings","ROIsMarked") # Image names/identifiers for image paths provided in image_files
# image_samples_common_identifier1 = c("Image0","Image0") #Name of a common identifier that links specific samples to a experiment/condition represented by a given image.

convert_to_spe <-function(omics_measurements_file, assay_name, metadata_file, meta_colname_sampleIDs, remove_samples=NULL, feature_colname, spatialCoords_colnames, spatialCoords_file=NULL, samples_common_identifier, image_files=NULL, image_ids=NULL, image_samples_common_identifier=NULL){
  dat = data.frame(read_excel(path=omics_measurements_file),check.names = FALSE)
  meta_dat = data.frame(read_excel(path=metadata_file),check.names = FALSE)
  if(!is.null(remove_samples)){
    remove_sample_colnums = which(colnames(dat) %in% remove_samples)
    dat = dat[,-remove_sample_colnums]
    meta_dat = meta_dat[!meta_dat[,meta_colname_sampleIDs] %in% remove_samples,]
  }
  # Separate sample columns and feature meta data columns in dat
  sample_colnums = which(colnames(dat) %in% meta_dat[,meta_colname_sampleIDs])
  sample_colnames = colnames(dat)[sample_colnums]
  dat_samples_only = dat[,sample_colnames]
  features_info = dat[,-c(sample_colnums)] # to be specified as rowData for SPE
  # The list of samples specified in the data and metadata to be stored in the SPE object should be exactly the same.
  # Keep rows in meta data that have a corresponding sample ID in the omics measurements file
  meta_dat_keep = meta_dat[meta_dat[,meta_colname_sampleIDs] %in% sample_colnames,] # To be specified as colData for SPE
  rownames(meta_dat_keep) = meta_dat[,meta_colname_sampleIDs]
  if (is.null(spatialCoords_file)){
    spatialCoords_dat = as.matrix(meta_dat_keep[,spatialCoords_colnames])
  }else{
    spatialCoords_dat = as.matrix(data.frame(read_excel(spatialCoords_file),check.names = FALSE))
  }
  spe.out <-SpatialExperiment(assays=list(as.matrix(dat_samples_only)),
                              colData=meta_dat_keep,
                              rowData=features_info,
                              spatialCoords = spatialCoords_dat,
                              sample_id = samples_common_identifier)
  names(assays(spe.out)) = c(assay_name)
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
      spe.out = addImg(spe.out,
                       sample_id = image_samples_common_identifier[i],
                       image_id = image_ids[i],
                       imageSource = image_files[i],
                       scaleFactor = NA_real_,
                       load = TRUE)
    }
  }
  return(spe.out)
}
