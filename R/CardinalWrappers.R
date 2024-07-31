#library(Cardinal)

#### SpaMTP Cardinal SSC clustering segmentation ###############################################################################################################################################################################


#' Adds Cardinal ssc segmentation annotation to m/z count data object
#'    - When ssc is run an object is returned which cannot plot raw m/z values
#'
#' @param data A Cardinal Object containing the raw/binned m/z count data.
#' @param data_ssc A Cardinal Object containing the ssc segmentation results. Note: Cardinal's spatialShrunkenCentroids() must be run to generate this object.
#' @param resolution An integer defining the ssc segmentation resolution to add to the raw Cardinal object (default = 25).
#'
#' @returns A Cardinal Object containing the raw m/z counts and the relative ssc segmentation annotations per pixel.
#' @export
#'
#' @examples
#' # ssc_data <- spatialShrunkenCentroids(CardinalObj, ...)
#' # new_CardinalObj <- add_ssc_annotation(CardinalObj, ssc_data, resolution =25)
add_ssc_annotation <- function(data, data_ssc, resolution = 25){

  message(paste0("Getting cluster segments for resolution (s) = ", resolution))
  data_bin <- data
  cluster_idx <- which(Cardinal::modelData(data_ssc)[["s"]] == resolution)

  classes <- Cardinal::resultData(data_ssc)[[cluster_idx]][[1]]
  pixel_data <- Cardinal::pixelData(data_bin)


  pixel_data[["ssc"]] <- classes

  Cardinal::pixelData(data_bin) <- pixel_data
  return(data_bin)

}

########################################################################################################################################################################################################################

