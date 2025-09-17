#### SpaMTP Cardinal SSC clustering segmentation ###############################################################################################################################################################################


#' Adds Cardinal ssc segmentation annotation to m/z count data object
#'
#' Add Spatial Shrunked Centroid (ssc) results from a non-intensity based Cardinal Object into the pixelData slot of a provided Cardinal Object containing intensity values for each m/z (e.g. MSImagingExperiment).
#' This function is additional functionality that can only be run on `Cardinal` Objects ONLY.
#'
#' @param data A Cardinal Object containing the raw/binned m/z count data.
#' @param data_ssc A Cardinal Object containing the ssc segmentation results. Note: Cardinal's spatialShrunkenCentroids() must be run to generate this object.
#' @param resolution An integer (< v3.6.0) or string (>= v3.6.0) defining the ssc segmentation resolution to add to the raw Cardinal object.
#'
#' @returns A Cardinal Object containing the raw m/z counts and the relative ssc segmentation annotations per pixel.
#' @export
#'
#' @examples
#' # ssc_data <- Cardinal::spatialShrunkenCentroids(CardinalObj, ...)
#' # new_CardinalObj <- add_ssc_annotation(CardinalObj, ssc_data, resolution ="r=2,k=8,s=32")
add_ssc_annotation <- function(data, data_ssc, resolution){

  data_bin <- data

  if (check_cardinal_version()){

    if(class(data_ssc) != 'ResultsList'){
      Cardinal::pixelData(data_bin)[["ssc"]] <- data_ssc$class
    } else {
      message(paste0("Getting cluster segments for resolution ", resolution))

      Cardinal::pixelData(data_bin)[["ssc"]] <- data_ssc@listData[[resolution]]$class
    }

  } else {
    message(paste0("Getting cluster segments for resolution (s) = ", resolution))
    cluster_idx <- which(Cardinal::modelData(data_ssc)[["s"]] == resolution)

    classes <- Cardinal::resultData(data_ssc)[[cluster_idx]][[1]]
    pixel_data <- Cardinal::pixelData(data_bin)


    pixel_data[["ssc"]] <- classes

    Cardinal::pixelData(data_bin) <- pixel_data

  }

  return(data_bin)

}

########################################################################################################################################################################################################################

