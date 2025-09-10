############ Adapted functions from Cardinal for binning matrix objects ############

#' Spectral binning of intensity values stored in a Matrix object, converted from matter.
#'
#' @param matrix matter matrix object containing the intensity values to be binned.
#' @param ref A vector of reference mass-to-charge (m/z) values for binning. If left unspecified, mass range will be used to generate reference peaks.
#' @param index Character string specifying the column name for the m/z values (default = "mz").
#' @param method A character vector specifying the binning method. Options include `"sum"`, `"mean"`, `"max"`, `"min"`. If not specified default method used is "sum".
#' @param tolerance Numeric value specifying the tolerance for m/z matching (default = `NA`).
#'
#' @return Matrix object containing the binned intensity values matching the provided reference list.
#' @export
#'
#' @examples
#' #Helper function for binning data in Matrix format
spectral_binning <- function(matrix, ref, index, method = c("sum", "mean", "max", "min"), tolerance) {
  # Ensure method is valid
  method <- match.arg(method)

  # Validate inputs
  if (!is.matrix(matrix)) {
    stop("'matrix' must be a numeric matrix")
  }
  if (is.null(ref) || length(ref) == 0) {
    stop("'ref' (reference bins) must be provided")
  }
  if (is.null(tolerance) || tolerance <= 0) {
    stop("'tolerance' must be a positive number")
  }

  # Initialize binned matrix
  binned_matrix <- matrix(0, nrow = length(ref), ncol = ncol(matrix))
  #return(binned_matrix)
  # Perform binning

  if(names(tolerance) == "relative"){
    for (i in seq_along(ref)) {
      # Calculate relative tolerance bounds
      lower_bound <- ref[i] - (ref[i] * tolerance)
      upper_bound <- ref[i] + (ref[i] * tolerance)

      # Identify rows (spectral features) within the current bin
      in_bin <- which(index >= lower_bound & index <= upper_bound)

      if (length(in_bin) > 0) {
        # Apply the chosen method to bin the data
        binned_matrix[i, ] <- switch(method,
                                     sum = colSums(matrix[in_bin, , drop = FALSE]),
                                     mean = colMeans(matrix[in_bin, , drop = FALSE]),
                                     max = apply(matrix[in_bin, , drop = FALSE], 2, max),
                                     min = apply(matrix[in_bin, , drop = FALSE], 2, min))
      }
    }
  } else{

    for (i in seq_along(ref)) {

      # Define the bin range for the current reference value
      lower_bound <- ref[i] - tolerance
      upper_bound <- ref[i] + tolerance

      # Identify rows (spectral features) within the current bin
      in_bin <- which(index >= lower_bound & index <= upper_bound)

      if (length(in_bin) > 0) {
        # Apply the chosen method to bin the data
        binned_matrix[i, ] <- switch(method,
                                     sum = colSums(matrix[in_bin, , drop = FALSE]),
                                     mean = colMeans(matrix[in_bin, , drop = FALSE]),
                                     max = apply(matrix[in_bin, , drop = FALSE], 2, max),
                                     min = apply(matrix[in_bin, , drop = FALSE], 2, min))
      }
    }
  }

  # Set row names to the reference values
  rownames(binned_matrix) <- as.character(ref)

  return(binned_matrix)
}




#' Bin SpaMTP Object
#'
#' Bins m/z peaks of a specified SpaMTP object to reduce dimensionality and remove noise.
#'
#' @param data SpaMTP Seurat Object containing the intensity data to be binned
#' @param resolution Numeric value indicating the relative bin size to use when binning the intensity data.
#' @param units Character string defining the relative units of the provided resolution size. Values can be either 'ppm' or 'mz' (default = 'ppm').
#' @param assay Character string defining the SpaMTP assay to extract intensity values from (default = "Spatial").
#' @param slot Character string defining the assay slot containing the intensity values (default = "counts").
#' @param method A character vector specifying the binning method. Input values must be either `"sum"`, `"mean"`, `"max"` or `"min"` (default = "sum").
#' @param return.only.mtx Boolean value indicating whether to return only the binned intensity matrix. If `FALSE`, a SpaMTP object will be returned with the binned values stored in an assay called `binned` (default = FALSE).
#'
#' @return Either a binned intensity matrix or a SpaMTP object contatining the binned intensity matrix stored in the `binned` assay.
#' @export
#'
#' @examples
#' #BinSpaMTP(spamtp.obj, resolution = 10, units = "ppm", return.only.mtx = TRUE)
BinSpaMTP <- function(data, resolution, units = "ppm", assay = "Spatial",slot = "counts", method = c("sum"), return.only.mtx = FALSE){

  orignal_bin_size <- resolution
  min_ref <- min(data[[assay]]@meta.data$raw_mz)
  max_ref <- max(data[[assay]]@meta.data$raw_mz)

  if (units == "ppm"){
    resolution <- 1e-6 * resolution
    bin_class <- "relative"
  } else if ( units == "mz"){
    bin_class <- "absolute"
  } else {
    stop("Incorrect units provided! units must be either 'ppm' or 'mz'. Please use either option for binning")
  }


  if(bin_class == "relative"){
    ref <- seq_rel(min_ref, max_ref, by=resolution)
  } else {
    ref <- seq.default(min_ref, max_ref, by=resolution)
  }

  if ( method %in% c("sum", "mean", "max", "min") ) {
    tol <- 0.5 * resolution
  } else {
    tol <- 2 * resolution
  }

  names(tol) <- bin_class
  index <- data[[assay]]@meta.data$raw_mz

  if (length(ref) < length(index)){
    mtx <- spectral_binning(matrix=as.matrix(data[[assay]][slot]), ref = ref, index = index, method = method, tolerance = tol)
    colnames(mtx) <- rownames(data@meta.data)
  } else {
    warning("Bin size is too small to bin m/z values together. Currently, with a bin size of ",
            orignal_bin_size, " ", units, " No m/z values will be binned together ... The original intensity matrix will be used!")
    mtx <- data[[assay]][slot]
  }


  if (return.only.mtx){
    return(mtx)
  } else {
    obj <- SeuratObject::CreateAssayObject(counts = mtx)
    data[["binned"]] <- obj
    return(data)
  }
}




