#' Loads in spatial metabolomic data directly to a SpaMTP Seurat Object
#'
#' @param name Character string of the object name. This should match the filename.
#' @param path Character string defining the directory path of the file. This should not include the file name.
#' @param mass.range Vector of numeric values indicating the mass range to use for the imported data (default = NULL).
#' @param resolution Numeric value defining the the accuracy to which the m/z values will be binned after reading. This value can be in either "ppm" or "mz" depending on the units type specified (default = 10).
#' @param units Character string defining the resolution value unit type, either c("ppm", "mz") (default = "ppm")
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE)
#' @param assay Character string describing the name of the new assay which stores the imported data (default = "Spatial").
#' @param bin_package Character string defining the package used to bin the imported data. Options are either "SpaMTP" or "Cardinal" (default = "SpaMTP").
#' @param ... Additional arguments passed to the \code{readMSIData} function.
#'
#' @return A new SpaMTP Seurat object contain the imported spatial metabolic intensity values
#' @export
#'
#' @examples
#' # data <-loadSM(name = "run1", folder = "/Documents/SpaMTP_test_data/", mass.range = c(160,1500), resolution = 10, assay = "Spatial")
LoadSM <- function (name, path, mass.range = NULL, resolution = 10, units = "ppm", verbose = TRUE, assay = "Spatial", bin_package = "SpaMTP", ...){

  if (check_cardinal_version()){
    file_name <- paste0(path, name)
    data <- Cardinal::readImzML(file = file_name, mass.range = mass.range, resolution = resolution, units = units, verbose = verbose, ...)
    if (!is.null(mass.range)| !is.null(resolution)){
      if (bin_package == "Cardinal"){
        verbose_message(message_text = "Binning data using Cardinal's m/z bin method .... ", verbose = verbose)
        warning("If data loading/conversion is taking a long time try chaning bin_method = 'SpaMTP'... This function speads up matrix conversion for data with identical m/z values for each pixel!")
        data <- Cardinal::bin(data, mass.range = mass.range, resolution = resolution, units = units, verbose = verbose, ...)
        data <- CardinalToSeurat(data, name, verbose = verbose, assay = assay)
      } else if (bin_package == "SpaMTP"){
        verbose_message(message_text = "Binning data using SpaMTP's m/z bin method .... ", verbose = verbose)
        mtx <- bin_cardinal(data, mass.range = mass.range, resolution = resolution, units = units, ...)
        data <- BinnedCardinalToSeurat(data, mtx, verbose = verbose, assay = assay)
      } else {
        stop("bin_package value is incorrect! bin_package must be either 'SpaMTP' or 'Cardinal'")
      }
    } else {
      data <- CardinalToSeurat(data, name, verbose = verbose, assay = assay)
    }
  } else {
    if (bin_package == "Cardinal"){
      verbose_message(message_text = "Binning data using Cardinal's m/z bin method .... ", verbose = verbose)
      warning("If data loading/conversion is taking a long time try chaning bin_method = 'SpaMTP'... This function speads up matrix conversion for data with identical m/z values for each pixel!")
      data <- Cardinal::readImzML(name,folder = path, mass.range =  mass.range, resolution = resolution, ...)
      data <- CardinalToSeurat(data, name, verbose = verbose, assay = assay)
    } else if (bin_package == "SpaMTP"){
      verbose_message(message_text = "Binning data using SpaMTP's m/z bin method .... ", verbose = verbose)
      data <- Cardinal::readImzML(name,folder = path, mass.range =  NULL, resolution = NULL, ...)
      mtx <- bin_cardinal(data, units = units, mass.range = mass.range, resolution = resolution, ...)
      data <- BinnedCardinalToSeurat(data, mtx, verbose = verbose, assay = assay)

    } else {
      stop("bin_package value is incorrect! bin_package must be either 'SpaMTP' or 'Cardinal'")
    }
  }
  return(data)
}



#' Read a Spatial Metabolomics image matrix file (.csv format)
#'
#' @param mtx.file Character string defining the path of the spatial metabolomic image matrix .csv file
#' @param assay Character string of the Seurat object assay name to store the relative intensity data (default = "Spatial").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param feature.start.column Numeric value defining the start index containing the x, y and m/z value columns within the table (default = 1).
#' @param mz.prefix Character string matching the prefix string in front of each m/z name (deafult = NULL).
#' @param project.name Character string defining the name of the sample to be assigned as orig.idents (default = "SpaMTP").
#'
#' ### Details
#'   NOTE: This file must be in a format similar to the one below:
#'
#'          A data.frame: 5 × 5
#'        x	   y    mz1   mz2   mz3
#'      <int> <int> <dbl> <dbl> <dbl>
#'    1	 0	   1	   0	   0	   11
#'    2	 0	   2	   0	   0	   0
#'    3	 0	   3	   0	   0	   0
#'    4	 0	   4	   20	   0	   0
#'    5	 0	   5	   0	   0	   0
#'
#'  Where:
#'      - The first 2 columns are labeled x and y containing the respective x/y spatial coordinates
#'      - The next columns are then the respective m/z values and their intensities for each spatial pixel
#'
#' @return A SpaMTP Seurat class object containing the intensity values in the counts slot of the designated assay
#' @export
#'
#' @examples
#' # msi_data <- ReadSM_mtx("~/Documents/msi_mtx.csv")
ReadSM_mtx <- function(mtx.file, assay = "Spatial", verbose = TRUE, feature.start.column = 1, mz.prefix = NULL, project.name = "SpaMTP"){

  #verbose_message(message_text = "Convering mtx to SpaMTP Seurat object .... ", verbose = verbose)

  if (!file.exists(mtx.file)) {
    stop("Expression matrix file missing. Expecting matrix.csv")
  }

  data <- as.data.frame(data.table::fread(mtx.file))

  if (feature.start.column > 0){
    verbose_message(message_text = paste0("Spliting matrix data from column ", feature.start.column," onwards .... "), verbose = verbose)
    data <- data[,feature.start.column:dim(data)[2]]
  }

  if ("x" %in% colnames(data) && "y" %in% colnames(data)) {
    coords <- data[c("x","y")]
    coords$x_coord <- coords$x
    coords$y_coord <- coords$y
    coords <- coords[c("x_coord","y_coord")]

    data <-  data[, !(names(data) %in% c("x", "y"))]

  } else {
    stop("X and Y Tissue coordinates not found. Expecting columns 'x' and 'y'")
  }

  barcodes <- paste0(coords$x_coord, "_", coords$y_coord)

  rownames(data) <- barcodes
  rownames(coords) <- barcodes

  if(!(is.null(mz.prefix))){
    colnames(data) <- gsub(mz.prefix, "mz-", colnames(data))
  } else{
    colnames(data) <- paste0("mz-", colnames(data))
  }

  if(is.null(project.name)){
    project.name <- "SpaMTP"
  }

  seuratobj <- Seurat::CreateSeuratObject(t(data), assay = assay, project = project.name)

  verbose_message(message_text = "Adding Pixel Metadata ....", verbose = verbose)

  for (name in colnames(coords)){
    seuratobj <- Seurat::AddMetaData(seuratobj,col.name = name, metadata = coords[[name]])
  }

  verbose_message(message_text = "Creating Centroids for Spatial Seurat Object ....", verbose = verbose)

  ## Add spatial data
  cents <- SeuratObject::CreateCentroids(data.frame(x = coords[["x_coord"]], y = coords[["y_coord"]], cell = rownames(coords)))

  segmentations.data <- list(
    "centroids" = cents
  )

  coords <- SeuratObject::CreateFOV(
    coords = segmentations.data,
    type = c("centroids"),
    molecules = NULL,
    assay = assay
  )
  seuratobj[["fov"]] <- coords


  metadata <- data.frame("raw_mz" = sapply(strsplit(rownames(seuratobj), "-"), function(x) as.numeric(x[[2]])))
  rownames(metadata) <- rownames(seuratobj)


  seuratobj[[assay]] <- Seurat::AddMetaData(object = seuratobj[[assay]],
                                            metadata = metadata,
                                            col.name = 'raw_mz')


  seuratobj[[assay]]@meta.data$mz_names <- rownames(seuratobj)

  return(seuratobj)
}



#' Converts a Cardinal Object into a Seurat Object
#'
#' @param data A Cardinal Object that is being converted into a Seurat Object.
#' @param mtx Matrix object containing
#' @param run_name A character string defining the run name of the Cardinal data to be converted to a Seurat Object
#' @param seurat.coord A Data.Frame containing two columns titled 'X_new' and 'Y_new' specifying the pixel coordinates of each data point. This is only required if mapping Spatial Metabolic data with a H&E image was done externally, and the SM coordinates need to change to align correctly to the ST data. Else, set to `NULL` (default = NULL).
#' @param assay Character string containing the name of the assay (default = "Spatial").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A Seurat Object containing the mz count data of the supplied Cardinal Object
#' @export
#'
#' @examples
#' # CardinalToSeurat(CardinalObj, run_name = "run_1")
BinnedCardinalToSeurat <- function(data, mtx, run_name,  assay = "Spatial", verbose = TRUE ){

  verbose_message(message_text = "Convering Cardinal object to Seurat object .... ", verbose = verbose)
  run_data <- data
  #run_data <- Cardinal::subsetPixels(data, Cardinal::run(data) == paste0(run_name)) #subset broken under Cardinal >3.6
  sparse_matrix <- mtx
  pixel_data <- Cardinal::pixelData(run_data)

  if ("x" %in% colnames(data.frame(pixel_data)) & "y" %in% colnames(data.frame(pixel_data))){
    pixel_data_df <- data.frame(pixel_data)
    pixel_data[["x_coord",]] <- pixel_data_df$x
    pixel_data[["y_coord",]] <- pixel_data_df$y
    Cardinal::pixelData(run_data) <- pixel_data
  } else {
    warning("There is no column called 'x' and 'y' in pixelData(CardinalObject)")
    stop("x and y pixel columns do not exist")
  }

  verbose_message(message_text = "Generating Seurat Barcode Labels from Pixel Coordinates .... ", verbose = verbose)

  x_coord <- Cardinal::pixelData(run_data)[["x_coord",]]
  y_coord <- Cardinal::pixelData(run_data)[["y_coord",]]
  spot_name <- paste0(x_coord,"_",y_coord)


  colnames(sparse_matrix)<- spot_name
  rownames(sparse_matrix)<- paste0("mz-", rownames(sparse_matrix))

  verbose_message(message_text = "Constructing Seurat Object ....", verbose = verbose)


  seuratobj <- Seurat::CreateSeuratObject(sparse_matrix, assay = "Spatial")

  verbose_message(message_text = "Adding Pixel Metadata ....", verbose = verbose)

  seuratobj <- Seurat::AddMetaData(seuratobj,col.name = "sample", metadata = Cardinal::run(run_data))

  for (name in names(Cardinal::pixelData(run_data))){
    seuratobj <- Seurat::AddMetaData(seuratobj,col.name = name, metadata = Cardinal::pixelData(run_data)[[name,]])
  }

  verbose_message(message_text = "Creating Centroids for Spatial Seurat Object ....", verbose = verbose)

  ## Add spatial data
  cents <- SeuratObject::CreateCentroids(data.frame(x = c(Cardinal::pixelData(run_data)[["x_coord",]]), y = c(Cardinal::pixelData(run_data)[["y_coord",]]), cell = c(spot_name)))


  segmentations.data <- list(
    "centroids" = cents
  )

  coords <- SeuratObject::CreateFOV(
    coords = segmentations.data,
    type = c("centroids"),
    molecules = NULL,
    assay = "Spatial"
  )

  seuratobj[["fov"]] <- coords


  metadata <- data.frame("raw_mz" = sapply(strsplit(rownames(seuratobj), "-"), function(x) as.numeric(x[[2]])))
  rownames(metadata) <- rownames(seuratobj)


  seuratobj[["Spatial"]] <- Seurat::AddMetaData(object = seuratobj[["Spatial"]],
                                                metadata = metadata,
                                                col.name = 'raw_mz')


  seuratobj[["Spatial"]]@meta.data$mz_names <- rownames(seuratobj)

  return(seuratobj)
}




############ Adapted functions from Cardinal for binning matrix objects ############


#' This function bins spectral data generated by cardinal into a Matrix object (rather than a matter object)
#'
#' @param x The spectral data to be binned, typically a matrix or data frame.
#' @param ref A vector of reference mass-to-charge (m/z) values for binning. If left unspecified, mass range will be used to generate reference peaks.
#' @param spectra Character string specifying the column name for spectral intensity values (default = "intensity").
#' @param index Character string specifying the column name for the m/z values (default = "mz").
#' @param method A character vector specifying the binning method. Options include `"sum"`, `"mean"`, `"max"`, `"min"`, `"linear"`, `"cubic"`, `"gaussian"`, and `"lanczos"`. If not specified default method used is "sum".
#' @param resolution Numeric value specifying the resolution for binning. If assigned as `NA`, the resolution will be automatically caculated based on the data (default = NA).
#' @param tolerance Numeric value specifying the tolerance for m/z matching (default = `NA`).
#' @param units Character string specifying the units for resolution. Options are `"ppm"` or `"mz"` (default = "ppm").
#' @param mass.range Numeric vector of length two specifying the lower and upper bounds of the mass range for binning. If set to NULL will use the entire range to bin (default = NULL).
#'
#' @return A binned spectral data object with the same structure as the input, but with modified spectral data according to the specified binning parameters.
#' @export
#'
#' @import Cardinal
#' @import matter
#'
#' @examples
#' # used as a helper function for LoadSM
bin_cardinal <-	function(x, ref,
                         spectra = "intensity", index = "mz",
                         method = c("sum", "mean", "max", "min",
                                    "linear", "cubic", "gaussian", "lanczos"),
                         resolution = NA, tolerance = NA, units = c("ppm", "mz"),
                         mass.range = NULL)
{
  if ( !is.null(mass.range) ) {
    if ( is.na(resolution) )
      stop("setting 'mass.range' requires setting 'resolution'")
    ref <- mass.range
  }
  if ( !missing(ref) ) {
    if ( is(ref, "MSImagingExperiment") || is(ref, "MassDataFrame") )
      ref <- mz(ref)
  }
  if ( missing(units) && !missing(resolution) )
    units <- get_units_from_names(resolution, units)

  if (!units %in% c("ppm", "mz")){
    stop("incorrect unit value! units must = either 'ppm' or 'mz'. Please change units accordingly ...")
  }

  units <- ifelse(units == "ppm", "relative", "absolute")

  if ( !is.na(resolution) )
    resolution <- switch(units,
                         relative=1e-6 * resolution,
                         absolute=resolution)
  if ( !is.na(tolerance) )
    tolerance <- ifelse(units == "relative",
                        1e-6 * tolerance,
                        tolerance)

  ans <- bin_SpectralImagingExperiment(x = x, ref=ref, spectra=spectra, index=index,
                                       method=method, resolution=resolution, tolerance=tolerance, units=units)


  return(ans)
}



#' This function bins spectral data generated by cardinal into a Matrix object (rather than a matter object)
#'

#' @param mass.range Numeric vector of length two specifying the lower and upper bounds of the mass range for binning. If set to NULL will use the entire range to bin (default = NULL).
#'



#' Function used to generate a binned intensity matrix based on reference peaks
#'
#' @param x The spectral data to be binned, typically a matrix or data frame.
#' @param ref A vector of reference mass-to-charge (m/z) values for binning. If left unspecified, mass range will be used to generate reference peaks.
#' @param spectra Character string specifying the column name for spectral intensity values (default = "intensity").
#' @param index Character string specifying the column name for the m/z values (default = "mz").
#' @param method A character vector specifying the binning method. Options include `"sum"`, `"mean"`, `"max"`, `"min"`, `"linear"`, `"cubic"`, `"gaussian"`, and `"lanczos"`. If not specified default method used is "sum".
#' @param resolution Numeric value specifying the resolution for binning. If assigned as `NA`, the resolution will be automatically caculated based on the data (default = NA).
#' @param tolerance Numeric value specifying the tolerance for m/z matching (default = `NA`).
#' @param units Character string specifying the units for resolution. Options are `"ppm"` or `"mz"` (default = "ppm").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed.
#'
#' @return Matrix object which has binned intensity values.
#' @export
#'
#' @import Cardinal
#' @import matter
#'
#' @examples
#' #Helper function for binning data in Matrix format
bin_SpectralImagingExperiment <- function(x, ref,
                                          spectra = "intensity", index = NULL,
                                          method = c("sum", "mean", "max", "min",
                                                     "linear", "cubic", "gaussian", "lanczos"),
                                          resolution = NA, tolerance = NA, units = c("relative", "absolute"),
                                          verbose = getCardinalVerbose())
{
  method <- match.arg(method)
  if ( length(index) > 1L )
    stop("more than 1 'index' array not allowed")
  if ( missing(units) && !missing(resolution) )
    units <- get_units_from_names(resolution, units)
  units <- match.arg(units)
  snm <- spectra
  inm <- index
  spectra <- spectra(x, snm)

  if ( is.null(inm) ) {
    inm <- "index"
    domain <- seq_len(nrow(x))
  } else {
    domain <- data.frame(Cardinal::featureData(x))[[inm]]
    if ( is.null(domain) )
      stop("index ", inm, " not found")
  }
  if ( is.sparse(spectra) ) {

    index <- atomindex(spectra)
    #spectra <- atomdata(spectra)
  } else {
    index <- rep.int(list(domain), ncol(x))
  }

  res.ref <- switch(units, relative="x", absolute="abs")
  if ( is.na(resolution) ) {
    if ( missing(ref) || is.null(ref) ) {
      res <- estres(domain, ref=res.ref)
    } else {
      res <- estres(ref, ref=res.ref)
    }
  } else {
    res <- setNames(resolution, units)
    if ( !missing(ref) && !is.null(ref) )
      ref <- switch(units,
                    relative=seq_rel(min(ref), max(ref), by=res),
                    absolute=seq.default(min(ref), max(ref), by=res))
  }

  if ( missing(ref) || is.null(ref) ) {
    from <- floor(min(domain))
    to <- ceiling(max(domain))
    ref <- switch(units,
                  relative=seq_rel(from, to, by=res),
                  absolute=seq.default(from, to, by=res))
  }
  if ( is.na(tolerance) ) {
    if ( method %in% c("sum", "mean", "max", "min") ) {
      tol <- 0.5 * res
    } else {
      tol <- 2 * res
    }
    from <- round(min(ref), digits=4L)
    to <- round(max(ref), digits=4L)
    #.Log("binning ", snm, " from ", inm, " ", from, " to ", to,
    #	" with ", units, " resolution ", round(res, digits=6L),
    #	message=verbose)
  } else {
    label <- if (length(ref) != 1L) "references" else "reference"
    #.Log("binning ", snm, " to ", length(ref), " ", inm, " ", label,
    #	" with ", units, " tolerance ", round(tolerance, digits=6L),
    #	message=verbose)
    tol <- setNames(tolerance, units)
  }

  all_identical <- all(sapply(index, identical, index[[1]]))

  if (all_identical) {
    message("All m/z values are identical between pixels ... ")
    index <-  index[[1]]
  } else {
    stop("Some m/z values are different between pixels ... SpaMTP will use the first list of indexes to bin data! If this is unwanted please set binning_function = 'Cardinal'")
  }

  verbose_message(message_text = "Converting matter matrix to Matrix ... ", verbose = verbose)

  mat <- as.matrix(Cardinal::spectra(x))

  if ( length(index) <= length(ref)){
    rownames(mat) <- index
    return(mat)
  } else {
    return(spectral_binning(matrix=mat, ref = ref, index = index, method = method, tolerance = tol))
  }
}




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




#' Bins m/z peaks of a SpaMTP object.
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

   mtx <- spectral_binning(matrix=as.matrix(data[[assay]][slot]), ref = ref, index = index, method = method, tolerance = tol)
   colnames(mtx) <- rownames(data@meta.data)

   if (return.only.mtx){
     return(mtx)
   } else {
     obj <- SeuratObject::CreateAssayObject(counts = mtx)
     data[["binned"]] <- obj
     return(data)
   }
}




