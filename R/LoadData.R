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



