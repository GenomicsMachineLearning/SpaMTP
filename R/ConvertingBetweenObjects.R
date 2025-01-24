#### SpaMTP Cardinal to Seurat Functions ###############################################################################################################################################################################


#' Converts a Cardinal Object into a Seurat Object
#'
#' @param data A Cardinal Object that is being converted into a Seurat Object.
#' @param multi.run Boolean indicating if there are multiple runs within the imported data. If `TRUE`, an index will be added to the pixel names per run, and an individual FOV will be generated per run in the Seurat Object (default = FALSE).
#' @param seurat.coord A Data.Frame containing two columns titled 'X_new' and 'Y_new' specifying the pixel coordinates of each data point. This is only required if mapping Spatial Metabolic data with a H&E image was done externally, and the SM coordinates need to change to align correctly to the ST data. Else, set to `NULL` (default = NULL).
#' @param assay Character string containing the name of the assay (default = "Spatial").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A Seurat Object containing the mz count data of the supplied Cardinal Object
#' @export
#'
#' @examples
#' # CardinalToSeurat(CardinalObj, run_name = "run_1", seurat.coord = NULL)
CardinalToSeurat <- function(data, multi.run = FALSE, seurat.coord = NULL, assay = "Spatial", verbose = TRUE){

  verbose_message(message_text = "Convering Cardinal object to Seurat object .... ", verbose = verbose)

  #run_data <- Cardinal::subsetPixels(data, Cardinal::run(data) == paste0(run_name)) #subset broken under Cardinal >3.6

  sparse_matrix <- Cardinal::spectra(data)
  pixel_data <- Cardinal::pixelData(data)

  if (!(is.null(seurat.coord))){
    verbose_message(message_text = "Convering Cardinal Coordinates to Seurat Visium Coordinates specified in the seurat.coord file .... ", verbose = verbose)
    pixel_data[["x_coord",]] <- seurat.coord$X_new # changes coordinates to matched Visium Object
    pixel_data[["y_coord",]] <- seurat.coord$Y_new
    Cardinal::pixelData(data) <- pixel_data
  } else{
    if ("x" %in% colnames(data.frame(pixel_data)) & "y" %in% colnames(data.frame(pixel_data))){
      pixel_data_df <- data.frame(pixel_data)
      pixel_data[["x_coord",]] <- pixel_data_df$x
      pixel_data[["y_coord",]] <- pixel_data_df$y
      Cardinal::pixelData(data) <- pixel_data
    } else {
      warning("There is no column called 'x' and 'y' in pixelData(CardinalObject)")
      stop("x and y pixel columns do not exist")
    }
  }

  verbose_message(message_text = "Generating Seurat Barcode Labels from Pixel Coordinates .... ", verbose = verbose)


  x_coord <- Cardinal::pixelData(data)[["x_coord",]]
  y_coord <- Cardinal::pixelData(data)[["y_coord",]]
  spot_name <- paste0(x_coord,"_",y_coord)

  if(multi.run){
    spot_name <- paste(spot_name, as.numeric(as.factor(Cardinal::run(data))), sep = "-")
  }

  colnames(sparse_matrix)<- spot_name
  rownames(sparse_matrix)<- paste("mz-", data.frame(Cardinal::featureData(data))[["mz"]], sep = "")

  verbose_message(message_text = "Constructing Seurat Object ....", verbose = verbose)

  mat <- matter::as.matrix(sparse_matrix)


  seuratobj <- Seurat::CreateSeuratObject(mat, assay = "Spatial")

  verbose_message(message_text = "Adding Pixel Metadata ....", verbose = verbose)

  seuratobj <- Seurat::AddMetaData(seuratobj,col.name = "sample", metadata = Cardinal::run(data))

  for (name in names(Cardinal::pixelData(data))){
    seuratobj <- Seurat::AddMetaData(seuratobj,col.name = name, metadata = Cardinal::pixelData(data)[[name,]])
  }

  verbose_message(message_text = "Creating Centroids for Spatial Seurat Object ....", verbose = verbose)

  ## Add spatial data
  if(multi.run){

    for (i in 1:length(unique(Cardinal::run(data)))){
      idx <- which(Cardinal::run(data) == unique(Cardinal::run(data))[i])
      cents <- SeuratObject::CreateCentroids(data.frame(x = c(Cardinal::pixelData(data)[["x_coord",]][idx]),
                                                        y = c(Cardinal::pixelData(data)[["y_coord",]][idx]),
                                                        cell = c(spot_name)[idx]))
      segmentations.data <- list(
        "centroids" = cents
      )

      coords <- SeuratObject::CreateFOV(
        coords = segmentations.data,
        type = c("centroids"),
        molecules = NULL,
        assay = "Spatial"
      )

      seuratobj[[paste0("fov.",i)]] <- coords

    }

  }else{
    cents <- SeuratObject::CreateCentroids(data.frame(x = c(Cardinal::pixelData(data)[["x_coord",]]), y = c(Cardinal::pixelData(data)[["y_coord",]]), cell = c(spot_name)))


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
  }


  metadata <- data.frame("raw_mz" = sapply(strsplit(rownames(seuratobj), "-"), function(x) as.numeric(x[[2]])))
  rownames(metadata) <- rownames(seuratobj)


  seuratobj[["Spatial"]] <- Seurat::AddMetaData(object = seuratobj[["Spatial"]],
                                                metadata = metadata,
                                                col.name = 'raw_mz')


  seuratobj[["Spatial"]]@meta.data$mz_names <- rownames(seuratobj)

  return(seuratobj)
}

#' Converts a Cardinal Object into a Seurat Object
#'
#' @param data A Cardinal Object that is being converted into a Seurat Object.
#' @param mtx Matrix object containing
#' @param multi.run Boolean indicating if there are multiple runs within the imported data. If `TRUE`, an index will be added to the pixel names per run, and an individual FOV will be generated per run in the Seurat Object (default = FALSE).
#' @param seurat.coord A Data.Frame containing two columns titled 'X_new' and 'Y_new' specifying the pixel coordinates of each data point. This is only required if mapping Spatial Metabolic data with a H&E image was done externally, and the SM coordinates need to change to align correctly to the ST data. Else, set to `NULL` (default = NULL).
#' @param assay Character string containing the name of the assay (default = "Spatial").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A Seurat Object containing the mz count data of the supplied Cardinal Object
#' @export
#'
#' @examples
#' # CardinalToSeurat(CardinalObj, run_name = "run_1")
BinnedCardinalToSeurat <- function(data, mtx, multi.run = FALSE,  assay = "Spatial", verbose = TRUE ){

  verbose_message(message_text = "Convering Cardinal object to Seurat object .... ", verbose = verbose)

  #run_data <- Cardinal::subsetPixels(data, Cardinal::run(data) == paste0(run_name)) #subset broken under Cardinal >3.6

  pixel_data <- Cardinal::pixelData(data)

  if ("x" %in% colnames(data.frame(pixel_data)) & "y" %in% colnames(data.frame(pixel_data))){
    pixel_data_df <- data.frame(pixel_data)
    pixel_data[["x_coord",]] <- pixel_data_df$x
    pixel_data[["y_coord",]] <- pixel_data_df$y
    Cardinal::pixelData(data) <- pixel_data
  } else {
    warning("There is no column called 'x' and 'y' in pixelData(CardinalObject)")
    stop("x and y pixel columns do not exist")
  }

  verbose_message(message_text = "Generating Seurat Barcode Labels from Pixel Coordinates .... ", verbose = verbose)

  x_coord <- Cardinal::pixelData(data)[["x_coord",]]
  y_coord <- Cardinal::pixelData(data)[["y_coord",]]
  spot_name <- paste0(x_coord,"_",y_coord)

  if(multi.run){
    spot_name <- paste(spot_name, as.numeric(as.factor(Cardinal::run(data))), sep = "-")
  }


  colnames(mtx)<- spot_name
  rownames(mtx)<- paste0("mz-", rownames(mtx))

  verbose_message(message_text = "Constructing Seurat Object ....", verbose = verbose)


  seuratobj <- Seurat::CreateSeuratObject(mtx, assay = "Spatial")

  verbose_message(message_text = "Adding Pixel Metadata ....", verbose = verbose)

  seuratobj <- Seurat::AddMetaData(seuratobj,col.name = "sample", metadata = Cardinal::run(data))

  for (name in names(Cardinal::pixelData(data))){
    seuratobj <- Seurat::AddMetaData(seuratobj,col.name = name, metadata = Cardinal::pixelData(data)[[name,]])
  }

  verbose_message(message_text = "Creating Centroids for Spatial Seurat Object ....", verbose = verbose)

  ## Add spatial data
  if(multi.run){

    for (i in 1:length(unique(Cardinal::run(data)))){
      idx <- which(Cardinal::run(data) == unique(Cardinal::run(data))[i])
      cents <- SeuratObject::CreateCentroids(data.frame(x = c(Cardinal::pixelData(data)[["x_coord",]][idx]),
                                                        y = c(Cardinal::pixelData(data)[["y_coord",]][idx]),
                                                        cell = c(spot_name)[idx]))
      segmentations.data <- list(
        "centroids" = cents
      )

      coords <- SeuratObject::CreateFOV(
        coords = segmentations.data,
        type = c("centroids"),
        molecules = NULL,
        assay = "Spatial"
      )

      seuratobj[[paste0("fov.",i)]] <- coords

    }

  }else{
    cents <- SeuratObject::CreateCentroids(data.frame(x = c(Cardinal::pixelData(data)[["x_coord",]]), y = c(Cardinal::pixelData(data)[["y_coord",]]), cell = c(spot_name)))


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
  }


  metadata <- data.frame("raw_mz" = sapply(strsplit(rownames(seuratobj), "-"), function(x) as.numeric(x[[2]])))
  rownames(metadata) <- rownames(seuratobj)


  seuratobj[["Spatial"]] <- Seurat::AddMetaData(object = seuratobj[["Spatial"]],
                                                metadata = metadata,
                                                col.name = 'raw_mz')


  seuratobj[["Spatial"]]@meta.data$mz_names <- rownames(seuratobj)

  return(seuratobj)
}



#' Converts a Seurat object to a Cardinal object, including annotations and metadata
#'
#' @param data Seurat object being converted.
#' @param assay Character string defining the Seurat Object assay name to pull intensity count data from (default = "Spatial").
#' @param slot Character string defining which slot from the Seurat Object assay to gather intensity data from (default = "counts").
#' @param run_col Character string describing the Seurat meta.data column where the run identities are stored (default = NULL).
#' @param feature.metadata Boolean value of whether the Seurat Object contains annotations stored in the feature metadata slot of the specified assay (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A Cardinal object containing intensity values and feature metadata (annotations)
#' @export
#'
#' @examples
#' # cardinal.obj <- ConvertSeuratToCardinal(SeuratObject, feature.metadata = TRUE)
ConvertSeuratToCardinal <- function(data, assay = "Spatial", slot = "counts", run_col = NULL, feature.metadata = FALSE, verbose = TRUE){

  verbose_message(message_text = "Gathering Intensity, m/z values and metadata from Seurat Object ...", verbose = verbose)

  mzs <- unlist(lapply(rownames(data[[assay]]@features), function(x) as.numeric(stringr::str_split(x, "mz-")[[1]][2])))
  coord <- SeuratObject::GetTissueCoordinates(data)

  if (!(is.null(run_col))){
    run <- data@meta.data[[run_col]]
  } else {
    run <- factor("run0")
  }

  pdata <- Cardinal::PositionDataFrame(run = run,
                             coord=coord[c("x","y")],
                             Seurat_metadata = data@meta.data[c(colnames(data@meta.data))]) #adds all other columns

  if (feature.metadata){
    fdata <- Cardinal::MassDataFrame(mz=mzs, feature_metadata = data[[assay]]@meta.data[c(colnames(data[[assay]]@meta.data))])
  } else {
    fdata <- Cardinal::MassDataFrame(mz=mzs)
  }

  mat <- data[[assay]][slot]
  rownames(mat) <- NULL
  colnames(mat) <- NULL

  verbose_message(message_text = "Converting intensity matrix and Generating Cardinal Object ...", verbose = verbose)

  if (check_cardinal_version()){
    cardinal.obj <- Cardinal::MSImagingExperiment(spectraData= matter::sparse_mat(Matrix::as.matrix(mat)),
                                                  featureData=fdata,
                                                  pixelData=pdata)
  } else {
    cardinal.obj <- Cardinal::MSImagingExperiment(imageData= matter::sparse_mat(Matrix::as.matrix(mat)),
                                                  featureData=fdata,
                                                  pixelData=pdata)

  }

  return(cardinal.obj)
}
########################################################################################################################################################################################################################


