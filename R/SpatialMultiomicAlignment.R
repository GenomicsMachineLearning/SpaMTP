#### SpaMTP MALDI to Visium Spot Merging Functions  ####################################################################################################################################################################



#' This function computes the coordinates of a square's four corners based on a given center point and width.
#'
#' @param center_x Numeric value defining the x-coordinate of the square's center.
#' @param center_y Numeric value defining the y-coordinate of the square's center.
#' @param width Numeric value indicating the width of the square.
#' @param name Character string specifying the label assigned to the square for identification.
#'
#' @return A data frame with columns `Selection`, `X`, and `Y` representing the square's name and corner coordinates.
#' @export
#'
#' @examples
#' get_square_coordinates(center_x = 5, center_y = 5, width = 4, name = "MySquare")
get_square_coordinates <- function(center_x, center_y, width, name) {
  # Calculate half width
  half_width <- width / 2

  # Calculate coordinates of the four corners
  x_coords <- c(center_x - half_width, center_x + half_width, center_x + half_width, center_x - half_width, center_x - half_width)
  y_coords <- c(center_y - half_width, center_y - half_width, center_y + half_width, center_y + half_width, center_y - half_width)

  # Return the coordinates as a matrix
  square_df <- data.frame(
    Selection = rep(name, n = 4),
    X = x_coords,
    Y = y_coords

  )
  return(square_df)
}


#' Maps Spatial Metabolomic (MALDI) data to corresponding Spatial Transcriptomics data and coordinates.
#'
#' @param SM.data A SpaMTP Seurat object representing the Spatial Metabolomics data.
#' @param ST.data A Seurat object representing the Spatial Transcriptomics data.
#' @param ST.hires Boolean string defining if the ST data is at a higher resolution compared to the SM pixel data. For example, generally Visium data will be lower res whereas Xenium/single-cell resolution spatial data will be a higher resolution (default = FALSE).
#' @param SM.assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values (default = "Spatial").
#' @param ST.assay Character string specifying the current assay to use to extract transcriptional data from (default = "Spatial").
#' @param SM.fov Character string of the image fov associated with the spatial metabolomic data (default = "fov").
#' @param ST.image Character string matching the image name associated with the ST data such as 'fov' or 'slice1' object (default = "slice1").
#' @param ST.scale.factor Character string defining the image resolution associated with the Visium image pixel data. If `NULL` the full-res coordinates will be used and no scaling will be performed. Note: This parameter is only required for aligning lowres ST data (default = "hires").
#' @param SM.pixel.width Numeric value defining the width of each SM pixel. If set to `NULL`, the median pixel width will be calculated based on the distance between each pixel (default = NULL).
#' @param overlap.threshold Numeric value defining the overlap proportion threshold for a SM pixel to be associated with a ST spot. For example, if res_increase = 0.2 then pixels that have at least 20% area overlap with the respective visium spot will be assigned a match. Note: This parameter is only required for aligning lowres ST data (default = 0.2).
#' @param annotations Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat object contains annotations assigned to m/z values (default = TRUE).
#' @param add.metadata Boolean defining whether to add the current metadata stored in the SM object to the new mapped multi-omic SpaMTP object (default = TRUE)
#' @param merge.unique.metadata Boolean indicating whether to summaries duplicated metadata terms to only store unique values in the metadata. Note: This parameter is only required for aligning lowres ST data, and `add.metadata` must be set to `TRUE` for this functionality to be implemented (default = TRUE).
#' @param map.data Boolean indicating whether to map normalised/additional data stored in the `data` slot of the SpaMTP assay. Note: this process is computationally expensive with large datasets (default = FALSE).
#' @param new_SPT.assay Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated transcriptomics data (default = "SPT").
#' @param new_SPM.assay Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated metabolomic data (default = "SPM").
#' @param verbose Boolean value indicating whether to print informative progression update messages and progress bars (default = TRUE).
#'
#' @return A SpaMTP Seurat object with Spatial Metabolomic data mapped to equivalent Spatial Transcripomics coordinates (Visium spots/Xenium cells).
#' @export
#'
#' @examples
#'
#' ## Mapping MALDI data to equivalent Visium spots
#' # MapSpatialOmics(VisiumObj, SeuratObj, ST.scale.factor = "hires", SM.assay = "Spatial", ST.assay = "Spatial")
#'
#' #' ## Mapping MALDI data to equivalent Xenium cells
#' # MapSpatialOmics(VisiumObj, SeuratObj, SM.assay = "Spatial", ST.assay = "Xenium")
MapSpatialOmics <- function(SM.data, ST.data, ST.hires = FALSE,
                            SM.assay = "Spatial", ST.assay = "Spatial",
                            SM.fov = "fov", ST.image = "slice1",
                            ST.scale.factor = "hires",
                            SM.pixel.width = NULL,
                            overlap.threshold = 0.2,
                            annotations = TRUE,
                            add.metadata = TRUE,
                            merge.unique.metadata = TRUE,
                            map.data = FALSE,
                            new_SPT.assay = "SPT",
                            new_SPM.assay = "SPM",
                            verbose = FALSE) {


  if (ST.hires){
    verbose_message(message_text = "Running `MapSpatialOmics` in hires mode! This is normally used for single-cell spatial data (Xenium)... \n", verbose = verbose)
    verbose_message(message_text = "Inputs of `ST.scale.factor`, `overlap.threshold` and `merge.unique.metadata` are not used for hires mode, and will be ignored ... \n", verbose = verbose)

    mapped.data <- hiresMapping(SM.data = SM.data,
                                ST.data = ST.data,
                                SM.assay = SM.assay,
                                ST.assay = ST.assay,
                                SM.fov = SM.fov,
                                ST.image = ST.image,
                                SM.pixel.width = SM.pixel.width,
                                annotations = annotations,
                                add.metadata = add.metadata,
                                map.data = map.data,
                                new_SPT.assay = new_SPT.assay,
                                new_SPM.assay = new_SPM.assay,
                                verbose = verbose)

  } else {
    verbose_message(message_text = "Running `MapSpatialOmics` in lowres mode! This is normally used for bin/spot-based spatial data (Visium)... \n", verbose = verbose)
    verbose_message(message_text = "Inputs of `ST.scale.factor`, `overlap.threshold` and `merge.unique.metadata` are required for lowres mode, please ensure they are included ... \n", verbose = verbose)

    mapped.data <- lowresMapping(SM.data = SM.data,
                                  ST.data = ST.data,
                                  SM.assay = SM.assay,
                                  ST.assay = ST.assay,
                                  SM.fov = SM.fov,
                                  ST.image = ST.image,
                                  ST.scale.factor = ST.scale.factor,
                                  SM.pixel.width = SM.pixel.width,
                                  overlap.threshold = overlap.threshold,
                                  annotations = annotations,
                                  add.metadata = add.metadata,
                                  merge.unique.metadata = merge.unique.metadata,
                                  map.data = map.data,
                                  new_SPT.assay = new_SPT.assay,
                                  new_SPM.assay = new_SPM.assay,
                                  verbose = verbose)

  }

  return(mapped.data)

}



#' Function used by MapSpatialOmics to align SM data to ST spots with lower resolution (i.e. Visium Spots)
#'
#' @param SM.data A SpaMTP Seurat object representing the Spatial Metabolomics data.
#' @param ST.data A Seurat object representing the Visium Spatial Transcriptomics data.
#' @param SM.assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values (default = "Spatial").
#' @param ST.assay Character string specifying the current assay to use to extract transcriptional data from (default = "Spatial").
#' @param SM.fov Character string of the image fov associated with the spatial metabolomic data (default = "fov").
#' @param ST.image Character string matching the image name associated with the ST data stored in the Visium object (default = "slice1").
#' @param ST.scale.factor Character string defining the image resolution associated with the Visium image pixel data. If `NULL` the full-res coordinates will be used and no scaling will be performed (default = "hires").
#' @param SM.pixel.width Numeric value defining the width of each SM pixel. If set to `NULL`, the median pixel width will be calculated based on the distance between each pixel (default = NULL).
#' @param overlap.threshold Numeric value defining the overlap proportion threshold for a SM pixel to be associated with a ST spot. For example, if res_increase = 0.2 then pixels that have at least 20% area overlap with the respective visium spot will be assigned a match (default = 0.2).
#' @param annotations Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat object contains annotations assigned to m/z values (default = TRUE).
#' @param add.metadata Boolean defining whether to add the current metadata stored in the SM object to the new mapped multi-omic SpaMTP object (default = TRUE)
#' @param merge.unique.metadata Boolean indicating whether to summaries duplicated metadata terms to only store unique values in the metadata. Note: `add.metadata` must be set to `TRUE` for this functionality to be implmented (default = TRUE).
#' @param map.data Boolean indicating whether to map normalised/additional data stored in the `data` slot of the SpaMTP assay. Note: this process is computationally expensive with large datasets (default = FALSE).
#' @param new_SPT.assay Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated transcriptomics data (default = "SPT").
#' @param new_SPM.assay Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated metabolomic data (default = "SPM").
#' @param verbose Boolean value indicating whether to print informative progression update messages and progress bars (default = TRUE).
#'
#' @import Seurat
#' @import SeuratObject
#'
#' @return A SpaMTP Seurat object with the Spatial Metabolomic data mapped to equivalent Spatial Transcripomics (Visium) spots.
#'
#' @examples
#' # Helper function for MapSpatialOmics
lowresMapping <- function(SM.data, ST.data,
                      SM.assay = "Spatial", ST.assay = "Spatial",
                      SM.fov = "fov", ST.image = "slice1",
                      ST.scale.factor = "hires",
                      SM.pixel.width = NULL,
                      overlap.threshold = 0.2,
                      annotations = TRUE,
                      add.metadata = TRUE,
                      merge.unique.metadata = TRUE,
                      map.data = FALSE,
                      new_SPT.assay = "SPT",
                      new_SPM.assay = "SPM",
                      verbose = FALSE){

  sample_data <- ST.data

  ## make polygons from coordinates
  df <- SeuratObject::GetTissueCoordinates(SM.data, image = SM.fov)

  if (is.null(SM.pixel.width)) {

    verbose_message(message_text = "No pixel width was provided, calculating meadian SM pixel width ... \n", verbose = verbose)

    diff_x <- diff(df$x)
    diff_y <- diff(df$y)

    # Combine the x and y differences into a data frame
    differences <- data.frame(diff_x, diff_y)

    # Filter out cases where the difference is zero
    non_zero_differences <- differences[differences$diff_x != 0 | differences$diff_y != 0,]

    # Choose the non-zero difference (maximum absolute value)
    new_widths <- apply(non_zero_differences, 1, function(row) max(abs(row)))
    SM.pixel.width <- stats::median(new_widths)
  }

  verbose_message(message_text = paste0("SM pixel width used for mapping = ", SM.pixel.width, "... \n"), verbose = verbose)

  polygon_list <- lapply(1:nrow(df), function(idx) {
    row <- df[idx,]
    get_square_coordinates(center_x = row$x, center_y = row$y, width = SM.pixel.width, name = row$cell)
  })

  # Combine the list of data frames into a single data frame
  verbose_message(message_text = "Generating polygon from SM pixel coordinates ... \n", verbose = verbose)

  polygons_df <- do.call(rbind, polygon_list)
  polygons <- SeuratObject::CreateSegmentation(polygons_df)


  SM.data[[SM.fov]][["annotation"]] <- polygons
  pol <- SeuratObject::GetTissueCoordinates(SM.data[[SM.fov]][["annotation"]])

  verbose_message(message_text = "Generating polygon from ST spot coordinates ... \n", verbose = verbose)

  ###Make list of polygons from annotations
  active_polygons <- SpatialPolygons(lapply(unique(pol$cell), function(cell_name) {
    cell_data <- pol[pol$cell == cell_name, ]
    Polygons(list(Polygon(cbind(cell_data$x, cell_data$y))), cell_name)}),
    1:length(unique(pol$cell)))

  empty_poly_df = data.frame(cell = unique(pol$cell))
  rownames(empty_poly_df) <- empty_poly_df$cell


  df_cells_spdf <- SpatialPolygonsDataFrame(active_polygons, empty_poly_df)

  polygon_df <- sf::st_as_sf(df_cells_spdf)

  verbose_message(message_text = "Assigning SM polygons to overlapping ST polygons ... \n", verbose = verbose)

  ## Find Cells within each polygon
  st_coordinates <- SeuratObject::GetTissueCoordinates(sample_data[[ST.image]][["centroids"]])
  st_coordinates$radius = sample_data[[ST.image]][["centroids"]]@radius * 0.5 #### 10X scalefactors_json.json file states @radius is actually = spot diameter

  if(!is.null(ST.scale.factor)){
    if(!ST.scale.factor %in% c("hires", "lowres")){
      stop("Invalid asignment of `ST.scale.factor`! values must be either 'hires' or 'lowres'. If fullres is required set `ST.scale.factor` = `NULL`.")
    } else {
      st_coordinates[c("x","y", "radius")] <- st_coordinates[c("x","y", "radius")] * sample_data[[ST.image]]@scale.factors[[ST.scale.factor]]
    }
  }


  points_sf <- sf::st_as_sf(st_coordinates, coords = c("x", "y"))
  points_sf <- sf::st_buffer(points_sf, dist = st_coordinates$radius)

  # Find intersecting polygons
  buffer_areas <- sf::st_area(points_sf)

  # Find intersecting polygons
  intersections <- sf::st_intersects(points_sf, polygon_df, sparse = FALSE)

  # Create a data frame with results
  result <- do.call(rbind, lapply(seq_len(nrow(points_sf)), function(i) {
    # Get the geometry of the buffered point
    point_geom <- points_sf$geometry[i]
    point_area <- buffer_areas[i]

    # Find polygons intersecting with the point
    intersecting_polygons <- polygon_df[intersections[i, ], ]

    # Calculate intersection areas
    intersection_areas <- sf::st_area(sf::st_intersection(point_geom, intersecting_polygons$geometry))

    # Calculate percentage of coverage
    coverage_percentage <- as.numeric(intersection_areas / point_area)

    # Filter polygons by coverage percentage
    valid_polygons <- intersecting_polygons$cell[coverage_percentage >= overlap.threshold]

    data.frame(
      point = st_coordinates$cell[i],
      polygons = paste(valid_polygons, collapse = ", ")
    )
  }))


  SpaMTP.obj <- AddMetaData(ST.data, result$polygons,"SPM_pixels")

  # remove cells with no SM data
  cells_with_na <- rownames(SpaMTP.obj@meta.data)[is.na(SpaMTP.obj@meta.data$SPM_pixels)]

  verbose_message(message_text = "Generating new Spatial Multi-Omic SpaMTP Seurat Object ... \n", verbose = verbose)


  # Subset the Seurat object based on cells with NA values
  SpaMTP.obj <- subset(SpaMTP.obj, cells = cells_with_na, invert = TRUE)

  mean_counts <- function(spm_pixels, count_matrix) {
    # Split SPM_pixels into individual pixel strings
    pixels <- strsplit(spm_pixels, ",\\s*")[[1]]

    if (length(pixels) > 1){
      row_means <- rowMeans(count_matrix[,pixels]) # Compute mean for each feature
    } else if (length(pixels) == 1){
      row_means <- count_matrix[,pixels]
    } else {
      row_means <- rep(0, length(rownames(count_matrix)))
    }
    return(row_means)
  }

  verbose_message(message_text = "Averaging SM intensity values per ST spot ... \n", verbose = verbose)

  if (map.data) {

    if ("data" %in% Layers(SM.data, assay = SM.assay)){

      mean_count_list <- lapply(SpaMTP.obj$SPM_pixels, mean_counts, count_matrix = SM.data[[SM.assay]]["counts"])
      mean_counts_df <- do.call(rbind, mean_count_list)

      verbose_message(message_text = "Averaging normalised SM intensity values stored in the `data` slot per ST spot ... \n", verbose = verbose)

      mean_data_list <- lapply(SpaMTP.obj$SPM_pixels, mean_counts, count_matrix = SM.data[[SM.assay]]["data"])
      mean_data_df <- do.call(rbind, mean_count_list)
      SpaMTP.obj[["SPM"]] <- CreateAssay5Object(counts = Matrix::t(mean_counts_df), data = Matrix::t(mean_data_df))

    } else {
      stop("data slot not present in SPM seurat object provided! Please run NormaliseSMData() or set map.data = FALSE")
    }
  } else {
    # Apply the function to each row in the reference data frame
    mean_count_list <- lapply(SpaMTP.obj$SPM_pixels, mean_counts, count_matrix = SM.data[[SM.assay]]["counts"])
    mean_counts_df <- do.call(rbind, mean_count_list)
    colnames(mean_counts_df) <- rownames(SM.data)
    SpaMTP.obj[[new_SPM.assay]] <- CreateAssay5Object(counts = Matrix::t(mean_counts_df))
  }


  #add metadata from SM

  add_metadata <- function(spm_pixels, SM_metadata) {
    # Split SPM_pixels into individual pixel strings
    pixels <- strsplit(spm_pixels, ",\\s*")[[1]]

    if (length(pixels) == 0){

      col_names <- colnames(SM_metadata)
      combined_df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(col_names)))
      colnames(combined_df) <- paste0(new_SPM.assay,"_",col_names)

    } else {
      sub.metadata <- SM_metadata[pixels,]
      # Apply the function to each column
      if (merge.unique.metadata){

        combined_row <- sapply(sub.metadata, function(x){
          paste(unique(x), collapse = ", ")
        })
      } else {
        combined_row <- sapply(sub.metadata, function(x){
          paste(unique(x), collapse = ", ")
        })
      }

      # Convert to data frame
      combined_df <- as.data.frame(t(combined_row), stringsAsFactors = FALSE)
      colnames(combined_df) <- paste0(new_SPM.assay,"_",colnames(sub.metadata))

    }
    return(combined_df)
  }


  if (add.metadata) {

    verbose_message(message_text = "Adding SM metadata to the new SpaMTP Seurat Object ... \n", verbose = verbose)

    metadata_list <- lapply(SpaMTP.obj$SPM_pixels, add_metadata, SM_metadata = SM.data@meta.data)
    metadata_df <- do.call(rbind, metadata_list)

    SpaMTP.obj@meta.data[colnames(metadata_df)] <- metadata_df[colnames(metadata_df)]
  }



  if (annotations) {

    verbose_message(message_text = "Adding metabolite annotation metadata to the new SpaMTP Seurat Object ... \n", verbose = verbose)

    ## adds m/z annotations to new object
    SpaMTP.obj[[new_SPM.assay]]@meta.data <- SM.data[[SM.assay]]@meta.data
  }

  SpaMTP.obj <- RenameAssays(SpaMTP.obj, assay.name = ST.assay, new.assay.name = new_SPT.assay)


  return(SpaMTP.obj)

}



#' Function used by MapSpatialOmics to align SM data to ST spots with higher resolution (i.e. Xenium cells)
#'
#' @param SM.data A SpaMTP Seurat object representing the Spatial Metabolomics data.
#' @param ST.data A Seurat object representing the Xenium Spatial Transcriptomics data.
#' @param SM.assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values (default = "Spatial").
#' @param ST.assay Character string specifying the current assay to use to extract transcriptional data from (default = "Spatial").
#' @param SM.fov Character string of the image fov associated with the spatial metabolomic data (default = "fov").
#' @param ST.image Character string matching the image name associated with the ST data stored in the Xenium object (default = "slice1").
#' @param SM.pixel.width Numeric value defining the width of each SM pixel. If set to `NULL`, the median pixel width will be calculated based on the distance between each pixel (default = NULL).
#' @param annotations Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat object contains annotations assigned to m/z values (default = TRUE).
#' @param add.metadata Boolean defining whether to add the current metadata stored in the SM object to the new mapped multi-omic SpaMTP object (default = TRUE)
#' @param map.data Boolean indicating whether to map normalised/additional data stored in the `data` slot of the SpaMTP assay. Note: this process is computationally expensive with large datasets (default = FALSE).
#' @param new_SPT.assay Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated transcriptomics data (default = "SPT").
#' @param new_SPM.assay Character string defining the assay name of the new overlaid SpaMTP Seurat object containing all updated metabolomic data (default = "SPM").
#' @param verbose Boolean value indicating whether to print informative progression update messages and progress bars (default = TRUE).
#'
#' @import Seurat
#' @import SeuratObject
#'
#' @return A SpaMTP Seurat object with the Spatial Metabolomic data mapped to equivalent Spatial Transcripomics (Xenium) cells.
#'
#' @examples
#' # Helper function for MapSpatialOmics
hiresMapping <- function(SM.data, ST.data,
                          SM.assay = "Spatial", ST.assay = "Xenium",
                          SM.fov = "fov", ST.image = "fov",
                          SM.pixel.width = NULL,
                          annotations = TRUE,
                          add.metadata = TRUE,
                          map.data = FALSE,
                          new_SPT.assay = "SPT",
                          new_SPM.assay = "SPM",
                          verbose = FALSE){

    sample_data <- ST.data

    ## make polygons from coordinates
    df <- GetTissueCoordinates(SM.data, image = SM.fov)

    if (is.null(SM.pixel.width)) {

      verbose_message(message_text = "No pixel width was provided, calculating meadian SM pixel width ... \n", verbose = verbose)
      diff_x <- diff(df$x)
      diff_y <- diff(df$y)

      # Combine the x and y differences into a data frame
      differences <- data.frame(diff_x, diff_y)

      # Filter out cases where the difference is zero
      non_zero_differences <- differences[differences$diff_x != 0 | differences$diff_y != 0,]

      # Choose the non-zero difference (maximum absolute value)
      new_widths <- apply(non_zero_differences, 1, function(row) max(abs(row)))
      SM.pixel.width <- median(new_widths)
    }

    verbose_message(message_text = paste0("SM pixel width used for mapping = ", SM.pixel.width, "... \n"), verbose = verbose)

    polygon_list <- lapply(1:nrow(df), function(idx) {
      row <- df[idx,]
      get_square_coordinates(center_x = row$x, center_y = row$y, width = SM.pixel.width, name = row$cell)
    })

    verbose_message(message_text = "Generating polygon from SM pixel coordinates ... \n", verbose = verbose)

    # Combine the list of data frames into a single data frame
    polygons_df <- do.call(rbind, polygon_list)
    polygons <- CreateSegmentation(polygons_df)

    sample_data[[ST.image]][["annotation"]] <- polygons
    pol <- GetTissueCoordinates(sample_data[[ST.image]][["annotation"]])


    ###Make list of polygons from annotations
    active_polygons <- SpatialPolygons(lapply(unique(pol$cell), function(cell_name) {
      cell_data <- pol[pol$cell == cell_name, ]
      Polygons(list(Polygon(cbind(cell_data$x, cell_data$y))), cell_name)}),
      1:length(unique(pol$cell)))

    empty_poly_df = data.frame(cell = unique(pol$cell))
    rownames(empty_poly_df) <- empty_poly_df$cell


    df_cells_spdf <- SpatialPolygonsDataFrame(active_polygons, empty_poly_df)

    verbose_message(message_text = "Assigning SM polygons to overlapping ST polygons ... \n", verbose = verbose)

    ## Find Cells within each polygon
    xenium_coordinates <- GetTissueCoordinates(sample_data[[ST.image]][["centroids"]])

    df_points_spdf <- SpatialPointsDataFrame(xenium_coordinates[, c("x", "y")], data = xenium_coordinates)
    result <- over(df_points_spdf, df_cells_spdf)

    xenium.data <- AddMetaData(ST.data, result$cell,"SPM_pixels")

    verbose_message(message_text = "Generating new Spatial Multi-Omic SpaMTP Seurat Object ... \n", verbose = verbose)


    # remove cells with no SM data
    cells_with_na <- rownames(xenium.data@meta.data)[is.na(xenium.data@meta.data$SPM_pixels)]

    # Subset the Seurat object based on cells with NA values
    xenium.data <- subset(xenium.data, cells = cells_with_na, invert = TRUE)

    # generate new SM counts matrix
    MALDI_df <- SM.data[[SM.assay]]["counts"][, xenium.data@meta.data$SPM_pixels]

    colnames(MALDI_df) <- rownames(xenium.data@meta.data)

    if (map.data) {

      verbose_message(message_text = "Mapping normalised SM data stored in `data` slot to new SpaMTP Seurat Object ... \n", verbose = verbose)

      MALDI_data_df <-  SM.data[[SM.assay]]["data"][, xenium.data@meta.data$SPM_pixels]
      colnames(MALDI_data_df) <- rownames(xenium.data@meta.data)
      xenium.data[[new_SPM.assay]] <- CreateAssay5Object(counts = MALDI_df, data = MALDI_data_df)
    } else {
      xenium.data[[new_SPM.assay]] <- CreateAssay5Object(counts = MALDI_df)
    }


    #add metadata from SM
    MALDI_metadata <- SM.data@meta.data[xenium.data@meta.data$SPM_pixels,]
    rownames(MALDI_metadata) <- rownames(xenium.data@meta.data)

    if (add.metadata) {
      verbose_message(message_text = "Adding SM metadata to the new SpaMTP Seurat Object ... \n", verbose = verbose)

      meta_data_colnames <- paste0(colnames(MALDI_metadata), "_", new_SPM.assay)
      xenium.data@meta.data[meta_data_colnames] <- MALDI_metadata
    }


    xenium.data <- SeuratObject::RenameAssays(object = xenium.data, assay.name = ST.assay, new.assay.name = new_SPT.assay, verbose = verbose)

    if (annotations) {
      verbose_message(message_text = "Adding metabolite annotation metadata to the new SpaMTP Seurat Object ... \n", verbose = verbose)

      ## adds m/z annotations to new object
      xenium.data[[new_SPM.assay]]@meta.data <- SM.data[[SM.assay]]@meta.data
    }


    return(xenium.data)

}


########################################################################################################################################################################################################################




#### SpaMTP Manual Alignment of ST and SM data ####################################################################################################################################################################
# Code below and some function have been modified from STUtility: https://github.com/jbergenstrahle/STUtility/tree/master


#' Shiny app allowing for manual alignment of SM and ST data coordinates
#'
#' @param sm.data SpaMTP Seurat Object containing SM data
#' @param st.data SpaMTP Seurat Object containing ST data
#' @param msi.pixel.multiplier Numeric value defining a scale.factor to multiple each SM pixel coordinates by (default = 20).
#' @param image.res Character string of the corresponding ST image scale factor to use (default = "lowres").
#' @param continous_cols Vector of colours to use for plotting continuous data. If NULL, the colour map "Reds" will be used (default = NULL).
#' @param catagorical_cols Vector of colours to use for plotting categorical data (default = NULL).
#' @param fov Character string matching the name of the SM FOV to use for plotting (default = "fov").
#' @param image.slice Character string matching the ST image slice name to use for plotting (default = "slice1").
#' @param shiny.host Character string of the shiny host network interface that the Shiny application will listen on when run (default = "0.0.0.0").
#' @param shiny.port Numeric 4 digit number defining the port that the Shiny application will listen to when run (default = 4698).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).
#'
#' @return A SpaMTP Seurat Object containing SM data with transformed coordinated to match the aligned ST data
#' @export
#'
#' @importFrom shiny runApp fluidPage modalDialog fluidRow column sliderInput checkboxInput selectInput actionButton plotOutput reactive
#' renderPlot eventReactive observe stopApp h4 numericInput HTML showModal
#' @importFrom shinyjs useShinyjs reset
#' @importFrom zeallot %<-%
#'
#'
#' @examples
#' # SM_Transformed <- AlignSpatialOmics(SM.data, ST.data)
AlignSpatialOmics <- function (
    sm.data,
    st.data,
    msi.pixel.multiplier = 20,
    image.res = "lowres",
    continous_cols = NULL,
    catagorical_cols = NULL,
    fov = "fov",
    image.slice = "slice1",
    shiny.host = "0.0.0.0",
    shiny.port = 4698,
    verbose = FALSE

) {

  options(shiny.host = shiny.host)
  options(shiny.port = shiny.port)

  #Get tissue coordinates from Seurat Objects

  ## ST cooridnates
  df <- GetTissueCoordinates(st.data)[c("x", "y")] * st.data@images[[image.slice]]@scale.factors[[image.res]]

  ## SM Coordinates
  df2 <- GetTissueCoordinates(sm.data)[c("x", "y")]
  df2$x <- df2$x * msi.pixel.multiplier / (st.data@images[[image.slice]]@scale.factors[["hires"]]/st.data@images[[image.slice]]@scale.factors[[image.res]])
  df2$y <- df2$y * msi.pixel.multiplier / (st.data@images[[image.slice]]@scale.factors[["hires"]]/st.data@images[[image.slice]]@scale.factors[[image.res]])
  df2 <- df2[c("x", "y")]

  # Calculate scatter for plotting
  df$pixel_x <- df$x
  df$pixel_y <- df$y

  sc1 <- df[c("x", "y")]
  rownames(sc1) <- NULL
  coords1 <- df[c("pixel_x", "pixel_y")]

  df2$pixel_x <- df2$x
  df2$pixel_y <- df2$y

  sc2 <- df2[c("x", "y")]
  rownames(sc2) <- NULL
  coords2<- df2[c("pixel_x", "pixel_y")]

  sc <- list("1" = list("scatter" = sc1, "coords" = coords1),
             "2" = list("scatter" = sc2, "coords" = coords2))


  arr <- st.data@images[[image.slice]]@image
  rotated_array <- aperm(arr, c(2, 1, 3))
  rotated_array <- rotated_array[ nrow(rotated_array):1, ,]
  color_matrix <- (as.raster(rotated_array))


  reference.index = 1
  align.index = 2
  scatters <- sc
  fixed.scatter <- scatters[[reference.index]]$scatter
  counter <- NULL
  coords.ls <- NULL
  transformations <-  list(diag(c(1, 1, 1)), diag(c(1, 1, 1)))
  tr.matrices <- lapply(transformations, function(x) diag(c(1, 1, 1)))

  id <- list("1" = dim(color_matrix), "2" = dim(color_matrix))
  image.dims <- id



  ui <- fluidPage(
    useShinyjs(),
    fluidRow(
      column(4,
             shiny::hr(),
             actionButton(inputId = "info", label = "Instructions"),
             shiny::hr(),
             fluidRow(
               column(width = 6, sliderInput(
                 inputId = "angle",
                 label = "Rotation angle",
                 value = 0, min = -120, max = 120, step = 0.1
               ))
             ),
             fluidRow(
               column(width = 6, sliderInput(
                 inputId = "shift_x",
                 label = "Move along x axis",
                 value = 0, min = -round(dim(color_matrix)[2]*(3/4)), max = round(dim(color_matrix)[2]*(3/4)), step = 1
               )),
               column(width = 6, sliderInput(
                 inputId = "shift_y",
                 label = "Move along y axis",
                 value = 0, min = -round(dim(color_matrix)[2]*(3/4)), max = round(dim(color_matrix)[2]*(3/4)), step = 1
               ))
             ),
             h4("stretch along blue axis:"),
             fluidRow(
               column(width = 6, sliderInput(
                 inputId = "stretch_angle1",
                 label = "angle",
                 value = 0, min = -180, max = 180, step = 0.1
               )),
               column(width = 6, sliderInput(
                 inputId = "stretch_factor1",
                 label = "stretch/squeeze",
                 value = 1, min = 0.1, max = 2, step = 0.01
               ))
             ),
             h4("stretch along red axis:"),
             fluidRow(
               column(width = 6, sliderInput(
                 inputId = "stretch_angle2",
                 label = "angle",
                 value = 0, min = -180, max = 180, step = 0.1
               )),
               column(width = 6, sliderInput(
                 inputId = "stretch_factor2",
                 label = "stretch/squeeze",
                 value = 1, min = 0.1, max = 2, step = 0.01
               ))
             ),
             fluidRow(
               column(4, numericInput(
                 inputId = "size_spot",
                 label = "SM spot size",
                 value = 0.5, min = 0, max = 5, step = 0.1
               )),
               column(4, numericInput(
                 inputId = "size_target",
                 label = "ST point size",
                 value = 0.3, min = 0, max = 5, step = 0.05
               )),
               column(4, selectInput(
                 inputId = "spot_shape",
                 label = "spot shape",
                 choices =   c("spot" = "circle",
                               "pixel" = "square")
               )),
               column(4, selectInput(
                 inputId = "sm_plot",
                 label = "SM plot feature",
                 choices =   setNames(colnames(sm.data@meta.data), colnames(sm.data@meta.data))
               )),
               column(4, selectInput(
                 inputId = "st_plot",
                 label = "ST plot feature",
                 choices =   setNames(colnames(st.data@meta.data), colnames(st.data@meta.data))
               ))
             ),
             fluidRow(

               column(4,  checkboxInput(inputId = "show_ST_img",
                                        label = "show image",
                                        value = TRUE)
               ),column(4,  checkboxInput(inputId = "show_ST_spots",
                                        label = "show ST spots",
                                        value = FALSE)
               ),
               column(4,  checkboxInput(inputId = "show_SM_spots",
                                        label = "show SM spots",
                                        value = TRUE)
               ),
               column(4,  checkboxInput(inputId = "flip_x",
                                        label = "mirror along x axis",
                                        value = FALSE)
               ),
               column(4,  checkboxInput(inputId = "flip_y",
                                        label = "mirror along y axis",
                                        value = FALSE)
               )

             ),
             #selectInput(inputId = "sample", choices = (1:length(scatters))[-reference.index],

            #             label = "Select sample", selected = reference.index),
             actionButton("myBtn", "Return aligned data")
      ),

      column(7, plotOutput("scatter")
      )
    )
  )

  server <- function(input, output) {

    rotation_angle <- reactive({
      rot_angle <- input$angle
      return(rot_angle)
    })

    translation_xy <- reactive({
      trxy <- c(input$shift_x, input$shift_y)
      return(trxy)
    })

    mirror_xy <- reactive({
      mirrxy <- c(input$flip_x, input$flip_y)
      return(mirrxy)
    })

    stretch_angle1 <- reactive({
      str_angle1 <- input$stretch_angle1
      return(str_angle1)
    })

    stretch_factor1 <- reactive({
      str_factor1 <- input$stretch_factor1
      return(str_factor1)
    })

    stretch_angle2 <- reactive({
      str_angle2 <- input$stretch_angle2
      return(str_angle2)
    })

    stretch_factor2 <- reactive({
      str_factor2 <- input$stretch_factor2
      return(str_factor2)
    })


    pt_size <- reactive({
      input$size_spot
    })

    st_feature_plot <- reactive({
      input$st_plot
    })

    sm_feature_plot <- reactive({
      input$sm_plot
    })

    pt_shape <- reactive({
      if (input$spot_shape == "square"){
        return(15)
      } else{
        return(19)
      }
    })

    pt_size_target <- reactive({
      input$size_target
    })


    pt_st_points <- reactive({
      input$show_ST_spots
    })

    pt_sm_points <- reactive({
      input$show_SM_spots
    })

    pt_st_img <- reactive({
      input$show_ST_img
    })



    coords_list <- reactive({

      # Obtain point set and spot pixel coordinates
      ls <- scatter.coords()
      scatter.t <- ls[[1]]; coords.t <- ls[[2]]

      # Set transformation parameters
      xt.yt <- translation_xy()
      xy.alpha <- rotation_angle()
      mirrxy <-  mirror_xy()
      str.alpha1 <- stretch_angle1()
      str.factor1 <- stretch_factor1()
      str.alpha2 <- stretch_angle2()
      str.factor2 <- stretch_factor2()

      # Apply reflections
      center <- apply(scatter.t, 2, mean)
      tr.mirror <- mirror(mirror.x = mirrxy[1], mirror.y = mirrxy[2], center.cur = center)

      # Apply rotation
      tr.rotate <- rotate(angle = -xy.alpha, center.cur = center)

      # Apply translation
      tr.translate <- translate(translate.x = xt.yt[1], translate.y = -xt.yt[2])

      # Apply stretch
      tr.stretch1 <- stretch(r = str.factor1, alpha = -str.alpha1, center.cur = center)
      tr.stretch2 <- stretch(r = str.factor2, alpha = -(str.alpha2 + 90), center.cur = center)

      # Combine transformations
      tr <- tr.stretch2%*%tr.stretch1%*%tr.translate%*%tr.rotate%*%tr.mirror


      # Apply transformations
      scatter.t <- t(tr%*%rbind(t(scatter.t), 1))[, 1:2]
      coords.t <- t(tr%*%rbind(t(coords.t), 1))[, 1:2]

      return(list(scatter = scatter.t, coords = coords.t, tr = tr, xylimits = image.dims[[align.index]]))
    })

    output$scatter <- renderPlot({

      coords.ls <<- coords_list()
      c(scatter.t, coords.t, tr, xylimit) %<-% coords.ls

      d <- round((sqrt(xylimit[1]^2 + xylimit[2]^2) - xylimit[2])/2)

      center <- apply(coords.t[, 1:2], 2, mean)

      arrows.1 <- function(x0, y0, length.ar, angle.ar, ...){

        angle.ar <- 2*pi*(-angle.ar/360)
        ab <- cos(angle.ar) * length.ar
        bc <- sign(sin(angle.ar)) * sqrt(length.ar^2 - ab^2)

        x1 <- x0 + ab
        y1 <- y0 + bc

        arrows(x0, y0, x1, y1, ...)
      }


      if (!is.null(continous_cols)){
        cont_pal <- continous_cols
      } else {
        cont_pal  <- RColorBrewer::brewer.pal("Reds", n = 9)
      }

      if (!is.null(catagorical_cols)){
        cat_pal <- catagorical_cols
      } else {
        cat_pal <- c("black", RColorBrewer::brewer.pal("Paired", n = 10))
      }


      if (is.numeric(sm.data@meta.data[[sm_feature_plot()]])){
        sm_cols <- cont_pal[as.numeric(cut(sm.data@meta.data[[sm_feature_plot()]],breaks = 9))]
      } else {
        sm_cols <- cat_pal[as.factor(sm.data@meta.data[[sm_feature_plot()]])]
      }

      if (is.numeric(st.data@meta.data[[st_feature_plot()]])){
        st_cols <- cont_pal[as.numeric(cut(st.data@meta.data[[st_feature_plot()]],breaks = 9))]
      } else {
        st_cols <- cat_pal[as.factor(st.data@meta.data[[st_feature_plot()]])]
      }


      if (pt_st_img()){
        plot(color_matrix)
      } else{

        plot(NULL, NULL, col = "white", xlim = c(0, dim(color_matrix)[1]), ylim = c(0, dim(color_matrix)[2]), xaxt = 'n', yaxt = 'n', ann = FALSE, bty = "n")
        #plot(fixed.scatter[, 1], fixed.scatter[, 2], col = st_cols, pch = as.numeric(pt_shape()), cex = pt_size_target(), xlim = c(0, dim(color_matrix)[1]), ylim = c(0, dim(color_matrix)[2]), xaxt = 'n', yaxt = 'n', ann = FALSE, bty = "n")
      }

      if (pt_st_points()){
        points(fixed.scatter[, 1], fixed.scatter[, 2], col = st_cols, pch = as.numeric(pt_shape()), cex = pt_size_target())
      }

      if (pt_sm_points()){
        points(coords.t[, 1], coords.t[, 2], col = sm_cols, pch = as.numeric(pt_shape()), cex = pt_size())
        arrows.1(x0 = center[1], y0 = center[2], angle.ar = stretch_angle1(), length.ar = 100*stretch_factor1(), lwd = 4, col = "blue")
        arrows.1(x0 = center[1], y0 = center[2], angle.ar = 90 + stretch_angle2(), length.ar = 100*stretch_factor2(), lwd = 4, col = "red")
      }

    }, height = 800, width = 800)

    scatter.coords <- eventReactive(align.index, {
      reset("angle"); reset("shift_x"); reset("shift_y"); reset("flip_x"); reset("flip_y"); reset("stretch_factor1"); reset("stretch_factor2"); reset("stretch_angle1"); reset("stretch_angle2")
      if (!is.null(counter)) {
        scatters[[counter]] <<- coords.ls[c(1, 2)]
        if (!is.null(tr.matrices[[counter]])) {
          tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
        } else {
          tr.matrices[[counter]] <<- coords.ls[[3]]
        }
      }
      scatter <- scatters[[as.numeric(align.index)]]$scatter
      coords <- scatters[[as.numeric(align.index)]]$coords
      counter <<- as.numeric(align.index)
      return(list(scatter, coords))
    })

    observe({
      if(input$myBtn > 0){
        if (!is.null(counter)) {
          scatters[[counter]] <<- coords.ls[c(1, 2)]
          if (!is.null(tr.matrices[[counter]])) {
            tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
            cat("Sample:", counter, "\n",  tr.matrices[[counter]][1, ], "\n", tr.matrices[[counter]][2, ], "\n", tr.matrices[[counter]][3, ], "\n\n")
          } else {
            tr.matrices[[counter]] <<- coords.ls[[3]]
          }
        }
        stopApp(tr.matrices)
      }
    })

    observeEvent(input$info, {
      showModal(modalDialog(
        title = "Instructions",
        HTML(
          "The alignment interface is modified from [STUtility](https://github.com/jbergenstrahle/STUtility/tree/master) to provide an interface for manually aligning SpaMTP (Seurat) Objects <br>",
             "This interface can be used for: <br>",
             "1. Aligning SM data to ST data coordinates.<br><br>",
             "2. Aligning SM data to a H&E image <br>",
             "<br>",
             "How to use: <br>",
             "- Adjust the coordinates of the SM data by changing the rotation, x-axis or y-axis position to match the provided reference dataset. <br>",
             "- Stretch the SM coordinates in either the direction matching the blue and/or red arrow. The angle of the arrows can also be change to match the required stretching direction. <br>",
             "- Once them SM data is aligned select the 'Return Aligned Data' button to generate a SpaMTP Seurat object with the adjusted coordinates. <br>",
             "<br>",
             "Note: If using in 'AlignSpatialOmics' mode then a SpaMTP object will be returned containing only the SM data with the ajusted coordinates. <br>",
             "For mapping SM data to corresponding ST spots, please run 'MapSpatialOmics()' with the original ST object and now updated SM object. <br>"
          ),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  }

  # Returned transformation matrices
  alignment.matrices <- runApp(list(ui = ui, server = server))
  alignment.matrices <- lapply(alignment.matrices, function(tr) {
    tr <- solve(tr)
    return(tr)
  })

  if (verbose) cat(paste("Finished image alignment. \n\n"))
  processed.ids <- which(unlist(lapply(alignment.matrices, function(tr) {!all(tr == diag(c(1, 1, 1)))})))

  # Raise error if none of the samples were processed
  if (length(processed.ids) == 0) stop("None of the samples were processed", call. = FALSE)

  # Obtain alignment matrix
  tr <- alignment.matrices[[processed.ids]]
  transformations[[processed.ids]] <- tr%*%transformations[[processed.ids]]

  map.rot.backward <- generate.map.affine(tr)
  map.rot.forward <- generate.map.affine(tr, forward = TRUE)


  # Warp pixels
  if (verbose) cat(paste0("Warping pixel coordinates for SM sample", " ... \n"))
  warped_xy <- sapply(setNames(as.data.frame(do.call(cbind, map.rot.forward(coords2$pixel_x, coords2$pixel_y))), nm = c("warped_x", "warped_y")), round, digits = 1)

  warped_mtx <- as.matrix(warped_xy)

  sm.data[[fov]][["centroids"]]@coords[,"x"] <- warped_mtx[,"warped_x"]  / st.data@images[[image.slice]]@scale.factors[[image.res]]
  sm.data[[fov]][["centroids"]]@coords[,"y"] <- warped_mtx[,"warped_y"]  / st.data@images[[image.slice]]@scale.factors[[image.res]]

  return(sm.data)
}

rotate <- function (
    angle,
    center.cur
) {
  alpha <- 2*pi*(angle/360)
  #center.cur <- c(center.cur, 0)
  #points(center.cur[1], center.cur[2], col = "red")
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.transf(center.cur[1], center.cur[2], alpha)%*%tr
  return(tr)
}


#' Creates a transformation matrix that translates an object
#' in 2D
#'
#' @param translate.x,translate.y translation of x, y coordinates

translate <- function (
    translate.x,
    translate.y
) {
  tr <- rigid.transl(translate.x, translate.y)
  return(tr)
}


#' Creates a transformation matrix that mirrors an object
#' in 2D along either the x axis or y axis around its
#' center of mass
#'
#' @param mirror.x,mirror.y Logical specifying whether or not an
#' object should be reflected
#' @param center.cur Coordinates of the current center of mass
#'

mirror <- function (
    mirror.x = FALSE,
    mirror.y = FALSE,
    center.cur
) {
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.refl(mirror.x, mirror.y)%*%tr
  tr <- rigid.transl(center.cur[1], center.cur[2])%*%tr
  return(tr)
}


#' Stretch along angle
#'
#' Creates a transformation matrix that stretches an object
#' along a specific axis
#'
#' @param r stretching factor
#' @param alpha angle
#' @param center.cur Coordinates of the current center of mass
#'

stretch <- function(r, alpha, center.cur) {
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.rot(alpha, forward = TRUE)%*%tr
  tr <- rigid.stretch(r)%*%tr
  tr <- rigid.rot(alpha, forward = FALSE)%*%tr
  tr <- rigid.transl(center.cur[1], center.cur[2])%*%tr
  return(tr)
}


#' Creates a transformation matrix for rotation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#'
#' @param alpha rotation angle
#' @param forward should the rotation be done in forward direction?
#'

rigid.rot <- function (
    alpha = 0,
    forward = TRUE
) {
  alpha <- 2*pi*(alpha/360)
  tr <- matrix(c(cos(alpha), ifelse(forward, -sin(alpha), sin(alpha)), 0, ifelse(forward, sin(alpha), -sin(alpha)), cos(alpha), 0, 0, 0, 1), nrow = 3)
  return(tr)
}


#' Creates a transformation matrix for rotation and translation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#' followed by a translation with an offset of (h, k). Points are assumed to be
#' centered at (0, 0).
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#' @param alpha rotation angle
#'

rigid.transf <- function (
    h = 0,
    k = 0,
    alpha = 0
) {
  tr <- matrix(c(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, h, k, 1), nrow = 3)
  return(tr)
}

#' Creates a transformation matrix for translation with an offset of (h, k)
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#'

rigid.transl <- function (
    h = 0,
    k = 0
) {
  tr <-  matrix(c(1, 0, 0, 0, 1, 0, h, k, 1), nrow = 3)
  return(tr)
}

#' Creates a transformation matrix for reflection
#'
#' Creates a transformation matrix for reflection where mirror.x will reflect the
#' points along the x axis and mirror.y will reflect thepoints along the y axis.
#' Points are assumed to be centered at (0, 0)
#'
#' @param mirror.x,mirror.y Logical: mirrors x or y axis if set to TRUE

rigid.refl <- function (
    mirror.x,
    mirror.y
) {
  tr <- diag(c(1, 1, 1))
  if (mirror.x) {
    tr[1, 1] <- - tr[1, 1]
  }
  if (mirror.y) {
    tr[2, 2] <- - tr[2, 2]
  }
  return(tr)
}

#' Creates a transformation matrix for stretching
#'
#' Creates a transformation matrix for stretching by a factor of r
#' along the x axis.
#'
#' @param r stretching factor

rigid.stretch <- function (
    r
) {
  tr <- matrix(c(r, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
}


#' Combines rigid tranformation matrices
#'
#' Combines rigid tranformation matrices in the following order:
#' translation of points to origin (0, 0) -> reflection of points
#' -> rotation by alpha degrees and translation of points to new center
#'
#' @param center.cur (x, y) image pixel coordinates specifying the current center of the tissue (stored in slot "tools" as "centers")
#' @param center.new (x, y) image pixel coordinates specifying the new center (image center)
#' @param alpha Rotation angle
#'
#' @inheritParams rigid.transf
#' @inheritParams rigid.transl
#' @inheritParams rigid.refl
#'
#' @examples
#' \dontrun{
#' library(imager)
#' library(tidyverse)
#' im <- load.image("https://upload.wikimedia.org/wikipedia/commons/thumb/f/fd/Aster_Tataricus.JPG/1024px-Aster_Tataricus.JPG")
#' d <- sRGBtoLab(im) %>% as.data.frame(wide="c")%>%
#'   dplyr::select(-x,-y)
#'
#' km <- kmeans(d, 2)
#'
#' # Run a segmentation to extract flower
#' seg <- as.cimg(abs(km$cluster - 2), dim = c(dim(im)[1:2], 1, 1))
#' plot(seg); highlight(seg == 1)
#'
#' # Detect edges
#' dx <- imgradient(seg, "x")
#' dy <- imgradient(seg, "y")
#' grad.mag <- sqrt(dx^2 + dy^2)
#' plot(grad.mag)
#'
#' # Extract points at edges
#' edges.px <- which(grad.mag > max(grad.mag[, , 1, 1])/2, arr.ind = TRUE)
#' points(edges.px, col = "green", cex = 0.1)
#'
#' # Apply transformations to point set
#' tr1 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(1200, 1200), alpha = 90)
#' tr2 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(500, 1200), mirror.x = T, alpha = 30)
#' tr3 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(1200, 500), mirror.y = T, alpha = 270)
#' plot(edges.px, xlim = c(0, 1700), ylim = c(0, 1700), cex = 0.1)
#' points(t(tr1%*%t(edges.px[, 1:3])), cex = 0.1, col = "red")
#' points(t(tr2%*%t(edges.px[, 1:3])), cex = 0.1, col = "yellow")
#' points(t(tr3%*%t(edges.px[, 1:3])), cex = 0.1, col = "blue")
#' }
#'
#' @export

combine.tr <- function (
    center.cur,
    center.new,
    alpha,
    mirror.x = FALSE,
    mirror.y = FALSE
) {

  alpha <- 2*pi*(alpha/360)
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])

  # reflect
  tr <- rigid.refl(mirror.x, mirror.y)%*%tr

  # rotate and translate to new center
  tr <- rigid.transf(center.new[1], center.new[2], alpha)%*%tr
}



generate.map.affine <- function (
    tr,
    forward = FALSE
) {

  if (forward) {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      xy <- t(solve(tr)%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  } else {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      xy <- t(tr%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  }
  return(map.affine)
}



#' Manually align an image (e.g. H&E, Immuno) to a SM SpaMTP dataset
#'
#' @param image_path Character string defining the full path of the image to align.
#' @param SpaMTP SpaMTP Seurat object to align the image to.
#' @param fov Character string defining the image fov that contains the SM data coordinates (default = "fov").
#' @param grey.scale Numeric value defining the grey scale to use for generating a tissue mask from the provided image (default = 0.5).
#' @param plot.greyscale Boolean indicating whether to display the grey scale plot used to generate the tissue mask (default = FALSE).
#' @param seed Integer value defining the seed to use for calculating random fake gene values for aligning the image (default = 123).
#' @param n.spots Integer specifying the number of fake spots to generate for the aligned image. If NULL the number of spots will match that of the SpaMTP object provided (default = NULL).
#' @param ... Additional inputs used by the AlignSpatialOmics function. Please see documentation or call ?AlignSpatialOmics for more infomation.
#'
#' @return A SpaMTP Seurat object with the provided image aligned to the SM data. The image is stored in the `@image$slice1` slot.
#'
#' @export
#'
#' @examples
#' # AddSMImage(image_path = "../HnE_image.png", SpaMTP = SpaMTP_obj)
AddSMImage <- function(image_path, SpaMTP, fov = "fov", grey.scale = 0.5, plot.greyscale = FALSE, seed = 123, n.spots = NULL, ...) {

  ## Load Image
  img <- magick::image_read(image_path)

  ## Convert Image to Array
  img_data <- magick::image_data(img)
  img_array <- as.integer(img_data)/255

  ## Get Tissue mask for fake gene counts matching tissue feature
  gray_img <- EBImage::channel(img_array, "gray")

  tissue_mask <- gray_img < grey.scale  # Adjust this threshold based on your image

  # Get tissue pixel coordinates
  tissue_coords <- which(tissue_mask, arr.ind = TRUE)

  if (plot.greyscale) {
    numeric_tissue_mask <- as.numeric(tissue_mask)

    # Reshape it back to matrix form
    numeric_tissue_mask <- matrix(numeric_tissue_mask, nrow = nrow(tissue_mask), ncol = ncol(tissue_mask))
    EBImage::display(numeric_tissue_mask, method = "raster", title = "Tissue Mask", all = TRUE)
  }

  ## Make Fake Seurat Object Coordinates
  image_height <- dim(img_array)[1]
  image_width <- dim(img_array)[2]

  set.seed(seed)  # for reproducibility
  if (is.null(n.spots)){
    n_spots <- dim(SpaMTP)[2]
  } else {
    n_spots <- as.numeric(n.spots)
  }

  sampled_coords <- tissue_coords[sample(1:nrow(tissue_coords), n_spots), ]
  x_coords <- sampled_coords[, 2]  # X coordinates (column index)
  y_coords <- sampled_coords[, 1]  # Y coordinates (row index)


  ## Create fake gene expression data
  n_genes <- 200  # Number of genes
  fake_expr_matrix <- matrix(rnorm(n_spots * n_genes), nrow = n_genes, ncol = n_spots)

  ## Create a data frame with spatial coordinates
  metadata <- data.frame(
    x = x_coords,
    y = y_coords,
    row.names = paste0("spot_", 1:n_spots)
  )

  ## Create Seurat Object
  seurat_obj <- Seurat::CreateSeuratObject(counts = fake_expr_matrix)
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata)

  ## Add Spatial Coordinates
  seurat_obj@meta.data$cell <- rownames(seurat_obj@meta.data)

  fake_coords <- seurat_obj@meta.data
  fake_coords$imagerow <- fake_coords$y
  fake_coords$imagecol <- fake_coords$x
  rownames(fake_coords) <- fake_coords$cell
  fake_coords$x <- fake_coords$x - min(fake_coords$x)
  fake_coords$y <- fake_coords$y - min(fake_coords$y)
  fake_fov <- CreateFOV(
    fake_coords[, c("imagerow", "imagecol")],
    type = "centroids",
    radius = SpaMTP@images[[fov]]$centroids@radius,
    assay = "Spatial",
    key = Key("slice", quiet = TRUE)
  )

  scale.factors <- Seurat::scalefactors(spot = 1, fiducial = 30, hires = 1, lowres = 1)

  ## Add Image to Seurat Object
  visium.fov <- new(
    Class = "VisiumV2",
    boundaries = fake_fov@boundaries,
    molecules = fake_fov@molecules,
    assay = fake_fov@assay,
    key = fake_fov@key,
    image = img_array,
    scale.factors = scale.factors
  )
  seurat_obj@images[["slice1"]] <- visium.fov

  ## Manually Align H&E Image
  aligned_SpaMTP <- AlignSpatialOmics(sm.data = SpaMTP, st.data = seurat_obj, image.slice = "slice1", fov = fov, ...)

  real_coords <- Seurat::GetTissueCoordinates(aligned_SpaMTP, image = fov)
  real_coords$imagerow <- real_coords$x
  real_coords$imagecol <- real_coords$y
  rownames(real_coords) <- real_coords$cell

  real_fov <- SeuratObject::CreateFOV(
    real_coords[, c("imagerow", "imagecol")],
    type = "centroids",
    radius = aligned_SpaMTP@images[[fov]]$centroids@radius,
    assay = "Spatial",
    key = SeuratObject::Key("slice", quiet = TRUE)
  )

  new.visium.fov <- methods::new(
    Class = "VisiumV2",
    boundaries = real_fov@boundaries,
    molecules = real_fov@molecules,
    assay = real_fov@assay,
    key = real_fov@key,
    image = img_array,
    scale.factors = scale.factors
  )

  aligned_SpaMTP@images[["slice1"]] <- new.visium.fov

  return(aligned_SpaMTP)
}

