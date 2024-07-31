library(Seurat)
library(RANN)
library(stats)
library(utils)
library(graphics)
library(imager)
library(shiny)
library(shinyjs)
library(magrittr)
library(zeallot)
library(RColorBrewer)



#### SpaMTP MALDI to Visium Spot Merging Functions  ####################################################################################################################################################################


#' Splits a pixel into multiple spots
#'     - This function takes a central point and a pixel radius and splits the pixel into multiple smaller spots based on the specified configuration.
#'
#' @param center_point A numeric vector of length 2 representing the (x, y) coordinates of the center point of the pixel.
#' @param pixel_radius The radius of the original pixel.
#' @param pseudo_n The number of smaller spots to split the pixel into. Must be either 4 or 9 for square pixels and 4 for a circular pixel.
#' @param pixel_shape A character string indicating the shape of the original pixel. Options are "square" (default) or "circle".
#' @param show_split_diagram A logical indicating whether to plot a diagram showing the original pixel and the split spots (default = FALSE).
#'
#' @returns A matrix with two columns representing the (x, y) coordinates of the centers of the smaller spots after splitting the pixel.
#' @export
#'
#' @examples
#' # Split a square pixel into 4 spots
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 4, pixel_shape = "circle", show_split_diagram = TRUE)
#'
#' # Split a circular pixel into 9 spots and display the split diagram
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 9, pixel_shape = "square", show_split_diagram = TRUE)
#'
#' # Split a square pixel into 16 spots and display the split diagram
#' split_pixel(center_point = c(0, 0), pixel_radius = 1, pseudo_n = 16, pixel_shape = "square", show_split_diagram = TRUE)
split_pixel <- function(center_point, pixel_radius, pseudo_n = 4, pixel_shape = "square", show_split_diagram = FALSE) {

  if (pseudo_n == 4 || ((pseudo_n == 9 || pseudo_n == 16) && pixel_shape == "square")) {

    x <- center_point[1]
    y <- center_point[2]

    if (pixel_shape == "circle") {

      R <- pixel_radius
      # Radius of each smaller circle
      r <- R * (sqrt(pseudo_n) - 1)
      dis <- r * sqrt(pseudo_n)
      # Coordinates of the centers of the smaller circles
      centers <- matrix(c(R - dis + x, R - dis + y, R - dis + x, -(R - dis) + y,
                          -(R - dis) + x, R - dis + y, -(R - dis) + x, -(R - dis) + y), ncol = 2, byrow = TRUE)

    } else {

      r <- (pixel_radius * 2) / (sqrt(pseudo_n) * 2)

      if (pseudo_n == 4) {
        centers <- matrix(c(x - r, y - r, x + r, y - r, x + r, y + r, x - r, y + r), ncol = 2, byrow = TRUE)
      } else if (pseudo_n == 9) {
        centers <- matrix(c(x - r * 0, y + r * 0, x - r * 0, y - r * 2, x - r * 2, y + r * 0,
                            x - r * 2, y + r * 2, x - r * 2, y - r * 2,
                            x + r * 2, y + r * 0, x + r * 0, y + r * 2,
                            x + r * 2, y - r * 2, x + r * 2, y + r * 2), ncol = 2, byrow = TRUE)
      } else {
        centers <- matrix(c(x - r * 1, y + r * 1, x + r * 1, y - r * 1, x + r * 1, y + r * 1,
                            x - r * 1, y - r * 1, x - r * 1, y - r * 3,
                            x - r * 1, y + r * 3, x - r * 3, y + r * 1,
                            x - r * 3, y - r * 1, x - r * 3, y - r * 3,
                            x - r * 3, y + r * 3, x + r * 1, y + r * 3,
                            x + r * 1, y - r * 3, x + r * 3, y + r * 1,
                            x + r * 3, y - r * 1, x + r * 3, y - r * 3,
                            x + r * 3, y + r * 3), ncol = 2, byrow = TRUE)
      }
    }



    if (show_split_diagram) {
      if (pixel_shape == "circle") {
        plot(c(x - pixel_radius * 2, x + pixel_radius * 2), c(y - pixel_radius * 2, y + pixel_radius * 2), type = "n", xlab = "X-axis", ylab = "Y-axis")
        graphics::symbols(centers[, 1], centers[, 2], circles = rep(r, pseudo_n), add = TRUE, inches = 1, col = "red", bg = "white")
        graphics::points(x, y, pch = 19, col = "blue")
        graphics::title(main = "Four Circles Inside a Larger Circle")

      } else {
        graphics::plot(c(x - pixel_radius * 2, x + pixel_radius * 2), c(y - pixel_radius * 2, y + pixel_radius * 2), type = "n", xlab = "X-axis", ylab = "Y-axis")
        graphics::rect(x - pixel_radius, y - pixel_radius, x + pixel_radius, y + pixel_radius, border = "blue", lty = 1, col = NA)
        graphics::symbols(centers[, 1], centers[, 2], circles = rep(r, pseudo_n), add = TRUE, inches = 1, col = "red", bg = "white")
        graphics::title(main = "Four Circles Inside a Square")
      }
    }

    return(centers)

  } else {
    warning("ERROR: pseudo_n must be either 4 or 9 for square pixels and 4 for circle")
    stop("n invalid integer")
  }
}



#' Increase the resolution of Spatial Metabolomics data for merging with Spatial Transcriptomic (Visium) data
#'    - This function takes a Spatial Metabolomics dataset and increases its resolution by generating pseudo-high-resolution spots. This is useful for merging the data with Spatial Transcriptomic (Visium) data that has a different resolution.
#'    - This function uses the split_pixel() function
#'
#' @param SM.data A Seurat Spatial Metabolomics object containing spatial metabolomic information. It should have 'x_coord' and 'y_coord' columns representing spot coordinates.
#' @param res_increase An integer specifying the factor by which the resolution should be increased. res_increase = 4, will generates 4 times as many spots in each dimension (default = 4).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A data frame with increased resolution metadata, containing additional columns 'new_x_coord', 'new_y_coord', 'new_barcode', and 'old_barcode'. The 'new_x_coord' and 'new_y_coord' columns represent the coordinates of the pseudo-high-resolution spots. The 'new_barcode' column is a unique identifier for each pseudo-high-resolution spot, and 'old_barcode' retains the original spot identifier.
#'
#' @examples
#' ## Increase resolution of MALDI dataset by a factor of 4
#' # increased_data <- increase_SM_res(SeuratObj, res_increase = 4)
increase_SM_res <- function(SM.data, res_increase = 4, verbose = TRUE) {

  verbose_message(message_text = "Increasing the resolution of MALDI Pixel Data ...\n", verbose = verbose)

  # Assuming SM.data@meta.data has columns 'x_coord' and 'y_coord'

  # Randomly sample 1000 indices
  sampled_indices <- sample(1:(nrow(SM.data@meta.data) - 1), size = 1000, replace = TRUE)

  # Initialize a list to store the distances
  distances <- c()

  # Calculate distances for each pair
  for (i in sampled_indices) {
    point1 <- c(SM.data@meta.data$x_coord[i], SM.data@meta.data$y_coord[i])
    point2 <- c(SM.data@meta.data$x_coord[i + 1], SM.data@meta.data$y_coord[i + 1])
    distances <- c(distances, sqrt(sum((point1 - point2)^2)) / 2)
  }

  # Find the median distance
  median_distance <- stats::median(distances)

  verbose_message(message_text = paste0("Median Distance Between MALDI Spots: ", median_distance, "\n"), verbose = verbose)

  # Generates new obs matrix with higher resolution
  new_meta_data <- data.frame(matrix(NA, nrow = nrow(SM.data@meta.data) * 4, ncol = ncol(SM.data@meta.data)))
  colnames(new_meta_data) <- colnames(SM.data@meta.data)

  total_spots = nrow(SM.data@meta.data)

  verbose_message(message_text = "Generating psuedo-highres MALDI data: ", verbose = verbose)


  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = total_spots, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")


  for (spot_idx in 1:total_spots) {

    sub_spot_idx <- (spot_idx * 4)-3
    pseudo_centers <- split_pixel(c(SM.data@meta.data$x_coord[spot_idx], SM.data@meta.data$y_coord[spot_idx]), median_distance, pseudo_n = 4, pixel_shape = "square", show_split_diagram = FALSE)

    for (idx in seq(0, (dim(pseudo_centers)[1])-1)) {
      new_meta_data[sub_spot_idx + idx, ] <- SM.data@meta.data[spot_idx, ]
      new_meta_data[sub_spot_idx + idx, "old_barcode"] <- rownames(SM.data@meta.data[spot_idx,])
      new_meta_data[sub_spot_idx + idx, c("new_x_coord", "new_y_coord")] <- pseudo_centers[idx+1, ]
      new_meta_data[sub_spot_idx + idx, "new_barcode"] <- paste0(pseudo_centers[idx+1, 1], "_", pseudo_centers[idx+1, 2])
    }
    utils::setTxtProgressBar(pb, spot_idx)
  }
  close(pb)
  return(new_meta_data)
}



#' Generate new MALDI counts matrix for equivalent Visium spots
#'    - This function takes an original Spatial Metabolomics counts matrix, along with a metadata table (`obs_x`) containing information about the correspondence between MALDI and Visium spots. It generates a new counts matrix where each MALDI spot has an associated aggregated count based on its corresponding Visium spots.
#'
#' @param original_SM A Seurat Spatial Metabolomics object containing the original counts matrix
#' @param obs_x A metadata table with information about the correspondence between MALDI and Visium spots. It should have columns 'Visium_spot' and 'MALDI_barcodes'.
#' @param assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values.
#' @param slots Vector of character strings describing which slots to pull the relative intensity values from (default = c("counts", "data")).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A data frame representing the new counts matrix for equivalent Visium spots,where each row corresponds to a Visium spot and columns correspond to m/z features.
#' @export
#'
#' @examples
#' ## Generate new MALDI counts matrix for equivalent Visium spots
#' # new_counts <- generate_new_SM_counts(SeuratObj, obs_x, assay = "Spatial")
generate_new_SM_counts <- function(original_SM, obs_x, assay, slots, verbose = TRUE) {

  verbose_message(message_text = "Merging MALDI counts ... ", verbose = verbose)

  data_list <- list()
  new_data_list <- list()

  for (slot in slots) {
    data_list[[slot]] <- t(as.data.frame(original_SM[[assay]][slot]))
    new_data_list[[slot]] <- data.frame(matrix(NA, nrow = 0, ncol = ncol(data_list[[slot]])))
    colnames(new_data_list[[slot]]) <- colnames(data_list[[slot]])
  }


  total_spots = nrow(obs_x)
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                              max = total_spots, # Maximum value of the progress bar
                              style = 3,    # Progress bar style (also available style = 1 and style = 2)
                              width = 50,   # Progress bar width. Defaults to getOption("width")
                              char = "=")

  for (spot_idx in 1:total_spots) {
    spot_barcode <- obs_x$Visium_spot[spot_idx]
    barcode_list <- unlist(strsplit(obs_x$MALDI_barcodes[spot_idx], ", "))

    if (all(barcode_list[1] == barcode_list)) {
      for (data_slot in names(data_list)){
        new_data_list[[data_slot]][spot_barcode, ] <- data_list[[data_slot]][barcode_list[1], ]
      }
    } else {
      for (data_slot in names(data_list)){
        mini_raw_count_dfs <- lapply(barcode_list, function(MALDI_barcodes) data_list[[data_slot]][MALDI_barcodes, ])
        combined_df <- do.call(rbind, mini_raw_count_dfs)
        new_data_list[[data_slot]][spot_barcode, ] <- colMeans(combined_df)
      }
    }
    utils::setTxtProgressBar(pb, spot_idx)
  }
  close(pb)

  return(new_data_list)
}



#' Converts and aggregates Spatial Metabolomic (MALDI) data to corresponding Spatial Transcriptomics (Visium) spots.
#'    - This function uses generate_new_SM_counts() and increase_SM_res()
#'
#' @param SM.data A Seurat object representing the Spatial Metabolomics data.
#' @param ST.data A Seurat object representing the Visium Spatial Transcriptomics data.
#' @param img_res Character string defining the image resolution associated with the Visium image pixel data (default = "hires").
#' @param res_increase Integer value defining the factor by which the resolution of MALDI spots should be increased before assignment. It should be either 4 or 9, see increase_SM_res() documentation for specifics (Default = NULL).
#' @param annotations Boolean value indicating if the Spatial Metabolomics (MALDI) Seurat object contains annotations assigned to m/z values (default = FALSE).
#' @param assay Character string defining the Seurat assay that contains the annotated counts and metadata corresponding to the m/z values (default = "Spatial").
#' @param slice Character string of the image slice name in the Visium object (default = "slice1").
#' @param slots Vector of character strings describing which slots from the Spatial Metabolomic Seurat object to adjust for the new overlayed object (default = c("counts", "data")).
#' @param new_SpM.assay Character string defining the assay name of the new overlayed Seurat object containing all updated metabolomic data (default = "SPM").
#' @param add.ST Boolean indicating whether to add the Spatial Metabolomic data from the Seurat Visium Object to the new updated Seruat Object (default = TRUE).
#' @param ST.assay Character string specifying the current assay to use to extract transcriptional data from (default = "Spatial").
#' @param ST.layers Vector of character strings defining the relative slots to extract from the initial Visium object, to add to the new Seurat Object if ST.assay == TRUE (default = c("counts", "data")).
#' @param new_SpT.assay Character string defining the assay name of the new overlayed Seurat object containing all updated transcriptomics data (default = "SPT").
#' @param verbose Boolean value indicating whether to print pregression update messages and progress bar (default = TRUE).
#'
#' @return A Seurat object with the Spatial Metabolomic data assigned to equivalent Spatial Transcripomics (Visium) spots.
#' @export
#'
#' @examples
#' ## Convert MALDI data to equivalent Visium spots
#' # MapSpatialOmics(VisiumObj, SeuratObj, img_res = "hires", new_library_id = "MALDI", res_increase = NULL)
MapSpatialOmics <- function(SM.data, ST.data, res_increase = NULL, annotations = FALSE, assay = "Spatial", slots = c("counts"), img_res = "hires", slice = "slice1", new_SpM.assay = "SPM", add.ST = TRUE, ST.assay = "Spatial", ST.layers = c("counts"), new_SpT.assay = "SPT", verbose = TRUE) {

  SM.coords <- GetTissueCoordinates(SM.data)

  scale.factor <- ST.data@images[[slice]]@scale.factors[[img_res]]

  SM.metadata <- SM.data@meta.data
  SM.metadata$x_coord <- SM.coords[,"x"] * scale.factor
  SM.metadata$y_coord <- SM.coords[,"y"] * scale.factor

  SM.data@meta.data[c("x_coord", "y_coord")] <- SM.metadata[c("x_coord", "y_coord")]

  new_SM_metadata <- SM.metadata

  ## Increases resolution of the SM data (This will shift the centroid position close to each Visium spot)
  if (!is.null(res_increase)) {
    if (res_increase == 4 || res_increase == 9) {
      new_SM_metadata <- increase_SM_res(SM.data, res_increase = 4, verbose = verbose)
    } else {
      stop("Error: res_increase must be either 4 or 9\n")
    }
  } else {
    new_SM_metadata$new_x_coord <- new_SM_metadata$x_coord
    new_SM_metadata$new_y_coord <- new_SM_metadata$y_coord
    new_SM_metadata$old_barcode <- rownames(new_SM_metadata)
  }

  verbose_message(message_text = "Assigning MALDI to Visium Spots ... \n", verbose = verbose)

  ## Get coordinates of ST data
  image_data <- data.frame(ST.data@images[[slice]]$centroids@coords) * scale.factor

  ## gets radius of visium spot in pixels size
  radius <- ST.data@images[[slice]]$centroids@radius * scale.factor

  ## find which pixels fall into each spot radius
  new_coords <- as.matrix(new_SM_metadata[, c("new_x_coord", "new_y_coord")])
  query_coords <- as.matrix(image_data[, c("x", "y")])
  v_points <- RANN::nn2(new_coords,query_coords, treetype = "kd",searchtype = "radius", radius = radius)$nn.idx


  ## Constructing a df to store new metadata based on binned pixels
  obs_ <- data.frame(
    index = rownames(new_SM_metadata),
    nFeature_Spatial = new_SM_metadata[[paste0("nFeature_",assay)]] ,
    old_barcode = new_SM_metadata$old_barcode,
    Visium_spot = "Not_assigned")


  ## Generating a new df which contains the corresponding visium spot for each SM pixel
  new_df <- lapply(1:nrow(v_points), function(i){
    lapply(v_points[i,], function(x){
      if (as.numeric(x) != 0){
        if (obs_[x,"Visium_spot"] == "Not_assigned"){
          obs_[x,"Visium_spot"] <- rownames(ST.data@meta.data)[i]
          obs_[x,]
        } else {
          obs_[length(obs_$Visium_spot)+1,] <- obs_[x,]
          obs_[length(obs_$Visium_spot)+1,] <- rownames(ST.data@meta.data)[i]
          obs_[length(obs_$Visium_spot)+1,]
        }
      }
    })
  })

  obs_ <- dplyr::bind_rows(
    new_df
  )

  ## Tidying up df to create Seurat object
  obs_x <- obs_[, c("Visium_spot", paste0("nFeature_",assay), "old_barcode")]

  obs_x <- obs_x %>%
    dplyr::group_by(Visium_spot) %>%
    dplyr::summarize(MALDI_barcodes = toString(unique(old_barcode)))

  obs_x <- data.frame(obs_x)
  counts_x <- generate_new_SM_counts(SM.data, obs_x, assay = assay, slots = slots, verbose = verbose)

  verbose_message(message_text = "Generating new MALDI Seurat Object ... ", verbose = verbose)

  seuratobj <- Seurat::CreateSeuratObject(counts = t(counts_x[[names(counts_x)[1]]]), assay = new_SpM.assay) #creates a new assay with the spatial metabolomics

  if (names(counts_x)[1] != "counts"){
    seuratobj[[new_SpM.assay]][names(counts_x)[1]] <- seuratobj[[new_SpM.assay]]
    seuratobj[[new_SpM.assay]]$counts <- NULL
  }

  if (length(names(counts_x)) > 1){
    for (slot_name in names(counts_x)[2:length(names(counts_x))]){
      seuratobj[[new_SpM.assay]][slot_name] <- t(counts_x[[slot_name]])
    }
  }


  rownames(obs_x) <- rownames(seuratobj[[new_SpM.assay]]@cells)
  seuratobj@meta.data[colnames(obs_x)] <- obs_x

  ## add additional visium meta.data
  seuratobj@meta.data[colnames(ST.data@meta.data)] <- ST.data@meta.data[colnames(ST.data@meta.data)][rownames(ST.data@meta.data) %in% rownames(seuratobj@meta.data),]

  vis_subset <- subset(ST.data, cells = rownames(seuratobj@meta.data))
  seuratobj <- subset(seuratobj, cells = rownames(vis_subset@meta.data))
  seuratobj@images[["slice1"]] <- vis_subset@images[[slice]]



  if (annotations){
    seuratobj[[new_SpM.assay]]@meta.data <- SM.data[[assay]]@meta.data
  }

  if (add.ST){
    visium_sub <- subset(ST.data, cells = intersect(colnames(seuratobj), colnames(ST.data)))

    # add Transcriptomic assay to combined object
    counts_layer <- NULL
    data_layer <- NULL

    for (layer in ST.layers) {
      if ( layer == "counts" ) {
        counts_layer <- visium_sub[[ST.assay]]$counts
      }  else if (layer == "data") {
        data_layer <- visium_sub[[ST.assay]]$data
      } else {
        stop("ST.layers given does not match requirments. Must be 'counts' and/or 'data'")
      }
    }
    st_assay <- CreateAssay5Object(counts = counts_layer, data = data_layer)
    seuratobj[[new_SpT.assay]] <- st_assay


  }
  return(seuratobj)

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
               column(4,  checkboxInput(inputId = "flip_x",
                                        label = "mirror along x axis",
                                        value = FALSE)
               ),
               column(4,  checkboxInput(inputId = "flip_y",
                                        label = "mirror along y axis",
                                        value = FALSE)
               ),
               column(4,  checkboxInput(inputId = "show_ST_spots",
                                        label = "show ST spots",
                                        value = FALSE)
               ),
               column(4,  checkboxInput(inputId = "show_ST_img",
                                        label = "show image",
                                        value = TRUE)
               ),
               column(4,  checkboxInput(inputId = "show_SM_spots",
                                        label = "show SM spots",
                                        value = TRUE)
               )

             ),
             selectInput(inputId = "sample", choices = (1:length(scatters))[-reference.index],

                         label = "Select sample", selected = reference.index),
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

      return(list(scatter = scatter.t, coords = coords.t, tr = tr, xylimits = image.dims[[input$sample]]))
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

    scatter.coords <- eventReactive(input$sample, {
      reset("angle"); reset("shift_x"); reset("shift_y"); reset("flip_x"); reset("flip_y"); reset("stretch_factor1"); reset("stretch_factor2"); reset("stretch_angle1"); reset("stretch_angle2")
      if (!is.null(counter)) {
        scatters[[counter]] <<- coords.ls[c(1, 2)]
        if (!is.null(tr.matrices[[counter]])) {
          tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
        } else {
          tr.matrices[[counter]] <<- coords.ls[[3]]
        }
      }
      scatter <- scatters[[as.numeric(input$sample)]]$scatter
      coords <- scatters[[as.numeric(input$sample)]]$coords
      counter <<- as.numeric(input$sample)
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
        HTML("The selected sample is highlighted by its coordinates under the tissue <br>",
             "highlighted in red. only rigid transformations are allowed, meaning <br>",
             "rotation, shifts along x/y-axes and reflections.<br><br>",
             "1. Select sample that you want to align to a reference [default: 2]<br>",
             "2. Adjust transformation parameters to fit the sample image to the reference<br>",
             "3. Repeat 1-4 until all samples are aligned<br>",
             "4. Press the 'return aligned data' button to return results"),
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


