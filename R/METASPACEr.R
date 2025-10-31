#### SpaMTP METASPACE R CLIENT #################################################################################################################################################################################


#' @title Metaspace R Client Constructor
#'
#' @description Creates a client object with methods to interact with the METASPACE GraphQL API.
#'
#' @param host Character string specifying the METASPACE host URL (default = "https://metaspace2020.org").
#' @param api_key Optional. Your METASPACE API key for accessing private datasets. If accessing a public dataset this can be left as `NULL` (default = NULL).
#'
#' @return A list containing functions for querying METASPACE data.
#'
#' @importFrom httr add_headers POST GET status_code content
#' @importFrom jsonlite toJSON
#' @importFrom rlang %||%
#'
#' @export
#'
#' @examples
#' # ms <- metaspace_client()
metaspace_client <- function(host = "https://metaspace2020.org", api_key = NULL) {
  graphql_url <- paste0(host, "/graphql")


  ##  GraphQL executor
  query_graphql <- function(query, variables = list()) {
    hdr <- list(`Content-Type` = "application/json", `Source` = "api")
    if (!is.null(api_key)) hdr$Authorization <- paste0("Api-Key ", api_key)
    headers <- do.call(add_headers, hdr)

    resp <- POST(graphql_url,
                 body = toJSON(list(query = query, variables = variables), auto_unbox = TRUE),
                 headers, encode = "json")
    if (status_code(resp) != 200) stop("HTTP ", status_code(resp))
    res <- content(resp, as = "parsed")
    if (!is.null(res$errors)) {
      msg <- res$errors[[1]]$message
      if (status_code(resp) == 401)
        stop("Authentication required – get your key at https://metaspace2020.org/user/me")
      stop("GraphQL error: ", msg)
    }
    res$data
  }


  ##  1. annotation table (with pagination)
  get_dataset_results <- function(dataset_id, fdr = 0.1,
                                  database = c("HMDB", "v4"),
                                  include_chem_mods = FALSE,
                                  include_neutral_losses = FALSE,
                                  page_size = 100) {
    # ---- get DB ID
    db_res <- query_graphql('query { allMolecularDBs { id name version } }')
    db_id <- NULL
    for (db in db_res$allMolecularDBs)
      if (db$name == database[1] && db$version == database[2]) { db_id <- db$id; break }
    if (is.null(db_id)) stop("Database not found: ", paste(database, collapse = " "))

    # ---- construct filter
    filter <- list(fdrLevel = fdr, databaseId = db_id)
    if (!include_chem_mods) filter$hasChemMod <- FALSE
    if (!include_neutral_losses) filter$hasNeutralLoss <- FALSE

    # ---- GraphQL query with pagination
    q <- '
    query($filter: AnnotationFilter, $dFilter: DatasetFilter, $limit: Int, $offset: Int) {
      allAnnotations(filter: $filter, datasetFilter: $dFilter, limit: $limit, offset: $offset) {
        sumFormula
        adduct
        neutralLoss
        chemMod
        mz
        msmScore
        rhoSpatial
        rhoSpectral
        rhoChaos
        fdrLevel
        offSample
        possibleCompounds { name }
        isotopeImages { mz url minIntensity maxIntensity totalIntensity }
      }
    }
  '

    ann_list <- list()
    offset <- 0
    repeat {
      vars <- list(filter = filter, dFilter = list(ids = dataset_id),
                   limit = page_size, offset = offset)
      res <- query_graphql(q, vars)

      anns <- res$allAnnotations
      if (length(anns) == 0) break  # no more annotations

      # ---- convert annotations to data.frame
      page_df <- lapply(anns, function(a) {
        names_str <- paste(sapply(a$possibleCompounds, `[[`, "name"), collapse = "; ")
        int <- if (length(a$isotopeImages)) a$isotopeImages[[1]]$maxIntensity else NA
        data.frame(
          formula = a$sumFormula, adduct = a$adduct,
          chemMod = ifelse(is.null(a$chemMod), "", a$chemMod),
          neutralLoss = ifelse(is.null(a$neutralLoss), "", a$neutralLoss),
          mz = a$mz, msm = a$msmScore, fdr = a$fdrLevel,
          rhoSpatial = a$rhoSpatial, rhoSpectral = a$rhoSpectral,
          moc = a$rhoChaos, offSample = a$offSample,
          intensity = int, moleculeNames = names_str,
          isotopeImages = list(a$isotopeImages),
          stringsAsFactors = FALSE
        )
      })
      ann_list <- c(ann_list, page_df)

      offset <- offset + length(anns)
    }

    if (length(ann_list) == 0) {
      warning("No annotations found")
      return(data.frame())
    }

    do.call(rbind, ann_list)
  }



  ##  2. dataset metadata
  get_dataset_info <- function(dataset_id) {
    q <- '
      query($id:String!){dataset(id:$id){
        id name status uploadDT polarity ionisationSource
        analyzer{type resolvingPower(mz:400)}
        organism organismPart condition maldiMatrix isPublic
        adducts submitter{id name} group{id name}
        databases{id name version}
      }}'
    res <- query_graphql(q, list(id = dataset_id))
    if (is.null(res$dataset)) stop("Dataset not found")
    res$dataset
  }


  ##  3. PNG dimension helper
  .png_dim <- function(url) {
    if (status_code(HEAD(url)) != 200) stop("Cannot reach image")
    raw <- GET(url, config = list(range = "bytes=0-23"))$content
    w <- readBin(raw[17:20], "integer", size = 4, endian = "big")
    h <- readBin(raw[21:24], "integer", size = 4, endian = "big")
    c(width = w, height = h)
  }


  ##  4. robust pixel extractor (MODIFIED TO RETURN FULL MATRIX)
  get_ion_image_pixels <- function(dataset_id, sum_formula, adduct,
                                   fdr = 0.1, database = c("HMDB","v4"),
                                   chem_mod = NULL, neutral_loss = NULL,
                                   isotope_idx = 1, relative = TRUE,
                                   hotspot_clipping = FALSE) {
    # ---- DB id
    db_res <- query_graphql('query{allMolecularDBs{id name version}}')
    db_id  <- NULL
    for (db in db_res$allMolecularDBs)
      if (db$name == database[1] && db$version == database[2]) { db_id <- db$id; break }
    if (is.null(db_id)) stop("Database not found")

    # ---- filter
    filter <- list(sumFormula = sum_formula, adduct = adduct,
                   databaseId = db_id, fdrLevel = fdr)
    if (!is.null(chem_mod) && chem_mod != "")   filter$chemMod <- chem_mod
    if (!is.null(neutral_loss) && neutral_loss != "") filter$neutralLoss <- neutral_loss

    q <- '
      query($filter:AnnotationFilter!,$dFilter:DatasetFilter){
        allAnnotations(filter:$filter,datasetFilter:$dFilter){
          isotopeImages{url mz maxIntensity minIntensity totalIntensity}
        }}'
    vars <- list(filter = filter, dFilter = list(ids = dataset_id))
    res  <- query_graphql(q, vars)
    if (!length(res$allAnnotations)) stop("Annotation not found")
    img_info <- res$allAnnotations[[1]]$isotopeImages[[isotope_idx]]
    if (is.null(img_info$url)) stop("No image URL")

    # ---- download PNG
    img_raw <- content(GET(img_info$url), "raw")
    img     <- readPNG(img_raw)
    nchan   <- dim(img)[3]

    # ---- intensity channel
    # R PNG reads array as [row, col, channel]
    if (nchan %in% c(1, 2)) {
      intensity_data <- img[,,1]
    } else if (nchan %in% c(3, 4)) {
      intensity_data <- img[,,1]
    } else stop("Unsupported PNG: ", nchan, " channels")

    # ---- scale to real intensity
    lo <- img_info$minIntensity %||% 0
    hi <- img_info$maxIntensity %||% 1

    # Perform scaling (0-1 normalized to min/max intensity)
    intensities <- lo + intensity_data * (hi - lo)

    # ---- alpha mask (Apply mask to set non-data points to 0)
    if (nchan %in% c(2, 4)) {
      mask <- img[,,nchan]
      intensities[mask == 0] <- 0
    }

    # ---- relative / absolute (Normalization)
    if (relative && hi > 0) intensities <- intensities / hi
    # If relative=FALSE, intensities are absolute (raw counts).

    # ---- hotspot clipping (Always apply to the final intensity values)
    if (hotspot_clipping && hi > 0) { # Only clip if there's actual signal
      # Only consider non-zero pixels for quantile calculation
      positive_intensities <- intensities[intensities > 0]
      if (length(positive_intensities) > 0) {
        # Using type 7, R's default, but check if Python uses a different method (e.g., linear interpolation)
        q75 <- quantile(positive_intensities, 0.75, type = 7)
        intensities[intensities > q75] <- q75
      }
    }

    # ---- RETURN FULL INTENSITY MATRIX (2D array)
    return(intensities)
  }


  ##  5. Python-style all_annotation_images
  all_annotation_images <- function(dataset_id,
                                    fdr = 0.1,
                                    database = c("HMDB", "v4"),
                                    only_first_isotope = TRUE,
                                    relative = TRUE,
                                    hotspot_clipping = FALSE,
                                    parallel = TRUE) {
    ann_df <- get_dataset_results(dataset_id, fdr = fdr, database = database)
    if (nrow(ann_df) == 0) stop("No annotations found at FDR ", fdr)

    get_single_image <- function(row) {
      isotope_idx <- if (only_first_isotope) 1 else seq_along(row$isotopeImages[[1]])

      lapply(isotope_idx, function(idx) {
        img_df <- get_ion_image_pixels(
          dataset_id  = dataset_id,
          sum_formula = row$formula,
          adduct      = row$adduct,
          fdr         = fdr,
          database    = database,
          chem_mod    = if (row$chemMod == "") NULL else row$chemMod,
          neutral_loss = if (row$neutralLoss == "") NULL else row$neutralLoss,
          isotope_idx = idx,
          relative    = relative,
          hotspot_clipping = hotspot_clipping
        )
        if (!nrow(img_df)) return(NULL)
        img_df
      })
    }

    if (parallel) {
      plan(multisession)
      all_images <- future_lapply(seq_len(nrow(ann_df)), function(i) get_single_image(ann_df[i, , drop = FALSE]))
    } else {
      all_images <- lapply(seq_len(nrow(ann_df)), function(i) get_single_image(ann_df[i, , drop = FALSE]))
    }
    names(all_images) <- ann_df$mz
    all_images
  }

  ##  Return public API (unchanged)
  list(
    query                 = query_graphql,
    get_results           = get_dataset_results,
    get_info              = get_dataset_info,
    get_ion_image_pixels  = get_ion_image_pixels,
    all_annotation_images = all_annotation_images,
    .png_dim              = .png_dim
  )
}


#' @title Core METASPACE Data Retrieval
#'
#' @description The main user-facing function to download annotations and ion images for a METASPACE dataset.
#'
#' @param dataset_id Character string of the METASPACE dataset ID.
#' @param fdr Numeric, the FDR level threshold (default = 0.1).
#' @param database Character vector, specifying the database name and version (default = c("HMDB", "v4")).
#' @param include_images Logical, if TRUE, downloads all ion image matrices (default = TRUE).
#' @param api_key Optional character string, your METASPACE API key (default = NULL).
#' @param isotope_idx Integer, index of the isotope image to download for all annotations (1 for main peak) (default = 1).
#' @param relative Logical, if TRUE, normalizes pixel intensities (0-1 range) (default = TRUE).
#'
#' @return A list with two elements: \code{annotations} (data.frame) and \code{images} (list of 2D matrices).
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom httr HEAD status_code
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom dplyr select
#' @export
#'
#' @examples
#' # get_metaspace(dataset_id = "2020-12-07_03h16m14s", fdr = 0.1, database= c("HMDB", "v4"), relative = TRUE)
get_metaspace <- function(dataset_id,
                          fdr = 0.1,
                          database = c("HMDB", "v4"),
                          include_images = TRUE,
                          api_key = NULL,
                          isotope_idx = 1,
                          relative = TRUE) {

  sm <- metaspace_client(api_key = api_key)

  ann_df <- sm$get_results(dataset_id, fdr = fdr, database = database)
  if (nrow(ann_df) == 0) stop("No annotations at FDR ", fdr)
  ann_df <- ann_df[order(ann_df$mz), ]
  mz_names <- sprintf("%.5f", ann_df$mz)
  rownames(ann_df) <- mz_names

  images <- NULL
  if (include_images) {
    message("Downloading ", nrow(ann_df), " ion images (full matrices) …")
    mats <- vector("list", nrow(ann_df))
    names(mats) <- mz_names
    pb <- txtProgressBar(0, nrow(ann_df), style = 3)
    for (i in seq_len(nrow(ann_df))) {
      row <- ann_df[i, , drop = FALSE]
      # mats[[i]] will now be a 2D matrix of intensity values, including 0s.
      mats[[i]] <- sm$get_ion_image_pixels(
        dataset_id   = dataset_id,
        sum_formula  = row$formula,
        adduct       = row$adduct,
        fdr          = fdr,
        database     = database,
        chem_mod     = if (row$chemMod == "") NULL else row$chemMod,
        neutral_loss = if (row$neutralLoss == "") NULL else row$neutralLoss,
        isotope_idx  = isotope_idx,
        relative     = relative
      )
      setTxtProgressBar(pb, i)
    }
    close(pb)
    images <- mats
  }

  list(annotations = ann_df, images = images)

}


#' @title Convert METASPACE Images to Feature Matrix
#'
#' @description Transforms the list of 2D ion image matrices from \code{get_metaspace} into a single matrix suitable for spatial analysis packages.
#'
#' @param metaspace_data A list object returned by \code{get_metaspace}.
#' @param transform Logical. If TRUE, transposes each image matrix before flattening (e.g., to align coordinate systems) (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A numeric matrix where rows are molecular features (m/z) and columns are spatial pixels (named 'x_y').
#'
#' @export
#'
#' @examples
#' # metaspace_to_feature_matrix(mtx, transform = TRUE)
metaspace_to_feature_matrix <- function(metaspace_data, transform = FALSE, verbose = TRUE) {

  if (is.null(metaspace_data$images) || length(metaspace_data$images) == 0) {
    stop("Input 'metaspace_data' does not contain image matrices. Ensure include_images=TRUE was used.")
  }


  if (transform){
    for (i in seq_along(metaspace_data$images)) {
      metaspace_data$images[[i]] <- t(metaspace_data$images[[i]])
    }
  }

  # 1. Extract first image to get dimensions and create column names
  first_image <- metaspace_data$images[[1]]



  if (!is.matrix(first_image)) {
    stop("Image data must be a matrix (full pixel grid). Ensure get_ion_image_pixels returns a matrix.")
  }

  # Get dimensions (rows = y, cols = x)
  h <- nrow(first_image)
  w <- ncol(first_image)

  # Create column names (pixel coordinates 'x_y')
  # We iterate x from 1 to w, and y from 1 to h
  x_coords <- rep(1:w, each = h)
  y_coords <- rep(1:h, times = w)
  pixel_names <- paste0(x_coords, "_", y_coords)

  # 2. Flatten all image matrices and combine into the final matrix
  # The final matrix will be (Features x Pixels)

  # Initialize the final matrix list
  flat_images_list <- vector("list", length(metaspace_data$images))
  names(flat_images_list) <- names(metaspace_data$images)

  message("Combining ", length(metaspace_data$images), " annotations into a feature matrix...")

  for (i in seq_along(metaspace_data$images)) {
    mat <- metaspace_data$images[[i]]

    # Check for consistency
    if (nrow(mat) != h || ncol(mat) != w) {
      stop(paste("Image matrix", i, "has inconsistent dimensions."))
    }

    flat_vector <- as.vector(mat)
    flat_images_list[[i]] <- flat_vector
  }

  # Combine list of vectors into a matrix (Features x Pixels)
  feature_matrix <- do.call(rbind, flat_images_list)

  # Assign column names (pixel coordinates)
  colnames(feature_matrix) <- pixel_names

  verbose_message(message_text = paste0("Conversion complete. Matrix dimensions: ", paste(dim(feature_matrix), collapse = " rows x "), " columns"), verbose = verbose)

  return(feature_matrix)
}


#' @title Load METASPACE Data and Create Seurat Object
#'
#' @description Downloads METASPACE data, converts it to a feature matrix, and packages it into a Seurat object with spatial metadata.
#'
#' @param dataset_id Character string of the METASPACE dataset ID.
#' @param fdr Numeric, the FDR level threshold (default = 0.1).
#' @param database Character vector, specifying the database name and version (default = c("HMDB", "v4")).
#' @param api_key Optional character string, your METASPACE API key (default = NULL).
#' @param relative Logical, if TRUE, normalizes pixel intensities (0-1 range) (default = FALSE).
#' @param transform Boolean specifiying whether to transform/flip the image (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A Seurat object initialized with the METASPACE data in the 'Spatial' assay.
#'
#' @importFrom tidyr separate
#' @importFrom SeuratObject CreateCentroids CreateFOV
#' @importFrom Seurat CreateSeuratObject AddMetaData
#' @importFrom matter as.matrix
#' @importFrom png readPNG
#' @importFrom dplyr %>% select
#'
#' @export
#'
#' @examples
#' # Load_METASPACE(dataset_id = "2020-12-07_03h16m14s", fdr = 0.1, database= c("HMDB", "v4"), relative = TRUE)
Load_METASPACE <- function(dataset_id, fdr = 0.1, database = c("HMDB", "v4"), api_key= NULL, relative = FALSE, transform = FALSE, verbose = TRUE){

  verbose_message(message_text = paste0("Downloading dataset [",dataset_id, "] from METASPACE... "), verbose = verbose)

  data <- get_metaspace(dataset_id = dataset_id,fdr=fdr, database=database, include_images = TRUE,api_key= api_key,relative= relative)

  verbose_message(message_text = "Gathering metabolite intensity values per pixel ...", verbose = verbose)

  sparse_matrix <- metaspace_to_feature_matrix(data, transform = transform, verbose = verbose)
  features <- rownames(sparse_matrix)
  pixels <- colnames(sparse_matrix)

  pixel_data <- data.frame(pixel = pixels) %>% separate(col = pixel, into = c("x", "y"), sep = "_", convert = TRUE)
  pixel_data$x_coord <- pixel_data$x
  pixel_data$y_coord <- pixel_data$y

  verbose_message(message_text = "Generating Seurat Barcode Labels from Pixel Coordinates .... ", verbose = verbose)

  rownames(sparse_matrix)<- paste("mz-", features, sep = "")

  verbose_message(message_text = "Constructing Seurat Object ....", verbose = verbose)

  mat <- matter::as.matrix(sparse_matrix)

  seuratobj <- Seurat::CreateSeuratObject(mat, assay = "Spatial")

  verbose_message(message_text = "Adding Pixel Metadata ....", verbose = verbose)

  pixel_metadata <- pixel_data
  pixel_metadata$metaspace_dataset  <-  dataset_id
  pixel_metadata$fdr_threshold <- fdr
  pixel_metadata$fdr_threshold <- paste(database, collapse = "-")

  seuratobj <- Seurat::AddMetaData(seuratobj, metadata = pixel_metadata)

  verbose_message(message_text = "Creating Centroids for Spatial Seurat Object ....", verbose = verbose)

  ## Add spatial data

  cents <- SeuratObject::CreateCentroids(data.frame(x = pixel_data$x, y = pixel_data$y, cell = c(pixels)))


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

  metaspace_metadata_names <-  c("formula","adduct","chemMod","neutralLoss","mz","msm",
                                 "fdr","rhoSpatial","rhoSpectral","moc","offSample","intensity","moleculeNames")

  seuratobj[["Spatial"]]@meta.data[metaspace_metadata_names] <- data$annotations[metaspace_metadata_names]

  return(seuratobj)
}
