#' Generates PCA analysis results for a SpaMTP Seurat Object
#'
#' This function run PCA analaysis on a SpaMTP Seurat Object.
#' The user can provide a bin/resolution size to increase the bin size and reduce the dimensionality/noise of the SM dataset prior to calculating PCAs.
#'
#' @param SpaMTP SpaMTP Seurat class object that contains spatial metabolic information.
#' @param npcs is an integer value to indicated preferred number of PCs to retain (default = 30).
#' @param variance_explained_threshold Numeric value defining the explained variance threshold (default = 0.9).
#' @param assay Character string defining the SpaMTP assay to extract intensity values from (default = "SPM").
#' @param slot Character string defining the assay slot containing the intensity values (default = "counts").
#' @param show_variance_plot Boolean indicating weather to display the variance plot output by this analysis (default = FALSE).
#' @param bin_resolution Numeric value defining the resolution to use for binning m/z peaks. If set to `NULL`, no binning will be performed (default = NULL).
#' @param resolution_units Character string specifying the units of the `bin_resolution`. Either 'ppm' or 'mz' can be provided. `bin_resolution` must be provided for this parameter to be implemented (default = "ppm").
#' @param bin_method Character string defining the method to use for binning respective m/z peaks that fall within a bin. Options for this parameter can be one of "sum", "mean", "max" or "min". `bin_resolution` must be provided for this parameter to be implemented (default = "sum").
#' @param reduction.name Character string indicating the name associated with the PCA results stored in the output SpaMTP Seurat object (default = "pca").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#'
#' @return SpaMTP object with pca results stored in the
#' @export
#'
#' @examples
#' ## For running PCA on un-adjusted peak bin sizes
#' # spamtp_obj <- RunMetabolicPCA(spamtp_obj, npcs = 50)
#'
#' ## For running PCA on increased peak bin sizes
#' # spamtp_obj <- RunMetabolicPCA(spamtp_obj, npcs = 50, bin_resolution = 50)
RunMetabolicPCA <- function(SpaMTP,
                            npcs = 30,
                            variance_explained_threshold = 0.9,
                            assay = "SPM",
                            slot = "counts",
                            show_variance_plot= FALSE,
                            bin_resolution = NULL,
                            resolution_units = "ppm",
                            bin_method = "sum",
                            reduction.name = "pca",
                            verbose = TRUE)
{
  verbose_message(message_text = "Running PCA Analysis ... ", verbose = verbose)

  npcs <- as.integer(npcs)

  if(!is.null(bin_resolution)){
    verbose_message(message_text = paste0("Binned SpaMTP object to a bin size of ", bin_resolution), verbose = verbose)
    if (!resolution_units %in% c("ppm", "mz")){
      stop("Incorrect value assigned to 'resolution_units'... value must be either 'ppm' or 'mz', please change accordingly")
    }
    if ( !bin_method %in% c("sum", "mean", "max", "min")){
      stop("Incorrect value assigned to 'bin_method'... value must be either 'sum', 'mean', 'max' or 'min', please change accordingly")
    }

    data <- BinSpaMTP(data = SpaMTP,
                      resolution = bin_resolution,
                      units = resolution_units,
                      assay = assay,
                      slot = slot,
                      method = bin_method, return.only.mtx = FALSE)

    data_mtx <- data[["binned"]]["counts"]
    temp_assay <- "binned"

  } else {
    data_mtx <- SpaMTP[[assay]][slot]
    data <- SpaMTP
    temp_assay <- assay
  }

  # PCA analysis
  verbose_message(message_text = "Scaling original matrix", verbose = verbose)


  mass_matrix = Matrix::t(data_mtx)

  verbose_message(message_text = "Running the principal component analysis ... " , verbose = verbose)

  # Runing PCA

  resampled_mat_standardised = as.matrix(Matrix::t(
    Matrix::t(mass_matrix) - Matrix::colSums(mass_matrix) / nrow(mass_matrix)
  ))

  verbose_message(message_text = "Computing the covariance" , verbose = verbose)
  cov_mat <- t(resampled_mat_standardised) %*% resampled_mat_standardised / (nrow(resampled_mat_standardised) - 1)

  verbose_message(message_text = "Computing the eigenvalue/eigenvectors", verbose = verbose)
  eigen_result <- eigen(cov_mat)
  gc()
  # Extract eigenvectors and eigenvalues
  eigenvectors <- eigen_result$vectors
  eigenvalues <- eigen_result$values

  verbose_message(message_text = "Computing PCA", verbose = verbose)

  pc = pbapply::pblapply(1:npcs, function(i) {
    temp = resampled_mat_standardised[, 1] * eigenvectors[1, i]
    for (j in 2:ncol(resampled_mat_standardised)) {
      temp = temp + resampled_mat_standardised[, j] * eigenvectors[j, i]
    }
    return(temp)
  })
  pc = do.call(cbind, pc)
  colnames(pc) = paste0("PC", 1:npcs)
  # make pca object
  eigenvectors <- eigenvectors[,1:npcs]
  colnames(eigenvectors) = paste0("PC", 1:npcs)
  #rownames(eigenvectors) = colnames(mass_matrix)
  eigenvalues <- eigenvalues[1:npcs]

  pca = list(
    sdev = sqrt(eigenvalues),
    rotation = eigenvectors,
    center = Matrix::colSums(mass_matrix) / nrow(mass_matrix),
    scale = FALSE,
    x = pc
  )
  pca = list_to_pprcomp(pca)

  verbose_message(message_text = "PCA finished!", verbose = verbose)

  rm(mass_matrix)
  gc()
  eigenvalues = pca$sdev ^ 2
  # Step 5: Compute Principal Components
  # Choose number of principal components, k
  # if not input, use scree test to help find retained components

  if (show_variance_plot) {
    if (!is.null(variance_explained_threshold)) {
      tryCatch({
        cumulative_variance = cumsum(eigenvalues) / sum(eigenvalues)
        threshold = variance_explained_threshold  # Example threshold



        par(mfrow = c(1, 1))
        par(mar = c(2, 2, 1, 1))
        # Plot cumulative proportion of variance explained
        plot(
          cumulative_variance,
          type = 'b',
          main = "Cumulative Variance Explained",
          xlab = "Number of Principal Components",
          ylab = "Cumulative Proportion of Variance Explained"
        )

        # Add a horizontal line at the desired threshold

        abline(h = threshold,
               col = "red",
               lty = 2)


        # Find the number of principal components to retain based on the threshold
        retained =  which(cumulative_variance >= threshold)[1] - 1
      },
      error = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      },
      warning = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      })

    } else{
      # if threshold not inputted, use Kaiser's criterion
      verbose_message(message_text = "Both variance_explained_threshold and npcs not inputted, use Kaiser's criterion for determination", verbose = verbose)


      plot(
        eigenvalues,
        type = 'b',
        main = "Scree Plot",
        xlab = "Principal Component",
        ylab = "Eigenvalue"
      )

      # Add a horizontal line at 1 (Kaiser's criterion)
      abline(h = 1,
             col = "red",
             lty = 2)

      # Add a vertical line at the elbow point
      elbow_point <- which(diff(eigenvalues) < 0)[1]
      abline(v = elbow_point,
             col = "blue",
             lty = 2)
      retained = length(which(eigenvalues >= 1))
    }
  }


  SpaMTP_pca <- pca

  rownames(SpaMTP_pca$rotation) <- SeuratObject::Features(data[[temp_assay]])
  rownames(SpaMTP_pca$x) <- rownames(data@meta.data)

  SpaMTP_pcas <- SeuratObject::CreateDimReducObject(embeddings = SpaMTP_pca$x, loadings = SpaMTP_pca$rotation, assay = assay, key = "pca_", stdev = SpaMTP_pca$sdev)

  SpaMTP[[reduction.name]] <- SpaMTP_pcas

  return(SpaMTP)
}






#' Perform Dimensionality Reduction using Graph-Regularised PCA on Spatial Data
#'
#' Computes a graph-regularised PCA using spatial coordinates and scaled expression data. A k-nearest neighbour (k-NN) graph is computed using spatial locations and used to regularise the PCA decomposition via a graph Laplacian. The result are stored in the `@reductions` section of the returned SpaMTP Seurat object.
#'
#' Note: This method has been adapted from [GraphPCA](https://doi.org/10.1186/s13059-024-03429-x) python package.
#'
#' @param data A SpaMTP Seurat object containing spatial data (feature data and spatial coordinates).
#' @param n_components Integer specifying the number of principal components to compute (default = 50).
#' @param assay Character string defining the name of the assay to use data from (default = "Spatial").
#' @param slot Character string defining the name of the slot to extract scaled data from (default = "scale.data").
#' @param image Character string matching the name of the image to use for extracting spatial coordinates (default = "slice1").
#' @param platform Character string matching either `"Visium"` or `"ST"` to determine how the k-NN graph is constructed. If "Visium" k-nns will handle the hexagon spot arrangement, including setting `n_neighbors` = 6, else "ST" assignment will set `n_neighbors` = 4 unless a value is specifically provided (default = "Visium").
#' @param lambda Numeric value defining the regularisation parameter that controls the influence of the graph Laplacian (default = 0.5).
#' @param n_neighbors Integer value specifying the number of spatial neighbours to use. If `NULL`, will default of 6 for "Visium" data and 4 for "ST" platforms (default = NULL).
#' @param include_self Boolean logical value indicating whether to include self-connections in the graph (default = FALSE).
#' @param alg Character string specifying the algorithm to use for nearest neighbour search (passed to `FNN::get.knn()`) (default ="kd_tree").
#' @param fast Boolean logical value stating whether to use fast approximate eigendecomposition via `RSpectra::eigs_sym()`. For large datasets this is recommended (default =TRUE).
#' @param graph_name Character string of the name to use for storing the computed spatial graph in `@graphs`(default ="SpatialKNN").
#' @param reduction_name Character string stating the name of the dimensionality reduction stored in `@reductions`(default ="SpatialPCA").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A SpaMTP Seurat object with a new graph stored in `@graphs` and spatially-aware PCA reduction values stored in `@reductions`.
#' @export
#'
#' @importFrom Matrix sparseMatrix
#'
#' @examples
#' # spamtp_obj <- RunSpatialGraphPCA(spamtp_obj, platform = "Visium")
RunSpatialGraphPCA <- function(data, n_components=50, assay = "Spatial", slot = "scale.data", image = "slice1", platform="Visium", lambda=0.5, n_neighbors=NULL, include_self = FALSE, alg = "kd_tree", fast = TRUE, graph_name = "SpatialKNN", reduction_name = "SpatialPCA", verbose = TRUE){

  if(!platform %in% c("Visium", "ST")){
    stop("Incorrect value for plantform! platform must be assigned either 'Visium' or 'ST'... If data is not Visium 'spot'-based then set platform to 'ST'")
  }

  if (!is.numeric(lambda)){
    stop("Incorrect value for lamda! Must be numeric value ... ")
  }

  Expr = t(data[[assay]][slot])


  if(is.null(n_neighbors)){

    if (platform == "Visium"){
      n_neighbors <- 6
      verbose_message(message_text = "No nearest neighbour value provided. Based on 'Visium' platform type setting n-neighbours = 6 ..." , verbose = verbose)
    }else{
      verbose_message(message_text = "No nearest neighbour value provided. Based on 'ST' platform type setting n-neighbours = 4 ..." , verbose = verbose)
      n_neighbors <- 4
    }
  } else {
    n_neighbors <- is.integer(n_neighbors)
    verbose_message(message_text = paste0("Using a nearest neighbour value = ",n_neighbors,"! NOTE: the recomended values are -> 'Visium' data: n_neighbors = 6  / 'ST' data: n_neighbors = 4 ...." ), verbose = verbose)
  }

  location = SeuratObject::GetTissueCoordinates(data, image = image)
  location = location[c("x", "y")]


  graph <- kneighbors_graph(location, n_neighbors = n_neighbors, platform = platform, include_self = include_self, alg = alg)
  graph <- 0.5 * (graph + t(graph))
  graphL <- igraph::graph_from_adjacency_matrix(graph, mode = "undirected", weighted = TRUE)
  graphL <- igraph::laplacian_matrix(graphL, normalized = FALSE)
  graphL <- Matrix::as.matrix(graphL)


  # Create identity matrix and add lambda * graphL
  n <- nrow(Expr)

  G <- sparseMatrix(i = 1:n, j = 1:n, x = rep(1, n)) + (lambda * graphL)

  X <- solve(G, Expr)


  rownames(graph) <- colnames(graph) <- colnames(data[[assay]][slot])
  data@graphs[[graph_name]] <- graph

  rm(G)
  rm(graph)
  rm(graphL)

  # 3. Eigendecomposition of C
  if (fast){
    eig <- RSpectra::eigs_sym(crossprod(Expr, X), k = n_components)
  } else {
    eig <- eigen(t(Expr) %*% X, symmetric = TRUE)
  }

  W <- eig$vectors[, seq_len(n_components), drop = FALSE]
  Z <- X %*% W

  rownames(Z) <- colnames(data[[assay]][slot])      # cells
  rownames(W) <- rownames(data[[assay]][slot])      # genes (or features)
  colnames(Z) <- paste0("PC_", 1:ncol(Z))
  colnames(W) <- paste0("PC_", 1:ncol(W))

  # Create the PCA reduction
  pca_reduction <- SeuratObject::CreateDimReducObject(
    embeddings = as.matrix(Z),
    loadings = W,
    assay = assay,
    key = "PC_"
  )

  # Store in the Seurat object
  data[[reduction_name]] <- pca_reduction

  return(data)

}


#' Construct a k-Nearest Neighbour Graph from Spatial Coordinates.
#'
#' Builds a sparse adjacency matrix representing the k-nearest neighbour relationships between spots or cells, using Euclidean distance between spatial coordinates.
#'
#' @param location A data frame or matrix with columns `x` and `y` representing spatial coordinates, and rownames match barcode/cell names.
#' @param n_neighbors Integer value defining the number of neighbours to consider.
#' @param platform Character string being either `"Visium"` or `"ST"`, specifying the platform of the dataset being analysed. Determines how distance thresholds and adjacency are calculated whereby "Visium" data is handled slightly differently to compensate for the hexagon shape of the spots.
#' @param include_self Boolean logical value stating whether to include self-connections (default = FALSE).
#' @param alg Character string specifying the algorithm to use for nearest neighbour search. Passed to `FNN::get.knn()` (default = "kd_tree").
#'
#' @return A sparse adjacency matrix (`dgCMatrix`) representing the k-NN graph.
#' @export
#'
#' @examples
#' ### HELPER FUNCTION
kneighbors_graph <- function(location, n_neighbors, platform, include_self = FALSE, alg = "kd_tree") {
  # Get the k-nearest neighbors using kd_tree (most similar to scikit-learn's default)
  knn_result <- FNN::get.knn(location, k = n_neighbors, algorithm = alg)

  if(platform == "Visium"){

    # Step 1: Flatten the distances and indices
    N <- nrow(knn_result$nn.index)

    # Reshape distances and column indices
    dists <- as.vector(t(knn_result$nn.dist))
    col_indices <- as.vector(t(knn_result$nn.index))

    # Create row indices (equivalent to np.repeat())
    row_indices <- rep(1:N, each = n_neighbors)

    # Apply neighbor correction if needed

    dist_cutoff <- median(dists) * 1.3  # Small amount of sway
    mask <- dists < dist_cutoff

    # Filter using the mask
    row_indices <- row_indices[mask]
    col_indices <- col_indices[mask]
    dists <- dists[mask]

    adjacency <- sparseMatrix(
      i = row_indices,
      j = col_indices,
      x = rep(1, length(row_indices)),
      dims = c(N, N),
      repr = "C"  # Ensures a "dgCMatrix" (like CSR in Python)
    )
  } else {

    # Create empty adjacency matrix
    n <- nrow(location)
    adjacency <- matrix(0, nrow = n, ncol = n)

    # Fill adjacency matrix with connections
    for (i in 1:n) {
      adjacency[i, knn_result$nn.index[i,]] <- 1
    }
    adjacency <- Matrix(adjacency, sparse = TRUE)


  }

  # Include self-loops if requested (in this case, no)
  if (include_self) {
    diag(adjacency) <- 1
  }


  return(adjacency)
}



#' Perform K-means clustering on a specified reduction
#'
#' This function runs K-means clustering on a specified reduction in a SpaMTP Seurat object and adds the cluster assignments to the object metadata.
#'
#' @param data A SpaMTP Seurat object containing the results from `RunSpatialGraphPCA()`.
#' @param reduction Character string stating the name of the reduction slot to use (default = "SpatialPCA").
#' @param cluster.name Character string of the name of the metadata column to store the cluster labels (default = "spatial_clusters").
#' @param centers Integer defining the number of clusters to form (default = 8).
#' @param iter.max Integer defining the maximum number of iterations allowed (default = 10).
#' @param nstart Integer stating the number of random sets to choose (default = 1).
#' @param algorithm Character string defining the K-means algorithm to use. One of `"Hartigan-Wong"`, `"Lloyd"`, `"Forgy"`, or `"MacQueen"` (default = "Hartigan-Wong").
#' @param trace Logical boolean indicating whether to produce tracing information on the progress of the algorithm (default = FALSE).
#' @param seed Integer of the random seed to use for reproducibility (default = 888).
#'
#' @return A SpaMTP Seurat object with a new metadata column containing the K-means cluster assignments.
#' @export
#'
#' @examples
#' # seurat_object <- GetKmeanClusters(spamtp_obj, reduction = "SpatialPCA", centers = 8, cluster.name = "test_clusters")
#' # SpatialDimPlot(spamtp_obj, group.by = "test_clusters")
GetKmeanClusters <- function(data, reduction = "SpatialPCA", cluster.name = "spatial_clusters", clusters = 8, iter.max = 10, nstart = 1, algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), trace = FALSE, seed = 888){

  if (!reduction %in% names(data@reductions)){
    stop("Reduction not present in SpaMTP Seruat Object! ", "'",reduction, "' was not found, please run names(data@reductions) to check possible reductions to use ....")
  }
  set.seed(seed)
  res <- stats::kmeans(data[[reduction]]@cell.embeddings, centers = centers, iter.max = iter.max, nstart = nstart, algorithm = algorithm,  trace = trace)

  data[[cluster.name]] <- res$cluster
  return(data)
}

