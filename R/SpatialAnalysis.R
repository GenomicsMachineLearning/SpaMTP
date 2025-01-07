#' Find top features and metabolites that are strongly correlated with a given feature
#'
#' @param data SpaMTP Seurat class object containing both Spatial Transcriptomic and Metabolic data assays.
#' @param mz Numeric string specifying the m/z to find correlated features for. One of `mz`, `gene` or `ident` must be provided, alternatives must be `NULL` (default = NULL).
#' @param gene Character string specifying the gene to find correlated features for. One of `mz`, `gene` or `ident` must be provided, alternatives must be `NULL` (default = NULL).
#' @param ident Character string defining the ident column in the data object's `@meta.data` slot to find correlated features for. One of `mz`, `gene` or `ident` must be provided, alternatives must be `NULL` (default = NULL).
#' @param SM.assay Character string specifying the name of the assay containing the spatial metabolomics (SM) data (default = "SPM").
#' @param ST.assay Character string specifying the name of the assay containing the spatial transcriptomics (ST) data. If NULL then only metabolites will be used (Default = NULL).
#' @param SM.slot Character string specifying the slot of the SM assay to use (default = "counts").
#' @param ST.slot Character string specifying the slot of the ST assay to use (default = "counts").
#' @param nfeatures Integer specifying the number of top correlated features to return (default = 10).
#'
#' @return A data frame containing the top correlated features with columns for the feature names and their correlation values.
#'
#' @export
#'
#' @examples
#' # result <- FindCorrelatedFeatures(data = SpaMTP, gene = "GeneX", nfeatures = 5)
FindCorrelatedFeatures <- function(data, mz = NULL, gene = NULL, ident = NULL, SM.assay = "SPM", ST.assay = NULL, SM.slot = "counts", ST.slot = "counts", nfeatures = 10){

  data_list <- list()
  met_counts <- data[[SM.assay]][SM.slot]

  if (!is.null(ST.assay)){
    tran_counts <- data[[ST.assay]][ST.slot]

    gene_mappings <- data.frame(gene = rownames(tran_counts))

    rownames(tran_counts) <- unlist(lapply(1:length(rownames(tran_counts)), function (x) {
      paste0("mz-",(round(as.numeric(gsub("mz-", "", rownames(met_counts)[length(rownames(met_counts))],))) + 100), x)
    }))

    gene_mappings$mz <- rownames(tran_counts)
    gene_mappings$raw_mz <- gsub("mz-", "", gene_mappings$mz)

    data[["tmp"]] <- SeuratObject::CreateAssay5Object(counts = rbind(met_counts, tran_counts))
    data[["tmp"]][SM.slot] <- data[["tmp"]]["counts"]
    SM.assay <- "tmp"
  } else{
    if (!is.null(ST.slot)){
      warning("`ST.slot` is not set to NULL, but `ST.assay` is ... Gene's will not be included in the analysis. If genes are to be included, please set `ST.assay` != NULL -> using the appropriate assay")
    }
  }

  if (is.null(mz) & ! is.null(gene) & is.null(ident)){  # for gene mapping
    idx <- which(gene_mappings$gene == gene)
    mz <- as.numeric(gene_mappings$raw_mz[idx])

  } else if (!is.null(mz) & is.null(gene) & is.null(ident)){ # for mz mapping
    mz <- mz
  } else if (!is.null(ident) & is.null(gene) & is.null(mz)){ # for ident mapping

    for (i in unique(data@meta.data[[ident]])){
      data_list[[i]] <- unlist(lapply(data@meta.data[[ident]], function(x) { ifelse(x == i, TRUE, FALSE)}))
    }

  } else {
    stop("Invalid input for 'mz = ', 'ident' = and 'gene = '... Only one inupt can be provided. Either mz, ident or gene, alternative must be set to NULL! Please check documentation ...")
  }

  data_cardinal <- ConvertSeuratToCardinal(data = data, assay = SM.assay, slot = SM.slot)

  if (!is.null(ident)){
    for (i in names(data_list)){
      data_list[[i]] <- Cardinal::colocalized(object = data_cardinal, ref = data_list[[i]], n = length(rownames(data[[SM.assay]][SM.slot])))
    }
  } else {
    data_list[["1"]] <- suppressWarnings(Cardinal::colocalized(data_cardinal, mz=mz, n = length(rownames(data[[SM.assay]][SM.slot]))))
  }

  if (!is.null(ST.assay)){
    for (i in names(data_list)){

      result <- data_list[[i]][order(data_list[[i]]$mz), ]
      result$features <- c(rownames(met_counts),gene_mappings$gene)
      result$modality <- c(rep("metabolite", length(rownames(met_counts))), c(rep("gene", length(gene_mappings$gene))))
      result <- result[c("features", colnames(result)[!colnames(result) %in% c("mz", "features")])]
      data_list[[i]] <- result
    }

  }

  for (i in names(data_list)){

    result <- data_list[[i]]
    result <- result[order(-abs(result$correlation)), ]
    result$ident <- i
    result$rank <- c(1:length(result$ident))
    data_list[[i]] <- result

  }

  results <- data.frame(do.call(rbind, data_list))

  if (!is.null(nfeatures)){
    results <- data.frame(results %>%
                            dplyr::group_by(ident) %>%
                            dplyr::slice_head(n = nfeatures) %>%
                            dplyr::ungroup())

  }


  if ( length(data_list) == 1){
    results$ident <- NULL
  }


  return(results)
}




#### Code is manipulated from Seurat ###########################################################################################

#' Compute the row variances for each m/z value
#'
#' @param x Matrix containing count values used for correlation
#'
#' @return Vector contating the row variances for each m/z value
#'
#' @examples
#' #HELPER FUNCTION
RowVar <- function(x) {
  .Call('_Seurat_RowVar', PACKAGE = 'Seurat', x)
}


#' Finds metabolites that display strong spatial patterns using MoransI
#'
#' @param object SpaMTP Seurat class object contating the intensity values for each m/z
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "SPM").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param image Character string defining the image to extract the tissue coordinates from (defualt = "slice1").
#' @param nfeatures Numeric values defining the top number of features to mark as the top spatially variable (default = 2000).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return Returns SpaMTP object containing the MoransI pvalue and rank stored in the assay feature meta.data
#' @export
#'
#' @examples
#' # SpaMTP.obj <- FindSpatiallyVariableMetabolites(SpaMTP)
FindSpatiallyVariableMetabolites <- function(object, assay = "SPM", slot = "counts",image = "slice1", nfeatures = 2000, verbose = TRUE){

  features <- rownames(x = object[[assay]])
  spatial.location <- GetTissueCoordinates(object = object[[image]])
  data <- GetAssayData(object = object, slot = slot)
  data <- as.matrix(x = data[features, ])
  data <- data[RowVar(x = data) > 0, ]
  svf.info <- RunMoransI(data = data, pos = spatial.location, verbose = verbose)
  colnames(x = svf.info) <- paste0("MoransI_", colnames(x = svf.info))
  var.name <- paste0("moransi", ".spatially.variable")
  var.name.rank <- paste0(var.name, ".rank")
  svf.info[[var.name]] <- FALSE
  svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[, 1])), , drop = FALSE]
  svf.info[[var.name]][1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
  svf.info[[var.name.rank]] <- 1:nrow(x = svf.info)
  object[[assay]][[names(x = svf.info)]] <- svf.info

  return(object)
}

#' Gets the top n number of spatially variable features
#'
#' @param object SpaMTP Seurat class object containing the intensity values for each m/z
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "SPM").
#' @param n Numeric value defining the top number of metabolites to return (default = 10).
#'
#' @return Vector containing m/z feature names corresponding to the top n number of spatially variable metabolites
#' @export
#'
#' @examples
#' # features <- GetSpatiallyVariableMetabolites(SpaMTP, n = 6)
GetSpatiallyVariableMetabolites <- function(object, assay = "SPM", n = 10){

  return(rownames(object[[assay]][["moransi.spatially.variable.rank"]]%>%arrange(moransi.spatially.variable.rank))[1:n])
}



#' Multi-Omic integration of Spatial Metabolomics and Transcriptomics data using Seurat's Weighted Nearest Neighbours function
#'
#' @param multiomic.data SpaMTP dataset contain Spatial Transcriptomics and Metabolomic datasets in two different assays
#' @param weight.list List containing the relative weightings for each modality, matching the reduction order. If NULL, weights will be automatically calculated else, two values must add to 1 (default = NULL).
#' @param reduction.list List containing character strings defining the reduction to use for each modality, in the order matching weight.list if applicable (default = list("spt.pca", "spm.pca")).
#' @param dims.list List containing the numeric range of principle component dimension to include for each modality (default = list(1:30,1:30)).
#' @param return.intermediate Boolean value indicating whether to store intermediate results in misc slot of SpaMTP Seurat class object (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param ... Additional arguments that can be parsed through Seurat's FindMultModalNeighbors function. For possible inputs please visit: https://www.rdocumentation.org/packages/Seurat/versions/5.0.3/topics/FindMultiModalNeighbors.
#'
#' @return SpaMTP Seurat class object containing a weighted nearest neighbours graph which integrates Metabolic and Transcriptomic modalities. This graph can be used for clustering.
#' @export
#'
#' @examples
#' # SpaMTP.obj <- MultiOmicIntegration(SpaMTP.obj, weight.list = list(0.5, 0.5), reduction.list =  list("spt.pca", "spm.pca"), dims.list = list(1:30, 1:30))
MultiOmicIntegration <- function (multiomic.data, weight.list = NULL, reduction.list =  list("spt.pca", "spm.pca"), dims.list = list(1:30, 1:30), return.intermediate = FALSE, verbose = FALSE, ...){

  if (is.null(weight.list)){
    mm.integration <- Seurat::FindMultiModalNeighbors(
      multiomic.data, reduction.list = reduction.list,
      dims.list = dims.list, return.intermediate = return.intermediate,verbose = verbose, ...)
  } else {
    mm.integration <- Seurat::FindMultiModalNeighbors(
      multiomic.data, reduction.list = reduction.list,
      dims.list = dims.list, return.intermediate = TRUE, verbose = verbose, ...)

    x <- rep(weight.list[[1]], length(names(mm.integration@misc$modality.weight@modality.weight.list[[reduction.list[[1]]]]))) ## Setting the SPM weights
    names(x) <- names(mm.integration@misc$modality.weight@modality.weight.list[[reduction.list[[1]]]])
    mm.integration@misc$modality.weight@modality.weight.list[[reduction.list[[1]]]] <- x

    x <- rep(weight.list[[2]], length(names(mm.integration@misc$modality.weight@modality.weight.list[[reduction.list[[2]]]]))) ## Setting the SPM weights
    names(x) <- names(mm.integration@misc$modality.weight@modality.weight.list[[reduction.list[[2]]]])
    mm.integration@misc$modality.weight@modality.weight.list[[reduction.list[[2]]]] <- x

    mm.integration <- Seurat::FindMultiModalNeighbors(
      multiomic.data, reduction.list = reduction.list,
      dims.list = dims.list, return.intermediate = return.intermediate, modality.weight = mm.integration@misc$modality.weight, verbose = verbose, ...)
  }

  return(mm.integration)
}



#' Helper function that generated PCA analysis results for a SpaMTP Seurat Object
#'
#' @param SpaMTP SpaMTP Seurat class object that contains spatial metabolic information.
#' @param npcs is an integer value to indicated preferred number of PCs to retain.
#' @param variance_explained_threshold Numeric value defining the explained variance threshold.
#' @param resampling_factor is a numerical value > 0, indicate how you want to resample the size of original matrix.
#' @param assay Character string defining the SpaMTP assay to extract intensity values from.
#' @param slot Character string defining the assay slot containing the intensity values.
#' @param show_variance_plot Boolean indicating weather to display the variance plot output by this analysis.
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed.
#'
#'
#' @return PCA object contating embeddings and projections
#'
#' @import dplyr
#'
#' @examples
#' # HELPER FUNCTION
getPCA <- function(SpaMTP,
                   npcs,
                   variance_explained_threshold,
                   assay,
                   slot,
                   show_variance_plot,
                   bin_resolution = NULL,
                   resolution_units = NULL,
                   bin_method = NULL,
                   verbose) {

}



#' Generates PCA analysis results for a SpaMTP Seurat Object
#'
#' @param SpaMTP SpaMTP Seurat class object that contains spatial metabolic information.
#' @param npcs is an integer value to indicated preferred number of PCs to retain (default = 30).
#' @param variance_explained_threshold Numeric value defining the explained variance threshold (default = 0.9).
#' @param resampling_factor is a numerical value > 0, indicate how you want to resample the size of original matrix (default = 1).
#' @param byrow is a boolean to indicates whether each column of the matrix is built byrow or bycol (default = FALSE).
#' @param assay Character string defining the SpaMTP assay to extract intensity values from (default = "SPM").
#' @param slot Character string defining the assay slot containing the intensity values (default = "counts").
#' @param show_variance_plot Boolean indicating weather to display the variance plot output by this analysis (default = FALSE).
#' @param reduction.name Character string indicating the name associated with the PCA results stored in the output SpaMTP Seurat object (default = "pca").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#'
#' @return SpaMTP object with pca results stored in the
#' @export
#'
#' @examples
#' # HELPER FUNCTION
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

  return(data)
}


#' Creates a pprcomp object based on an input list
#'
#' @param lst List containing PCA results
#'
#' @return A pprcomp object contating results from PCA analysis
#'
#' @examples
#' #HELPER FUNCTION
list_to_pprcomp <- function(lst) {
  # Create an empty object with class pprcomp
  obj <- structure(list(), class = "prcomp")
  # Assign components from the list to the object
  obj$sdev <- lst$sdev
  obj$rotation <- lst$rotation
  obj$center <- lst$center
  obj$scale <- lst$scale
  obj$x <- lst$x
  # Add other components as needed

  # Return the constructed pprcomp object
  return(obj)
}


#' Finds the index values of the m/z values with their respective GSEA result
#'
#' @param lst List containing relative mz analytes and pathways
#' @param value Value returned based on the GSEA results
#'
#' @return returns a vector of indices that match the relative GSEA results to the m/z list
#'
#' @examples
#' #HELPER FUNCTION
find_index <- function(lst, value) {
  indices <- which(sapply(lst, function(x) value %in% x))
  if (length(indices) == 0) {
    return(NULL)  # If value not found, return NULL
  } else {
    return(indices)
  }
}


#####################################################
### Bellow helper functions are sourced from https://stackoverflow.com/questions/11123152/function-for-resizing-matrices-in-r by Vyga.


#' Helper function to rescale a sampled matrix
#'
#' @param x Vector defining new matrix coordinates
#' @param newrange Vector defining range of old coordinates
#'
#' @return Vector containing rescaled coordinates
#'
#' @examples
#' #HELPER FUNCTION
rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

#' Helper function to resize a matrix back to its original layout after sampling
#'
#' @param mat matrix defining the sampled object
#' @param ndim number of dimentions to resize the sampled matrix to (default = dim(mat)).
#'
#' @return Returns a resized sampled matrix to match the dimentions of the original
#'
#' @examples
#' #HELPER FUNCTION
ResizeMat <- function(mat, ndim=dim(mat)){
  #if(!require(fields)) stop("`fields` required.")

  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)

  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)

  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))

  # interpolation
  ans[ncord] <- fields::interp.surface(obj, loc)

  ans
}

#######################################################################################################################

