#' Mult-Omic data integration
#'
#' This function performs multi-omic integration of Spatial Metabolomics and Spatial Transcriptomics data using Seurat's Weighted Nearest Neighbours function.
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




#' Create a singular multiomics assay by merging data from multiple assays.
#'
#' Combines the `scale.data` slots from multiple assays in a SpaMTP Seurat object into a single new assay.
#' Useful for integrating multiple modalities (e.g. transcriptomics, proteomics, metabolomics) that have already been scaled.
#'
#' @param SpaMTP A SpaMTP Seurat object that contains atleast two assays to be merged.
#' @param assays.to.merge A character vector specifying the names of assays whose `scale.data` slots should be merged. At least two assay names must be provided and both must contain the `scale.data` slot, for example: assays.to.merge = c("SPM", "SPT").
#' @param new.assay A character string specifying the name of the new assay to be created (default = "merged").
#' @param return.original Boolean value defining if the returned SpaMTP Seurat object will contain the original individual assays. If the data size is large it is recommended to set to `False` (default = TRUE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A SpaMTP Seurat object containing a new assay with the merged scaled data values.
#' @export
#'
#' @details
#' This function assumes that each specified assay has been processed with `Seurat::ScaleData()` (or an alternative scaling method), and that their `scale.data` slots contain numeric matrices. The merged assay will use the row-bound `scale.data` matrices as the `counts`, `data`, and `scale.data` slots.
#'
#' @examples
#' # merged_obj <- CreateMergedModalityAssay(SpaMTP = spamtp_obj, assays.to.merge = c("SPM", "SPT"),new.assay = "merged")
CreateMergedModalityAssay <- function(SpaMTP, assays.to.merge, new.assay = "merged", return.original = TRUE, verbose = FALSE){

  if(length(assays.to.merge) < 1){
    stop("Incorrect length of assays.to.merge! atleast two assay names must be provided to combine the scale.data slots. Please adjust assays.to.merge accordingly.")
  }

  for (assay in assays.to.merge){
    if(is.null(SpaMTP@assays[[assay]])){
      stop("Assay does not exist! The provided assay name is not present in the SpaMTP Seurat Obejct.")
    }
    if(is.null(SpaMTP@assays[[assay]]["scale.data"])){
      stop("No scale.data slot present in the ", assay, " assay! Please run Seurat::ScaleData() first!")
    }
  }
  scaled_data <- lapply(assays.to.merge, function(x){
    SpaMTP@assays[[x]]$scale.data
  })

  SpaMTP[[new.assay]] <- SeuratObject::CreateAssay5Object(counts = do.call(rbind, scaled_data))

  message("Warning: Restoring feature names to contain '_' ...")

  rownames(SpaMTP[[new.assay]]) <- gsub("-", "_", x = rownames(SpaMTP[[new.assay]]))

  verbose_message(message_text = "NOTE: the matrix containing merged scaled data has been assigned to the `$counts`, `$data` and `$scaled.data` slots. All matricies are the same ...", verbose = verbose)

  SpaMTP[[new.assay]]$data <- SpaMTP[[new.assay]]$counts
  SpaMTP[[new.assay]]$scale.data <- SpaMTP[[new.assay]]$counts

  if (!return.original){
    Seurat::DefaultAssay(SpaMTP) <- new.assay

    for (a in assays.to.merge){
      SpaMTP[[a]] <- NULL
      }
  }


  return(SpaMTP)
}
