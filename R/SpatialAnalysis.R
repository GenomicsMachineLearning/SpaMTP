#' Find top features and metabolites that are strongly correlated with a given feature
#'
#' This function returns a list of features ranked by highest Pearson correlation score.
#' The specified feature to corrleate against can be either a m/z value, gene or an ident (i.e. cluster).
#' For multi-omic data, both a metabolic and transcriptomic assay can be specified to calculate correlation of both metabolites and genes.
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

    if ("cor" %in% colnames(result)){
      names(result)[names(result) == "cor"] <- "correlation"
    }

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


#' Find Spatially Variable Metabolites
#'
#' Finds metabolites that display strong spatial patterns using MoransI.
#' Each m/z value is ranked by MoransI score and the results are stored in the SpaMTP Seurat Object feature metadata.
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

  DefaultAssay(object) <- assay
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

#' Get top spatially variable metabolites
#'
#' This function returns the names of the top n number of spatially variable features.
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



#######################################################################################################################
## Code below are helper functions
#######################################################################################################################



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
  obj$rotation <- lst$rotationf
  obj$center <- lst$center
  obj$scale <- lst$scale
  obj$x <- lst$x
  # Add other components as needed

  # Return the constructed pprcomp object
  return(obj)
}


#### Code is manipulated from Seurat ##################################################################################

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
