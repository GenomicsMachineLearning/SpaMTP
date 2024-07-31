library(Seurat)
library(DropletUtils)
library(data.table)

#### SpaMTP Saving Data Objects ########################################################################################################################################################################################

#' Saves Seurat Spatial Metabolomic object to files.
#'    - This can be used for saving data for transfer to python Anndata using scanpy.read_10x_mtx()
#'    - For saving in R saveRDS() is recommended
#'
#' @param data A Seurat Spatial Metabolomic Object being saved.
#' @param outdir Character string of the directory to save the mtx.mtx, barcode.tsv, features.tsv, barcode_metadata.csv and feature_metadata.csv in.
#' @param assay Character string defining the Seurat assay that contains the m/z count data (default = "Spatial").
#' @param slot Character string defining the Seurat assay slot that contains the m/z values directly (default = "counts").
#' @param annotations Boolean values defining if the Seurat Object contains annotations to be saved (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @export
#'
#' @examples
#' # saveSeuratData(SeuratObject, "~/Documents/seuratobj_files/", annotations = TRUE)
saveSeuratData <- function(data, outdir, assay = "Spatial", slot = "counts", annotations = FALSE, verbose = TRUE){

  verbose_message(message_text = paste0("Generating new directory to store output here: ", outdir), verbose = verbose)
  verbose_message(message_text = paste0("Writing ", slot," slot to matrix.mtx, barcode.tsv, genes.tsv"), verbose = verbose)
  DropletUtils::write10xCounts(data[[assay]][slot], path = outdir, overwrite = TRUE)

  verbose_message(message_text = "Writing @metadata slot to metadata.csv", verbose = verbose)
  data.table::fwrite(data@meta.data, paste0(outdir,"barcode_metadata.csv"))

  if (annotations){
    verbose_message(message_text = "Writing feature metadata annotations to feature_metadata.csv", verbose = verbose)
    data.table::fwrite(data[[assay]]@meta.data, paste0(outdir,"feature_metadata.csv"))
  }

}

########################################################################################################################################################################################################################
