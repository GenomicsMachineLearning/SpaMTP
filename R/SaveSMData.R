
#### SpaMTP Saving Data Objects ########################################################################################################################################################################################

#' Saves SpaMTP Object
#'
#' This function saves a SpaMTP Seurat Object into a standard single-cell/spatial file format.
#' This includes a filtered_feature_bc_matrix folder containing files storing the features, barcode/pixels and intensity matrix.
#' Metadata and sapatial files (such as scale factors and hires/lowres images) are also stored.
#'
#' @param data A Spatial Metabolomic SpaMTP Seurat Object being saved.
#' @param outdir Character string of the directory to save the mtx.mtx, barcode.tsv, features.tsv, barcode_metadata.csv and feature_metadata.csv in.
#' @param assay Character string defining the Seurat assay that contains the m/z count data (default = "Spatial").
#' @param slot Character string defining the Seurat assay slot that contains the m/z values directly (default = "counts").
#' @param image Character string defining the image stored within the SpaMTP Seurat object to save. If `NULL` no image will be saved (default = NULL).
#' @param annotations Boolean values defining if the Seurat Object contains annotations to be saved (default = FALSE).
#' @param verbose Boolean indicating whether to show informative processing messages. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' ### Details
#' * This can be used for saving data for transfer to python. Can be read in as Anndata using scanpy.read_10x_mtx().
#' * For saving in R saveRDS() is recommended.
#'
#' @export
#'
#' @examples
#' # saveSpaMTPData(SeuratObject, "../output", annotations = TRUE)
SaveSpaMTPData <- function(data, outdir, assay = "Spatial", slot = "counts", image = NULL, annotations = FALSE, verbose = TRUE){

  if (!dir.exists(outdir)) {
    verbose_message(message_text = paste0("Generating new directory to store output here: ", outdir), verbose = verbose)
    dir.create(outdir)
  } else {
    verbose_message(message_text = paste0("Directory already exists, storing output here: ", outdir), verbose = verbose)
  }

  verbose_message(message_text = paste0("Writing ", slot," slot to matrix.mtx, barcode.tsv, genes.tsv"), verbose = verbose)
  DropletUtils::write10xCounts(data[[assay]][slot], path = paste0(outdir,"/filtered_feature_bc_matrix/"), overwrite = TRUE)

  verbose_message(message_text = "Writing @metadata slot to metadata.csv", verbose = verbose)
  data.table::fwrite(data@meta.data, paste0(outdir,"/barcode_metadata.csv"))

  if (!is.null(image)){

    verbose_message(message_text ="Generating 'spatial' directory ... ", verbose = verbose)
    dir.create(paste0(outdir, "/spatial/"))

    if ("scale.factors" %in% slotNames(data@images[[image]])){
      scale.factors <- list("tissue_hires_scalef" = data@images[[image]]@scale.factors[["hires"]],
                          "tissue_lowres_scalef" = data@images[[image]]@scale.factors[["lowres"]],
                          "fiducial_diameter_fullres" = data@images[[image]]@scale.factors[["fiducial"]],
                          "spot_diameter_fullres" = data@images[[image]]@scale.factors[["spot"]]* 2)

      sfJSON <- jsonlite::toJSON(
        rapply(scale.factors, function(x) if (length(x) == 1L) jsonlite::unbox(x) else x,
               how = "replace"))

      verbose_message(message_text ="Writing scalefactors_json.json file ...", verbose = verbose)
      write(sfJSON, file = paste0(outdir, "/spatial/scalefactors_json.json"))
    }

    if ("image" %in% slotNames(data@images[[image]])){
      verbose_message(message_text ="Writing image to `spatial/tissue_lowres_image.png` ...", verbose = verbose)
      png::writePNG(data@images[[image]]@image, paste0(outdir, "/spatial/tissue_lowres_image.png"))
    }

    coords <- GetTissueCoordinates(data, image = image)
    coords$in_tissue <- 1

    x_coords <- sort(unique(coords[,"x"]))
    names(x_coords) <- 1:length(x_coords)
    coords$arrayrow <- names(x_coords)[match(coords$x, x_coords)]

    y_coords <- sort(unique(coords[,"y"]))
    names(y_coords) <- 1:length(y_coords)
    coords$arraycol <- names(y_coords)[match(coords$y, y_coords)]

    coords <- coords[c("cell", "in_tissue", "arrayrow", "arraycol", "x","y")]
    colnames(coords) <- NULL

    verbose_message(message_text ="Writing tissue coordinate file` ...", verbose = verbose)
    data.table::fwrite(coords, paste0(outdir,"/spatial/tissue_positions_list.csv"))

  }
  if (annotations){
    verbose_message(message_text = "Writing feature metadata annotations to feature_metadata.csv", verbose = verbose)
    data.table::fwrite(data[[assay]]@meta.data, paste0(outdir,"/feature_metadata.csv"))
  }

}

########################################################################################################################################################################################################################
