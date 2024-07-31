#' Helper function for suppressing function progress messages
#'
#' @param message_text Character string containing the message being shown
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the messsage will be suppressed (default = TRUE).
#'
verbose_message <- function(message_text, verbose) {
  if (verbose) {
    message(message_text)
  }
}



#Functiona modified from: https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat_v2.R

#'@importFrom magrittr %>% %<>%
NULL

#' Intermediate solution to \code{subset()}:
#' subset FOVs/centroids if selected cells are NOT found in each FOV
#' NOTE: some code parts and args are taken from SeuratObject

#' Function params/args:
#' @param object An S4 object or A \code{FOV} object
#' @param subset Logical expression indicating features/variables to keep
#' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
#' @param idents A vector of identity classes to keep
#' @param features A vector of feature names or indices to keep
#' @param Update.slots If to update slots of an object
#' @param Update.object If to update final object, default to TRUE.
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param ... Arguments passed to \code{subset()} and other methods
#'
#' @return A subset Seurat object
#' @export
#'
#' @examples
#' # sub <- subset_obt(seurat.obj, idents = "Sample1")
subset_SPM <- function(
    object = NULL,
    subset = NULL,
    cells = NULL,
    idents = NULL,
    features = NULL,
    Update.slots = TRUE,
    Update.object = TRUE,
    verbose = TRUE,
    ...)
{

  if (Update.slots) {
    verbose_message(message_text = "Updating object slots..", verbose = verbose)
    object %<>% UpdateSlots()
  }

  verbose_message(message_text = "Cloing object..", verbose = verbose)
  obj_subset <- object

  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) {
    cells <- Cells(obj_subset)[cells]
  }

  if (!missing(subset) || !is.null(idents)) {
    verbose_message(message_text = "Extracting cells matched to `subset` and/or `idents`", verbose = verbose)
  }

  if (class(obj_subset) == "FOV") {
    verbose_message(message_text = "object class is `FOV` ", verbose = verbose)
    cells <- Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    cells <-
      WhichCells(object = obj_subset,
                 cells = cells,
                 idents = idents,
                 expression = subset,
                 return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      WhichCells(object = obj_subset,
                 cells = cells,
                 idents = idents,
                 return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- Cells(obj_subset)
  }

  # added support for object class `FOV`
  if (class(obj_subset) == "FOV") {
    verbose_message(message_text = "Matching cells for object class `FOV`..", verbose = verbose)
    cells_check <- any(obj_subset %>% Cells %in% cells)
  } else {
    # check if cells are present in all FOV
    verbose_message(message_text = "Matching cells in FOVs..", verbose = verbose)

    cells_check <-
      lapply(Images(obj_subset) %>% seq,
             function(i) {
               any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)
             }) %>% unlist
  }

  if (all(cells_check)) {
    verbose_message(message_text = paste0("Cell subsets are found in all FOVs!", "\n",
                    "Subsetting object.."), verbose = verbose)


    obj_subset %<>% base::subset(cells = cells,
                                 idents = idents,
                                 features = features,
                                 ...)
    # subset FOVs
    verbose_message(message_text = "Subsetting FOVs..", verbose = verbose)

    fovs <-
      lapply(Images(obj_subset) %>% seq, function(i) {
        base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                     cells = cells,
                     idents = idents,
                     features = features,
                     ...)
      })
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }

  } else {
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <-
      lapply(Images(obj_subset) %>% seq, function(i) {
        if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {

          verbose_message(message_text = paste0("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i]), verbose = verbose)
          verbose_message(message_text = "Subsetting Centroids..", verbose = verbose)

          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells,
                       idents = idents,
                       features = features,
                       ...)
        }
      })
    # remove FOVs with no matching cells

    verbose_message(message_text = paste0("Removing FOVs where cells are NOT found: ", "\n",
                    paste0(Images(object)[which(!cells_check == TRUE)], "\n")), verbose = verbose)

    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }

    # subset final object
    verbose_message(message_text = "..subset final object", verbose = verbose)

    obj_subset %<>%
      base::subset(cells = cells,
                   idents = idents,
                   features = features,
                   ...)
  }

  if (Update.object && !class(obj_subset) == "FOV") {
    verbose_message(message_text = "Updating object..", verbose = verbose)

    obj_subset %<>% UpdateSeuratObject() }

  verbose_message(message_text = "Object is ready!", verbose = verbose)
  return(obj_subset)

}



########################################################################
