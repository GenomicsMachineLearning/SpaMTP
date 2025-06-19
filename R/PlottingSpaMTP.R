#### SpaMTP Seurat Plotting Functions #################################################################################################################################################################################

#' Finds the nearest m/z peak to a given value in a SpaMTP Object
#'
#' @param data Seurat Spatial Metabolomic object containing mz values
#' @param target_mz Numeric value defining the target m/z peak
#' @param assay Character string indicating which Seurat object assay to pull data form. If set to NULL will use current default assay based on Seurat's `DefaultAssay()` function (default = NULL).
#'
#' @returns String of the closest m/z value within the given dataset
#' @export
#'
#' @examples
#' # FindNearestMZ(SeuratObj, target_mz = 400.01)
FindNearestMZ <- function(data, target_mz, assay = NULL){
  if (is.null(assay)){
    assay <- DefaultAssay(data)
  }
  numbers <- as.numeric(gsub("mz-", "", SeuratObject::Features(data, assay = assay)))
  closest_number <- numbers[which.min(abs(numbers - target_mz))]
  return(paste0("mz-",closest_number))
}


#' Sums the intensity values of multiple m/z values into one
#'
#' This function inputs any vector of m/z values and combines their intensity values per pixel.
#' This 'binned' data is stored in the SpaMTP Object's `@meta.data` slot.
#'
#' @param data SpaMTP Seurat class object containing m/z intensities.
#' @param mzs Vector of m/z names to be binned together
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "Spatial").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param bin_name Character string defining the name of the meta.data column that stores the data (default = "Binned_Metabolites").
#'
#' @return Binned intensity value stored in barcode meta.data slot
#' @export
#'
#' @examples
#' # SpaMPT.obj <- BinMetabolites(SpaMPT.obj, mz = c('mz-740.471557617188','mz-784.528564453125','mz-897.603637695312'), bin_name = "Lipids")
BinMetabolites <- function(data, mzs, assay = "Spatial", slot = "data", bin_name = "Binned_Metabolites"){

  binned_counts <- bin.mz(data, mz_list = mzs, assay = assay, slot = slot) # bins m/z masses

  data[[bin_name]]<- binned_counts #adds to metadata

  return(data)
}

#' Bins multiple m/z values into one.
#'
#' This function is a helper function for `BinMetabolites`.
#'
#' @param data Seurat Spatial Metabolomic object containing mz values
#' @param mz_list Vector of m/z names to be binned into one value
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "Spatial").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param stored.in.metadata Boolean value indicating if the mz_list should be searched in the metadata DataFrame. If FALSE searches in Seurat object assay slot, if TRUE searches in metadata slot (default = FALSE).
#'
#' @returns Vector of binned mz counts for each spot/pixel
#'
#' @examples
#' # mz_values <- plusminus(SeuratObj, 448.2, 0.05)
#' # bin.mz(SeuratObj, mz_values)
bin.mz <- function(data, mz_list, assay = "Spatial", slot = "counts", stored.in.metadata = FALSE){
  data_copy <- data

  if (stored.in.metadata){
    metadata_counts <- data_copy@meta.data[mz_list]
    if (length(colnames(metadata_counts)) < 2) {
      stop("One or more genes not found in the assay meta.data.")
    }
    binned.data <- Matrix::rowSums(metadata_counts)

  } else{
    assay_counts <- data_copy[[assay]][slot]
    selected_genes <- assay_counts[mz_list, , drop = FALSE]
    if (is.null(selected_genes)) {
      stop("One or more genes not found in the assay counts.")
    }
    binned.data <- Matrix::colSums(selected_genes)
  }

  return(binned.data)
}


#' Identifies all mz peaks within a plus-minus range of the target_mz
#'
#'
#' @param data Seurat Spatial Metabolomic object containing mz values
#' @param target_mz Numeric value defining the target m/z peak
#' @param plus_minus Numeric value defining the range/threshold either side of the target peak to also be identified
#'
#' @returns Vector of m/z values within plus-minus range of the target mz value
#'
#' @examples
#' # plusminus(SeuratObj, target_mz = 400.01, plus_minus = 0.005)
plusminus <- function(data, target_mz, plus_minus){
  feature_list <- c()
  center <- FindNearestMZ(data, target_mz)

  center_value <- as.numeric(gsub("mz-", "", center))
  up_value <- center_value + as.numeric(plus_minus)
  low_value <- center_value - as.numeric(plus_minus)

  upper <- FindNearestMZ(data, up_value)
  lower <- FindNearestMZ(data, low_value)

  feature_list <- c(feature_list, center)
  if (!(upper %in% feature_list)){
    feature_list <- c(feature_list, upper)
  }
  if (!(lower %in% feature_list)){
    feature_list <- c(feature_list, lower)
  }

  if (length(feature_list) == 1){
    warning("No other m/z peaks in plusminus range -> increase to include more peaks")
  }
  return(feature_list)
}



#' Helper Function to generate merged counts within the plus minus range provided
#'
#' @param object Seurat Spatial Metabolomic object containing mz values.
#' @param mz_list Vector of numeric m/z values to plot (e.g. c(mz-400.1578, mz-300.1)).
#' @param plusminus Numeric value defining the range/threshold either side of each mz peak provided.
#'
#' @return List contatining a Seurat data object which contains the binned counts in the `@metadata` slot, the column name under where the data is stored and the title of this binned which will be plotted.
#'
#' @examples
#' ### HELPER FUNCTION ###
plot_plus_minus <-function(object, mz_list, plusminus){

  data_copy <- object
  col_names_to_plot <- c()
  plot_titles <- c()

  for (target_mz in mz_list){
    mz_integer <- as.numeric(strsplit(target_mz, "-")[[1]][2])

    meta_col <- bin.mz(data_copy, plusminus(data_copy, mz_integer, plusminus))

    col_name <- paste0(target_mz,"_plusminus_", plusminus)
    plot_name <- paste0("mz: ", round(mz_integer, 3)," \u00b1 ", plusminus)

    col_names_to_plot <- c(col_names_to_plot,col_name)
    plot_titles <- c(plot_titles,plot_name)
    data_copy[[col_name]] = meta_col
  }

  return(list( "data_copy" = data_copy,
               "col_names_to_plot" = col_names_to_plot,
               "plot_titles" = plot_titles))

}



#' Helper function to determine if a column contains categorical or continuous Data
#'
#' @param col data.frame column containing data of interest.
#'
#' @return Character string defining the data type contained withing the column
#'
#' @examples
#' #HELPER FUNCTION
check_column_type <- function(col) {
  if (is.factor(col) || is.character(col)) {
    return("Categorical")
  } else if (is.numeric(col)) {
    return("Continuous")
  } else {
    return("Unknown")
  }
}


#' Helper function for converting Seurat Class ggplots from spot to pixel layout
#'
#' @param plot ggplot object contating the doplot to be converted into pixel layout
#'
#' @return plot where spots are in pixel layout rather then spot
#' @export
#'
#' @examples
#' #pixelPlot(SpatialFeaturPlot(SpaMTP.obj, features = "nFeature_Spatial"))
pixelPlot <- function(plot){

  plots <- lapply(1:length(plot), function (i){

    image_data <- plot[[i]]$data
    data_col <- rlang::quo_text(plot[[i]]$mapping$fill)
    image_metadata <- ggplot_build(plot[[i]])
    size <- unique(image_metadata$data[[1]]$size)
    image_data$fill <- image_metadata$data[[1]]$fill

    if (check_column_type(image_data[[data_col]]) == "Categorical"){
      palette <- unique(image_data$fill)
      names(palette) <- unlist(unique(image_data[[data_col]]))

      custom_plot <- ggplot(image_data, aes(x = y, y = x)) +
        geom_point(aes(color = !!rlang::sym(data_col)), shape = 15, size = size) +  # Customize point size and appearance  # Define size range
        labs(title = plot[[i]]$labels$title, color = plot[[i]]$labels$title) & theme_void() &
        theme(plot.title = element_text(hjust = 0.5)) & scale_color_manual(values = palette)

    } else {
      custom_plot <- ggplot(image_data, aes(x = y, y = x)) +
        geom_point(aes(color = !!rlang::sym(data_col)), shape = 15, size = size) +  # Customize point size and appearance  # Define size range
        labs(title = plot[[i]]$labels$title, color = plot[[i]]$labels$title) & theme_void() &
        theme(plot.title = element_text(hjust = 0.5)) & scale_color_gradientn(colors = unique(image_data[order(image_data[[data_col]]),]$fill))
    }

    custom_plot

  })

  return(purrr::reduce(plots, `+`))
}




#' Plot expression of m/z values spatially
#'
#' This is for plotting SM data from a SpaMTP Seurat Object that does not contain an image(i.e. H&E image).
#' This function inherits from Seurat::ImageFeaturePlot().
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param mzs Vector of numeric m/z values to plot (e.g. c(400.1578, 300.1)). The function FindNearestMZ() is used to automatically find the nearest m/z value to the ones given.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param fov Character string of name of FOV to plot (default = NULL).
#' @param boundaries A vector of segmentation boundaries per image to plot (default = NULL).
#' @param cols Vector of character strings defining colours used for plotting (default = c("lightgrey", "firebrick1")).
#' @param size Numeric value defining the point size for cells/spots when plotting (default = 0.5).
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param split.by Character string defining a factor in the Seurat Object metadata to split the feature plot by (default = NULL).
#' @param molecules Vector of character strings describing molecules to plot (default = NULL).
#' @param mols.size Numeric value for the point size of molecules to plot (default = 0.1).
#' @param mols.cols A vector of colours for the molecules to plot (default = NULL).
#' @param nmols Integer of the max number of each molecule specified in 'molecules' to be plot (default = 1000).
#' @param alpha Numeric value between 0 and 1 defining the image alpha (default = 1).
#' @param border.color Character string specifying the colour of each spot/cell border (default = white).
#' @param border.size Numeric value for the thickness of the cell segmentation border (default = NULL).
#' @param dark.background Boolean value indicating if the plot background is coloured black (default = FALSE).
#' @param blend Boolean value indicating whether to scale and blend expression values to visualize coexpression of two features (default = FALSE).
#' @param blend.threshold Numeric value defining the color cutoff from weak signal to strong signal; ranges from 0 to 1 (default = 0.5).
#' @param crop Boolean value of whether to crop the plots to area with cells only (default = FALSE).
#' @param cells Vector of character strings defining a group of cells to plot (default = NUll; plots all cells).
#' @param scale Set color scaling across multiple plots; c("features", "all", "none").
#' @param overlap Overlay boundaries from a single image to create a single plot (default = FALSE).
#' @param axes Boolean defining if to keep axes and panel background (default = FALSE).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param coord.fixed Boolean value of it to plot cartesian coordinates with fixed aspect ratio (default = TURE).
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "Spatial").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param plot.pixel Boolean indicating if the plot should display pixel square shapes, if false will plot with spots (deafult = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified m/z values
#' @export
#'
#' @examples
#' # ImageMZPlot(SeuratObj, mzs = c(400.678, 300.1))
#' # ImageMZPlot(SeuratObj, mzs = c(400.678, 300.1), plusminus = 0.05)
ImageMZPlot <- function(object,
                        mzs,
                        plusminus = NULL,
                        fov = NULL,
                        boundaries = NULL,
                        cols = if (isTRUE(x = blend)) {
                          c("lightgrey", "#ff0000", "#00ff00")
                        } else {
                          c("lightgrey", "firebrick1")
                        },
                        size = 0.5,
                        min.cutoff = NA,
                        max.cutoff = NA,
                        split.by = NULL,
                        molecules = NULL,
                        mols.size = 0.1,
                        mols.cols = NULL,
                        nmols = 1000,
                        alpha = 1,
                        border.color = "white",
                        border.size = NULL,
                        dark.background = TRUE,
                        blend = FALSE,
                        blend.threshold = 0.5,
                        crop = FALSE,
                        cells = NULL,
                        scale = c("feature", "all", "none"),
                        overlap = FALSE,
                        axes = FALSE,
                        combine = TRUE,
                        coord.fixed = TRUE,
                        assay = "Spatial",
                        slot = "data",
                        plot.pixel = FALSE,
                        verbose = TRUE
){

  if (is.null(mzs)){
    stop("No mz values have been supplied")
  } else{

    mz_list <- c()
    for (target_mz in mzs){
      mz_string <- FindNearestMZ(object, target_mz, assay = assay)
      mz_list <- c(mz_list, mz_string)
    }
  }

  if (is.null(assay)){

    verbose_message(message_text =  "Seurat Assay set to NULL, default assay will be used", verbose = verbose)

    assay <- DefaultAssay(object)
  } else {
    DefaultAssay(object) <- assay
  }

  if (is.null(slot)){
    verbose_message(message_text =  "Default data slot will be used. If data slot is not present will use counts slot instead", verbose = verbose)

  } else if (slot != "data"){

    object[[assay]]$data <- object[[assay]][slot]
  }

  if (!(is.null(plusminus))){

    pl.plustmin.data <- plot_plus_minus(object, mz_list, plusminus)

    data_copy <- pl.plustmin.data$data_copy
    col_names_to_plot <- pl.plustmin.data$col_names_to_plot
    plot_titles <- pl.plustmin.data$plot_titles

    plot <- Seurat::ImageFeaturePlot(object = data_copy,
                                     features = col_names_to_plot,
                                     fov = fov,
                                     boundaries = boundaries,
                                     cols = cols,
                                     size = size,
                                     min.cutoff = min.cutoff,
                                     max.cutoff = max.cutoff,
                                     split.by = split.by,
                                     molecules = molecules,
                                     mols.size = mols.size,
                                     mols.cols = mols.cols,
                                     nmols = nmols,
                                     alpha = alpha,
                                     border.color = border.color,
                                     border.size = border.size,
                                     dark.background = dark.background,
                                     crop = crop,
                                     cells = cells,
                                     scale = scale,
                                     overlap = overlap,
                                     axes = axes,
                                     combine = combine,
                                     coord.fixed = coord.fixed
    )

    for (plot_idx in seq(1, length(plot_titles))){
      plot[[plot_idx]] <- plot[[plot_idx]] &
        ggplot2::ggtitle(plot_titles[[plot_idx]]) &
        ggplot2::labs(fill = plot_titles[[plot_idx]])
    }

  } else {

    plot <- Seurat::ImageFeaturePlot(object = object,
                                     features = mz_list,
                                     fov = fov,
                                     boundaries = boundaries,
                                     cols = cols,
                                     size = size,
                                     min.cutoff = min.cutoff,
                                     max.cutoff = max.cutoff,
                                     split.by = split.by,
                                     molecules = molecules,
                                     mols.size = mols.size,
                                     mols.cols = mols.cols,
                                     nmols = nmols,
                                     alpha = alpha,
                                     border.color = border.color,
                                     border.size = border.size,
                                     dark.background = dark.background,
                                     crop = crop,
                                     cells = cells,
                                     scale = scale,
                                     overlap = overlap,
                                     axes = axes,
                                     combine = combine,
                                     coord.fixed = coord.fixed
    )

    for (plot_idx in seq(1, length(mz_list))){
      mz_integer <- as.numeric(strsplit(mz_list[plot_idx], "-")[[1]][2])
      plot[[plot_idx]] <- plot[[plot_idx]] &
        ggplot2::ggtitle(paste0("mz: ",round(mz_integer,3)))&
        ggplot2::labs(fill = paste0("mz: ",round(mz_integer,3)))
    }
  }

  if (plot.pixel){
    plot <- pixelPlot(plot)
  }

  return(plot)

}



#' Plot expression of annotated metabolites spatially
#'
#' This is for plotting annotated SM data from a SpaMTP Seurat Object that does not contain an image(i.e. H&E image).
#' This function inputs the metabolite name such as 'Glutamine' rather then a specific m/z value.
#' For m/z values please use `ImageMZPlot()`.This function inherits off Seurat::ImageFeaturePlot().
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param metabolites Vector of metabolite names to plot (e.g. c("Glucose", "Glutamine")). The Seurat Object provided must contain annotations in the respective assay metadata.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param fov Character string of name of FOV to plot (default = NULL).
#' @param boundaries A vector of segmentation boundaries per image to plot (default = NULL).
#' @param cols Vector of character strings defining colours used for plotting (default = c("lightgrey", "firebrick1")).
#' @param size Numeric value defining the point size for cells/spots when plotting (default = 0.5).
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param split.by Character string defining a factor in the Seurat Object metadata to split the feature plot by (default = NULL).
#' @param molecules Vector of character strings describing molecules to plot (default = NULL).
#' @param mols.size Numeric value for the point size of molecules to plot (default = 0.1).
#' @param mols.cols A vector of colours for the molecules to plot (default = NULL).
#' @param nmols Integer of the max number of each molecule specified in 'molecules' to be plot (default = 1000).
#' @param alpha Numeric value between 0 and 1 defining the spot alpha (default = 1).
#' @param border.color Character string specifying the colour of each spot/cell border (default = white).
#' @param border.size Numeric value for the thickness of the cell segmentation border (default = NULL).
#' @param dark.background Boolean value indicating if the plot background is coloured black (default = FALSE).
#' @param blend Boolean value indicating whether to scale and blend expression values to visualize coexpression of two features (default = FALSE).
#' @param blend.threshold Numeric value defining the color cutoff from weak signal to strong signal; ranges from 0 to 1 (default = 0.5).
#' @param crop Boolean value of whether to crop the plots to area with cells only (default = FALSE).
#' @param cells Vector of character strings defining a group of cells to plot (default = NUll; plots all cells).
#' @param scale Set color scaling across multiple plots; c("features", "all", "none").
#' @param overlap Overlay boundaries from a single image to create a single plot (default = FALSE).
#' @param axes Boolean defining if to keep axes and panel background (default = FALSE).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param coord.fixed Boolean value of it to plot cartesian coordinates with fixed aspect ratio (default = TURE).
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "Spatial").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param column.name Character string defining the column name where the annotations are stored in the slot meta.data (default = "all_IsomerNames").
#' @param plot.exact Boolean value describing if to only plot exact matches to the metabolite search terms, else will plot all metabolites which contain serach word in name (default = TRUE).
#' @param plot.pixel Boolean indicating if the plot should display pixel square shapes, if false will plot with spots (deafult = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified metabolites.
#' @export
#'
#' @examples
#' # ImageMZPlot(SeuratObj, mzs = c("Glucose", "Glutamine"))
#' # ImageMZPlot(SeuratObj, mzs = c("Glucose", "Glutamine"), plusminus = 0.05)
ImageMZAnnotationPlot <- function(object,
                                  metabolites,
                                  plusminus = NULL,
                                  fov = NULL,
                                  boundaries = NULL,
                                  cols = if (isTRUE(x = blend)) {
                                    c("lightgrey", "#ff0000", "#00ff00")
                                  } else {
                                    c("lightgrey", "firebrick1")
                                  },
                                  size = 0.5,
                                  min.cutoff = NA,
                                  max.cutoff = NA,
                                  split.by = NULL,
                                  molecules = NULL,
                                  mols.size = 0.1,
                                  mols.cols = NULL,
                                  nmols = 1000,
                                  alpha = 1,
                                  border.color = "white",
                                  border.size = NULL,
                                  dark.background = TRUE,
                                  blend = FALSE,
                                  blend.threshold = 0.5,
                                  crop = FALSE,
                                  cells = NULL,
                                  scale = c("feature", "all", "none"),
                                  overlap = FALSE,
                                  axes = FALSE,
                                  combine = TRUE,
                                  coord.fixed = TRUE,
                                  assay = "Spatial",
                                  slot = "data",
                                  column.name = "all_IsomerNames",
                                  plot.exact = TRUE,
                                  plot.pixel = FALSE,
                                  verbose = TRUE

){

  if (is.null(assay)){
    stop("Seurat Assay set to NULL, please provide correct assay name")
  }

  multi_annotation <- FALSE
  multi_plusminus <- NULL
  mzs <- c()
  for (metabolite in metabolites){
    met.row <- SearchAnnotations(object,metabolite, assay = assay, column.name = column.name, search.exact = plot.exact)

    if (dim(met.row)[1] == 0){
      warning(paste("There are no entries for the metabolite: ", metabolite, " in the current annotation metadata ...",
                    "\n please refine your search term. You can check for annotations using SearchAnnotations()"))
      stop("n entries must be > 1")

    } else if (dim(met.row)[1] > 1){
      warning(paste("There are multiple m/z values assigned to the metabolite: ", metabolite,
                    "\n NOTE: Data from all m/z values will be merged ... "))
      multi_annotation <- TRUE
      mz_values <- c()

      for (row in 1:dim(met.row)[1]){

        all_annots <- unlist(strsplit(met.row[[column.name]][row], "; "))

        if (length(all_annots) != 1){
          warning(paste("There are another ", (length(all_annots)-1),
                        "annotations assocated with this metabolite ( ",metabolite,
                        " )... These being: ",
                        "\n", paste(list(all_annots)), "\n Use SearchAnnotations() to see ..."))
        }

        mz <- met.row$mz_names[row]
        mz_values <- c(mz_values, mz)
      }


      data_copy <- object
      col_names_to_plot <- mz_values

      stored.in.metadata <- FALSE #this is used to pull from metadata for multi-plusminus plots


      if (!(is.null(plusminus))){
        pl.plustmin.data <- plot_plus_minus(object, mz_values, plusminus)

        data_copy <- pl.plustmin.data$data_copy
        col_names_to_plot <- pl.plustmin.data$col_names_to_plot
        plot_titles <- pl.plustmin.data$plot_titles
        stored.in.metadata <- TRUE
      }

      binned.data <- bin.mz( data_copy, col_names_to_plot, stored.in.metadata = stored.in.metadata, assay = assay, slot = slot)

      col_name <- paste0(metabolite,"_binned")

      mzs <- c(mzs,col_name)

      data_copy[[col_name]] = binned.data
      multi_plusminus <- plusminus
      plusminus <- NULL



    } else {
      all_annots <- unlist(strsplit(met.row[[column.name]], "; "))

      if (length(all_annots) != 1){
        warning(paste("There are another ", (length(all_annots)-1),
                      "annotations assocated with this metabolite ( ",metabolite,
                      " )... These being: ",
                      "\n", paste(list(all_annots)), "\n Use SearchAnnotations() to see ..."))
      }

      mz <- met.row$raw_mz
      mzs <- c(mzs,mz)
      data_copy <- object
    }
  }


  if (multi_annotation){

    DefaultAssay(data_copy) <- assay

    if (is.null(slot)){
      verbose_message(message_text =  "Default data slot will be used. If data slot is not present will use counts slot instead", verbose = verbose)

    } else if (slot != "data"){

      data_copy[[assay]]$data <- data_copy[[assay]][slot]
    }

    plot <- Seurat::ImageFeaturePlot(data_copy,
                                     features =  mzs,
                                     fov = fov,
                                     boundaries = boundaries,
                                     cols = cols,
                                     size = size,
                                     min.cutoff = min.cutoff,
                                     max.cutoff = max.cutoff,
                                     split.by = split.by,
                                     molecules = molecules,
                                     mols.size = mols.size,
                                     mols.cols = mols.cols,
                                     nmols = nmols,
                                     alpha = alpha,
                                     border.color = border.color,
                                     border.size = border.size,
                                     dark.background = dark.background,
                                     blend = blend,
                                     blend.threshold = blend.threshold,
                                     crop = crop,
                                     cells = cells,
                                     scale = scale,
                                     overlap = overlap,
                                     axes = axes,
                                     combine = combine,
                                     coord.fixed = coord.fixed)
  } else {
    plot <- ImageMZPlot(data_copy,
                        mzs,
                        plusminus = plusminus,
                        fov = fov,
                        boundaries = boundaries,
                        cols = cols,
                        size = size,
                        min.cutoff = min.cutoff,
                        max.cutoff = max.cutoff,
                        split.by = split.by,
                        molecules = molecules,
                        mols.size = mols.size,
                        mols.cols = mols.cols,
                        nmols = nmols,
                        alpha = alpha,
                        border.color = border.color,
                        border.size = border.size,
                        dark.background = dark.background,
                        blend = blend,
                        blend.threshold = blend.threshold,
                        crop = crop,
                        cells = cells,
                        scale = scale,
                        overlap = overlap,
                        axes = axes,
                        combine = combine,
                        coord.fixed = coord.fixed,
                        assay = assay,
                        slot = slot,
                        verbose = verbose)
  }
  for (i in 1:length(plot)){

    plusmin_str <- ""

    if (!(is.null(multi_plusminus))){
      plusmin_str <- paste0(" \u00b1 ", multi_plusminus)
    }

    if (!(is.null(plusminus))){
      plusmin_str <- paste0(" \u00b1 ", plusminus)
    }

    plot[[i]] <- plot[[i]] &
      ggplot2::ggtitle(paste0(metabolites[i], plusmin_str))&
      ggplot2::labs(fill = paste0(metabolites[i], plusmin_str))

  }

  if (plot.pixel){
    plot <- pixelPlot(plot)
  }

  return(plot)

}



#' Plot expression of m/z values spatially for a Spatial SpaMTP Seurat Objects.
#'
#' This is for plotting SM m/z data from a SpaMTP Seurat Object that contains an image (i.e. H&E image).
#' This function inherits off Seurat::SpatialFeaturePlot(). Look here for more detailed documentation about inputs.
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param mzs Vector of numeric m/z values to plot (e.g. c(400.1578, 300.1)). The function FindNearestMZ() is used to automatically find the nearest m/z value to the ones given.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param images Character string of the name of the image to plot (default = NULL).
#' @param crop Boolean value indicating if to crop the plot to focus on only points being plotted (default = TRUE).
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "Spatial").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param keep.scale Character string describing how to handle the color scale across multiple plots. Check Seurat::SpatialFeaturePlot() for all options (default = "feature").
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param ncol Integer defining the number of columns if plotting multiple plots (default = NULL).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param pt.size.factor Numeric value defining the point size for spots when plotting (default = 1.6).
#' @param alpha Numeric value between 0 and 1 defining the spot alpha (default = 1).
#' @param image.alpha Numeric value between 0 and 1 defining the image alpha (default = 1).
#' @param stroke Numeric value describing the width of the border around the spot (default = 0.25).
#' @param interactive Boolean value of if to launch an interactive SpatialDimPlot or SpatialFeaturePlot session, see Seurat::ISpatialDimPlot() or Seurat::ISpatialFeaturePlot() for more details (default = FALSE).
#' @param information An optional dataframe or matrix of extra information to be displayed on hover (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified m/z values for Spatial data with H&E Image
#' @export
#'
#' @examples
#' # SpatialMZPlot(SeuratObj, mzs = c(400.678, 300.1))
#' # SpatialMZPlot(SeuratObj, mzs = c(400.678, 300.1), plusminus = 0.05)
SpatialMZPlot <- function(object,
                          mzs,
                          plusminus = NULL,
                          images = NULL,
                          crop = TRUE,
                          assay = "Spatial",
                          slot = "counts",
                          keep.scale = "feature",
                          min.cutoff = NA,
                          max.cutoff = NA,
                          ncol = NULL,
                          combine = TRUE,
                          pt.size.factor = 1.6,
                          alpha = c(1, 1),
                          image.alpha = 1,
                          stroke = 0.25,
                          interactive = FALSE,
                          information = NULL,
                          verbose = TRUE
){

  if (is.null(mzs)){
    stop("No mz values have been supplied")
  } else{

    mz_list <- c()
    for (target_mz in mzs){
      mz_string <- FindNearestMZ(object, target_mz, assay = assay)
      mz_list <- c(mz_list, mz_string)
    }
  }

  if (is.null(assay)){
    verbose_message(message_text =  "Seurat Assay set to NULL, default assay will be used", verbose = verbose)
    assay <- DefaultAssay(object)
  } else {
    DefaultAssay(object) <- assay
  }

  if (!(is.null(plusminus))){

    data_copy <- object
    col_names_to_plot <- c()
    plot_titles <- c()

    for (target_mz in mz_list){
      mz_integer <- as.numeric(strsplit(target_mz, "-")[[1]][2])

      meta_col <- bin.mz(data_copy, plusminus(data_copy, mz_integer, plusminus), assay = assay, slot = slot)

      col_name <- paste0(target_mz,"_plusminus_", plusminus)
      plot_name <- paste0("mz: ", round(mz_integer, 3)," \u00b1 ", plusminus)

      col_names_to_plot <- c(col_names_to_plot,col_name)
      plot_titles <- c(plot_titles,plot_name)
      data_copy[[col_name]] = meta_col
    }

    plot <- Seurat::SpatialFeaturePlot(object = data_copy,
                                       features = col_names_to_plot,
                                       images = images,
                                       crop = crop,
                                       slot = slot,
                                       keep.scale = keep.scale,
                                       min.cutoff = min.cutoff,
                                       max.cutoff = max.cutoff,
                                       ncol = ncol,
                                       combine = combine,
                                       pt.size.factor = pt.size.factor,
                                       alpha = alpha,
                                       image.alpha = image.alpha,
                                       stroke = stroke,
                                       interactive = interactive,
                                       information = information
    )

    for (plot_idx in seq(1, length(plot_titles))){
      plot[[plot_idx]] <- plot[[plot_idx]] &
        ggplot2::labs(fill = plot_titles[[plot_idx]])
    }

  } else {

    plot <- Seurat::SpatialFeaturePlot(object = object,
                                       features = mz_list,
                                       images = images,
                                       crop = crop,
                                       slot = slot,
                                       keep.scale = keep.scale,
                                       min.cutoff = min.cutoff,
                                       max.cutoff = max.cutoff,
                                       ncol = ncol,
                                       combine = combine,
                                       pt.size.factor = pt.size.factor,
                                       alpha = alpha,
                                       image.alpha = image.alpha,
                                       stroke = stroke,
                                       interactive = interactive,
                                       information = information
    )

    for (plot_idx in seq(1, length(mz_list))){
      mz_integer <- as.numeric(strsplit(mz_list[plot_idx], "-")[[1]][2])
      plot[[plot_idx]] <- plot[[plot_idx]] &
        ggplot2::labs(fill = paste0("mz: ",round(mz_integer,3)))
    }
  }


  return(plot)

}



#' Plot expression of metabolites in spatially from a Spatial SpaMTP Objects.
#'
#' This is for plotting values of annotated metabolites stored in a SpaMTP Seurat Object containing an image(i.e. H&E image).
#' This function inherits off Seurat::ImageFeaturePlot(). Look here for more detailed documentation about inputs.
#'
#' @param object Seurat Spatial Metabolomic Object to Visualise.
#' @param metabolites Vector of metabolite names to plot (e.g. c("Glucose", "Glutamine")). The Seurat Object provided must contain annotations in the respective assay metadata.
#' @param plusminus Numeric value defining the range/threshold either side of the target peak/peaks to be binned together for plotting (default = NULL).
#' @param images Character string of the name of the image to plot (default = NULL).
#' @param crop Boolean value indicating if to crop the plot to focus on only points being plotted (default = TRUE).
#' @param assay Character string indicating which Seurat object assay to pull data form (default = "Spatial").
#' @param slot Character string indicating the assay slot to use to pull expression values form (default = "counts").
#' @param keep.scale Character string describing how to handle the color scale across multiple plots. Check Seurat::SpatialFeaturePlot() for all options (default = "feature").
#' @param min.cutoff Vector of numeric value describing the minimum cutoff values for each m/z feature (default = NA).
#' @param max.cutoff Vector of numeric value describing the maximum cutoff values for each m/z feature (default = NA).
#' @param ncol Integer defining the number of columns if plotting multiple plots (default = NULL).
#' @param combine Boolean value stating if to combine plots into a single patchworked ggplot object (default = TRUE).
#' @param pt.size.factor Numeric value defining the point size for spots when plotting (default = 1.6).
#' @param alpha Numeric value between 0 and 1 defining the spot alpha (default = 1).
#' @param image.alpha Numeric value between 0 and 1 defining the image alpha (default = 1).
#' @param stroke Numeric value describing the width of the border around the spot (default. =0.25).
#' @param interactive Boolean value of if to launch an interactive SpatialDimPlot or SpatialFeaturePlot session, see Seurat::ISpatialDimPlot() or Seurat::ISpatialFeaturePlot() for more details (default = FALSE).
#' @param information An optional dataframe or matrix of extra infomation to be displayed on hover (default = NULL).
#' @param column.name Character string defining the column name where the annotations are stored in the slot meta.data (default = "all_IsomerNames").
#' @param plot.exact Boolean value describing if to only plot exact matches to the metabolite search terms, else will plot all metabolites which contain serach word in name (default = TRUE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A ggplot showing the Spatial representation of expression data of specified m/z values for Spatial data with H&E Image
#' @export
#'
#' @examples
#' # SpatialMZAnnotationPlot(SeuratObj, mzs = c("Glucose", "Glutamine"))
#' # SpatialMZAnnotationPlot(SeuratObj, mzs = c("Glucose", "Glutamine"), plusminus = 0.05)
SpatialMZAnnotationPlot <- function(object,
                                    metabolites,
                                    plusminus = NULL,
                                    images = NULL,
                                    crop = TRUE,
                                    assay = "Spatial",
                                    slot = "counts",
                                    keep.scale = "feature",
                                    min.cutoff = NA,
                                    max.cutoff = NA,
                                    ncol = NULL,
                                    combine = TRUE,
                                    pt.size.factor = 1.6,
                                    alpha = c(1, 1),
                                    image.alpha = 1,
                                    stroke = 0.25,
                                    interactive = FALSE,
                                    information = NULL,
                                    column.name = "all_IsomerNames",
                                    plot.exact = TRUE,
                                    verbose = TRUE

){

  if (is.null(assay)){
    stop("Seurat Assay set to NULL, please provide correct assay name")
  }

  mzs <- c()
  for (metabolite in metabolites){
    met.row <- SearchAnnotations(object,metabolite, assay = assay, column.name = column.name,search.exact = plot.exact)

    if (dim(met.row)[1] != 1){
      warning(paste("There are either none or multiple entries for the metabolite: ", metabolite,
                    "\n please check using SearchAnnotations() and FindDuplicateAnnotations"))
      stop("n entries != 1")
    }

    all_annots <- unlist(strsplit(met.row[[column.name]], "; "))

    if (length(all_annots) != 1){
      warning(paste("There are another ", (length(all_annots)-1),
                    "annotations assocated with this metabolite ( ",metabolite," )... These being: ",
                    "\n", paste(list(all_annots)), "\n Use SearchAnnotations() to see ..."))
    }

    mz <- met.row$raw_mz
    mzs <- c(mzs,mz)
  }

  plot <- SpatialMZPlot(object,
                        mzs,
                        plusminus = plusminus,
                        images = images,
                        crop = crop,
                        assay = assay,
                        slot = slot,
                        keep.scale = keep.scale,
                        min.cutoff = min.cutoff,
                        max.cutoff = max.cutoff,
                        ncol = ncol,
                        combine = combine,
                        pt.size.factor = pt.size.factor,
                        alpha = alpha,
                        image.alpha = image.alpha,
                        stroke = stroke,
                        interactive = interactive,
                        information = information,
                        verbose = verbose
  )

  for (i in 1:length(plot)){

    plusmin_str <- ""

    if (!(is.null(plusminus))){
      plusmin_str <- paste0(" \u00b1 ", plusminus)
    }

    plot[[i]] <- plot[[i]] &
      ggplot2::labs(fill = paste0(metabolites[i], plusmin_str))

  }

  return(plot)

}


#' Plot mass intensity spectra
#'
#' This function plots mean mass spectra intensity values for a given SpaMTP Seurat Object, and groups/splits by categories if supplied.
#'
#' @param data Seurat object containing data to be plot.
#' @param group.by Character string defining the name of the meta.data column to group the data by. Results from each group will be overlayed on the one plot (default = NULL).
#' @param split.by Character string defining the name of the meta.data column to group and split the data by. Results from each group will be plotted individually (default = NULL).
#' @param cols Vector of character strings defining the colours to be used to identify each group (default = NULL).
#' @param assay Character string defining the relative Seurat Object assay to pull the required intensity data from (default = "Spatial").
#' @param slot Character string defining the relative slot from the Seurat Object assay to pull the required intensity data from (default = "counts").
#' @param label.annotations Boolean value defining whether to plot metabolite annotations of the supplied mz.labels or metabolite.labels on the plot (default = FALSE).
#' @param annotation.column Character string defining the name of the feature meta.data column which contains the stored m/z annotations (default = "all_IsomerNames").
#' @param mz.labels Vector of character strings defining the m/z values to display on the plot (default = NULL).
#' @param metabolite.labels Vector of character strings defining the metabolite names to display on the plot (default = NULL).
#' @param xlab Character string describing the x-axis title (default = "m/z").
#' @param ylab Character string describing the y-axis title (default = "intensity").
#' @param mass.range Vector of numeric values defining the range of m/z values to include (default = NULL).
#' @param y.lim Vector of numeric values defining the range of intensity values to include (default = NULL).
#' @param labelCex Numeric values for character expansion factor. Seen graphics::text() for more details (default = 0.7).
#' @param labelAdj One or two values in `[0,1]` which specify the x and y adjustments for the label. Seen graphics::text() for more details (default = NULL).
#' @param labelOffset Value that controls the distance of the text label from the specified coordinates. Seen graphics::text() for more details (default = 0).
#' @param labelCol Character string defining the colour of the annotation labels (default = "#eb4034").
#' @param nlabels.to.show Numeric value defining the number of annotations to show per m/z (default = NULL).
#'
#' @return A mass spectrometry plot displaying mean intensity values
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' ## Plot mean of whole tissue section
#' # MassIntensityPlot(SeuratObj)
#'
#' ## Plot ssc segmentation groups on the same plot
#' # MassIntensityPlot(SeuratObj, group.by = "ssc")
#'
#' ## Plot mean of each ssc segmentation of separate plot with mz annotations
#' # MassIntensityPlot(SeuratObj, split.by= "ssc", mz.labels = c(329.166), plot.layout = c(5,2))
#'
#' ## Plot mean of each ssc segmentation of separate plot with metabolite annotations
#' # MassIntensityPlot(SeuratObj, split.by= "ssc", mz.labels = c(329.166), label.annotations = TRUE)
MassIntensityPlot <- function (data,
                               group.by = NULL,
                               split.by = NULL,
                               cols = NULL,
                               assay = "Spatial",
                               slot = "counts",
                               label.annotations = FALSE,
                               annotation.column = "all_IsomerNames",
                               mz.labels = NULL,
                               metabolite.labels = NULL,
                               xlab = "m/z",
                               ylab = "intensity",
                               mass.range = NULL,
                               y.lim = NULL,
                               labelCex = 5,
                               labelAdj = -1,
                               labelOffset = 0,
                               labelCol = "#eb4034",
                               nlabels.to.show = NULL){

  if (!(is.null(group.by))&!(is.null(split.by))){
    stop("'group.by' and 'split.by' cannot both be valid -> pick only one option to set = idents")
  }
  if (!(is.null(mz.labels))&!(is.null(metabolite.labels))){
    stop("'mz.labels' and 'metabolite.labels' cannot both be valid -> pick only one option to set = c(labels)")
  }

  if (!(is.null(annotation.column))){
    if (!(annotation.column %in% colnames(data[[assay]]@meta.data))){
      warning(paste("'",annotation.column," column not in object metadata. If data object does not have annotations set annotation.column = NULL"))
      stop("annotation.column does not exist")
    } else {
      if (!is.null(nlabels.to.show)){
        data[[assay]]@meta.data[[annotation.column]] <- labels_to_show(data[[assay]]@meta.data[[annotation.column]], n = nlabels.to.show)
      }
    }
  }

  if (!(is.null(group.by))|!(is.null(split.by))){

    metadata.column <- ifelse(!(is.null(group.by)), group.by, split.by)

    if (!(metadata.column %in% colnames(data@meta.data))){
      warning(paste("'",metadata.column," column not in object metadata. Pick and approprite group.by or split.by column"))
      stop("metadata columne does not exist")
    }


    run <- unique(data@meta.data[[metadata.column]])

    SeuratObject::Idents(data) <- metadata.column
    data_list <- list()
    for (ident in unique(run)){
      cells <- Seurat::WhichCells(data, idents = ident)


      if (length(cells) > 1){
        mean_data <- Matrix::rowMeans(data[[assay]][slot][,cells])
      } else {
        mean_data <- data[[assay]][slot]
      }

      data_list[[ident]] <- mean_data
    }

    means <- Matrix::as.matrix(as.data.frame(data_list))
    x_coord <- c(1:length(unique(run)))

  } else {
    means <- Matrix::as.matrix(Matrix::rowMeans(data[[assay]][slot]))
    run <- factor("Sample_Mean")
    x_coord <- c(1)
  }

  colnames(means) <- NULL

  mzs <- unlist(lapply(rownames(data[[assay]]@features), function(x) as.numeric(stringr::str_split(x, "mz-")[[1]][2])))

  if (!(is.null(mz.labels))){

    labels <- data[[assay]]@meta.data$raw_mz
    mz_list <- c()
    for (target_mz in mz.labels){
      mz_string <- FindNearestMZ(data, target_mz, assay = assay)
      mz_list <- c(mz_list, stringr::str_split(pattern = "mz-", string = mz_string)[[1]][2])
    }

    matching <- labels %in% as.numeric(mz_list)

    labels[!matching] <- NA


  } else if (!(is.null(metabolite.labels))){

    labels <- data[[assay]]@meta.data$raw_mz
    mzs <- c()
    for (metabolite in metabolite.labels){
      met.row <- SearchAnnotations(data, metabolite, assay = assay, column.name = annotation.column, search.exact = TRUE)

      if (dim(met.row)[1] != 1){
        warning(paste("There are either none or multiple entries for the metabolite: ", metabolite,
                      "\n please check using SearchAnnotations() and FindDuplicateAnnotations"))
        stop("n entries != 1")
      }

      all_annots <- unlist(strsplit(met.row[[annotation.column]], "; "))

      if (length(all_annots) != 1){
        warning(paste("There are another ", (length(all_annots)-1),
                      "annotations assocated with this metabolite ( ",metabolite," )... These being: ",
                      "\n", paste(list(all_annots)), "\n Use SearchAnnotations() to see ..."))
      }

      mz <- met.row$raw_mz
      mzs <- c(mzs,mz)

      matching <- labels %in% as.numeric(mzs)

      labels[!matching] <- NA

    }
  } else {
    labels <- rep(NA, length(data[[assay]]@meta.data$raw_mz))
  }


  if (label.annotations){
    feature.metadata <- data[[assay]]@meta.data

    non_na_indices <- !is.na(labels)
    labels[non_na_indices] <- feature.metadata[[annotation.column]][match(labels[non_na_indices], feature.metadata$raw_mz)]

  } else {
    labels <- unlist(lapply(labels, function(x) {if(!is.na(x)){paste0("mz-", x)} else {x}}))
  }

  df <- data.frame(means)
  df$y <- as.numeric(gsub("mz-", "", rownames(df)))
  df$annotation <- labels

  if("means" %in% colnames(df)){
    df$variable <- run
    df$value <- df$means
    df_long <- df[c("y", "annotation", "variable", "value")]
  } else{

    df_long <- df %>%
      pivot_longer(
        cols = starts_with("X"),   # Select columns X1, X2, X3
        names_to = "variable",     # Name for the new column that will store column names (X1, X2, X3)
        values_to = "value"        # Name for the new column that will store the values
      )

    df_long <- df_long %>%
      dplyr::mutate(variable = as.numeric(gsub("X", "", variable)))

    # Replace the indices with the corresponding names
    df_long <- df_long %>%
      dplyr::mutate(variable = ifelse(variable <= length(unique(run)), unique(run)[variable], NA))

    df_long[["variable"]] <- as.character(df_long[["variable"]])
  }

  if (!is.null(mass.range)){
    df_long <- df_long[df_long$y > mass.range[1] & df_long$y < mass.range[2],]
  }


  max_values_df <- df_long %>%
    group_by(y) %>%
    summarize(
      max_value = max(value, na.rm = TRUE),
      annotation = first(annotation),
      .groups = "drop"  # Drop grouping after summarization
    )

  max_values_df <- data.frame(max_values_df)

  if (!is.null(split.by)){
    p <- ggplot(df_long, aes(x = y, y = value, color = variable)) +
      geom_line() +
      geom_text(
        data = max_values_df,
        aes(x = y, y = max_value, label = annotation),
        inherit.aes = FALSE,  # Avoid inheriting `aes` from ggplot
        size = labelCex,  # Avoid inheriting `aes` from ggplot
        vjust = labelAdj,
        color = labelCol,
        hjust = labelOffset
      ) +
      labs(x = xlab, y = ylab, color = "Variable") +
      theme_classic() +
      facet_wrap(~ variable, ncol = 1, scales = "free_y")
  }else{
    p <- ggplot(df_long, aes(x = y, y = value, color = variable)) +
      geom_line() +
      geom_text(
        data = max_values_df,
        aes(x = y, y = max_value, label = annotation),
        inherit.aes = FALSE,
        size = labelCex,  # Avoid inheriting `aes` from ggplot
        vjust = labelAdj,
        color = labelCol,
        hjust = labelOffset

      ) +
      labs(x = xlab, y = ylab, color = "Variable") +
      theme_classic()
  }

  if (!is.null(y.lim)){
    p <- p + ylim(y.lim[1], y.lim[2])
  }

  if(!is.null(cols)){
    p <- p + scale_color_manual(values = cols)
  }

  return(p)
}



#' Generates a 3D spatial feature plot from a SpaMTP object
#'
#' This function generates an interactive plotly 3D spatial plot of 2 specified features.
#' These features can include a combination of metabolites and/or genes expression values, or numerical meta.data values.
#' Features must be numeric values. Please see ` SpaMTP Additional Features` vignette for instructions on how to plot catagorical data, or > 2 features.
#'
#' @param data A SpaMTP Seurat Object.
#' @param features A character vector specifying the features to be plotted.
#' @param assays A character vector specifying the Seurat Object assays to be used for plotting (default = c("SPT", "SPM")).
#' @param slots A character vector specifying the assay slots to be used for plotting (deafult = "counts").
#' @param between.layer.height A numeric value specifying the height between layers (default = 100).
#' @param names A character vector specifying custom names for the features. If NULL will use feature names (default = NULL).
#' @param size A numeric value specifying the size of markers in the plot (default = 3).
#' @param col.palette Character list defining the colour palettes to use for plotting. If only 1 is supplied than it will be used for both variables. Possible palettes to use are: Blackbody,Bluered,Blues,Cividis,Earth,Electric,Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,RdBu,Reds,Viridis,YlGnBu,YlOrRd (default = "Reds").
#' @param x.axis.label A character string specifying the label for the x-axis (default = "x").
#' @param y.axis.label A character string specifying the label for the y-axis (default = "y").
#' @param z.axis.label A character string specifying the label for the z-axis (default = "z").
#' @param show.x.ticks A logical value specifying whether to show ticks on the x-axis (default = FALSE).
#' @param show.y.ticks A logical value specifying whether to show ticks on the y-axis (default = FALSE).
#' @param show.z.ticks A logical value specifying whether to show ticks on the z-axis (default = FALSE).
#' @param show.image Character string specifying the image name to plot. If NULL then no image is plot (default = NULL).
#' @param plot.height Numeric value defining the height of the returned plot (default = 800).
#' @param plot.width Numeric value defining the width of the returned plot (default = 1500).
#' @param image.sf Character string defining the image scalefactor to use (default = "lowres").
#' @param downscale.image Numeric value defining the amount to downscale the image by. This parameter is only necessary when show.image != NULL. Downscaling the image is not necessary, however when rendering plotly to HTML excessively large images can cause issues whereby setting a downscaling value can improve this (default = NULL).
#'
#'
#' @return A 3D Plotly plot
#' @export
#'
#' @import plotly
#'
#' @examples
#' # Plot3DFeature(data = my_data, features = c("gene1", "gene2"), assays = c("SPT", "SPM"), show.image = "slice1")
Plot3DFeature <- function(data,
                          features,
                          assays = c("SPT", "SPM"),
                          slots = "counts",
                          between.layer.height = 100,
                          names= NULL,
                          size = 3,
                          col.palette = "Reds",
                          x.axis.label = "x",
                          y.axis.label = "y",
                          z.axis.label = "z",
                          show.x.ticks = FALSE,
                          show.y.ticks = FALSE,
                          show.z.ticks = FALSE,
                          show.image = NULL,
                          plot.height = 800,
                          plot.width = 1500,
                          image.sf = "lowres",
                          downscale.image = NULL
){

  ## handeling of inncorect input legnths
  if (length(features) < 1 | length(features) > 2){
    stop("Number of features supplied does not match required length. Either 1 or 2 features must be supplied")
  }
  if (length(assays) > 2){
    stop("Number of assays supplied does not match required length. 2 assays or less must be supplied")
  }
  if (length(slots) > 2){
    stop("Number of slots supplied does not match required length. 2 slots or less must be supplied")
  }

  ## Handling of only 1 submitted value
  if (length(features) ==  1 ){
    features <- c(features, features)
  }
  if (length(assays) ==  1 ){
    assays <- c(assays, assays)
  }
  if (length(slots) ==  1 ){
    slots <- c(slots, slots)
  }
  if (length(col.palette) ==  1 ){
    col.palette <- list(col.palette[[1]], col.palette[[1]])
  }
  if (length(names) ==  1 ){
    names <- c(names, names)
  }



  feature_data <- list()

  i <- 1
  default_names <- c()
  for (feature in features) {

    if (feature %in% colnames(data@meta.data)){
      default_names <- c(default_names, feature)
      feature_data[[i]] <- data@meta.data[[feature]]
    } else {
      if (length(assays) != 0){

        feature_data[[i]] <- tryCatch({data[[assays[i]]][slots[i]][feature,]},
                                      error = function(err){
                                        stop("The feature provided does not exist in the ", assays[i],": ",slots[i], " object. Please provide a value feature")})
      } else {
        stop("No assay supplied. The feature ", feature," either does not exist in @meta.data slot, or no matching assay was provided for the gene/feature. Please check SpaMTP object")
      }
      default_names <- c(default_names, paste0(feature,"_", assays[i]))
    }
    i <- i + 1
  }
  if (!(is.null(names))){
    if (length(names) == 2) {
      default_names <- names
    } else {
      warning("length of names must be == 2. Default names will be used instead ... ")
    }
  }


  # Create scatter3d traces for each layer
  trace1 <-
    plot_ly(GetTissueCoordinates(data),
            x = ~x,
            y = ~y,
            z = rep(0 + between.layer.height, times = dim(GetTissueCoordinates(data))[1]),
            type = "scatter3d",
            mode = "markers",
            height = plot.height,
            width = plot.width,
            name = default_names[1],
            marker = list(color = feature_data[[1]],
                          coloraxis = 'coloraxis', size = size)) %>%
    layout(scene = list(
      aspectmode = "data",
      xaxis = list(title = x.axis.label, showticklabels = show.x.ticks),
      yaxis = list(title = y.axis.label, showticklabels = show.y.ticks),
      zaxis = list(title = z.axis.label, showticklabels = show.z.ticks)))

  plot <-  trace1 %>% add_trace(GetTissueCoordinates(data),
                                x = ~x,
                                y = ~y,
                                z = rep(0, times = dim(GetTissueCoordinates(data))[1]),
                                type = "scatter3d",
                                mode = "markers",
                                name = default_names[2],
                                marker = list(color = feature_data[[2]],
                                              coloraxis = 'coloraxis2')) %>%
    layout(autosize = F,
           scene = list(
             aspectmode = "data",
             xaxis = list(title = x.axis.label, showticklabels = show.x.ticks),
             yaxis = list(title = y.axis.label, showticklabels = show.y.ticks),
             zaxis = list(title = z.axis.label, showticklabels = show.z.ticks)),
           coloraxis = list(colorbar = list(orientation = "v",
                                            xanchor ="right",
                                            x = 0,
                                            len = 0.5,
                                            title = list( side = "top",
                                                          text = default_names[1]
                                            )),
                            colorscale = col.palette[[1]]),
           coloraxis2 = list(colorbar = list(orientation = "v",
                                             xanchor ="left",
                                             len = 0.5,
                                             x = 0,
                                             title = list( side = "top",
                                                           text = default_names[2]
                                             )),
                             colorscale = col.palette[[2]])

    )

  if (!is.null(show.image)){

    color_matrix <- as.raster(data@images[[show.image]]@image)
    row_indices <- rep(seq_len(nrow(color_matrix)), each = ncol(color_matrix))
    col_indices <- rep(seq_len(ncol(color_matrix)), times = nrow(color_matrix))

    # Flatten the color matrix into a vector
    colors <- as.vector(color_matrix)

    # Create a data frame
    df <- data.frame(row = row_indices,
                     col = col_indices,
                     color = colors)

    if (!is.null(downscale.image)){

      grid_df <- df
      grid_size <- as.numeric(downscale.image)

      # Create a new column to identify grid cells
      grid_df$grid_row <- floor(grid_df$row / grid_size)
      grid_df$grid_col <- floor(grid_df$col / grid_size)

      # Aggregate points by grid cell
      df <- grid_df %>%
        dplyr::group_by(grid_row, grid_col) %>%
        dplyr::summarize(
          row = mean(row),
          col = mean(col),
          color = dplyr::first(color),  # Or use a method to choose color
          .groups = "drop"
        )
    }
    plot <- plot %>% add_trace(df,
                               x = df$row / data@images[[show.image]]@scale.factors[[image.sf]],
                               y = df$col / data@images[[show.image]]@scale.factors[[image.sf]],
                               z = rep(0 - between.layer.height, times = dim(df)[1]),
                               type = "scatter3d",
                               mode = "markers",
                               name = "H&E",
                               marker = list(color = df$color))

  }
  return(plot)

}




#' Generates interactive 3D spatial density plot for m/z values
#'
#' This function generates an interactive HTML file that can be used to generate spatial density kernals based on the intensity of any selected m/z value.
#'
#' @param object SpaMTP Seurat object contains spatial metabolomics data
#' @param assay Character string defining the Seurat assay that contains intensity values to plot (default = "SPM").
#' @param slot Character string defining the assay slot to extract the intensity values from (default = "counts").
#' @param folder Character string defining the directory path to save the density map HTML file to (default = getwd()).
#' @param ... Arguments passed to SpaMTP::AnnotateSeuratMALDI()
#'
#' @return Return a html contains the annotation/average intensity of peakbins/density of distribution of each peak
#' @export
#'
#' @examples
#' # DensityMap(SpaMTP.obj)
DensityMap = function(object, assay = "SPM", slot = "counts", folder = getwd(),...){

  mass_matrix = Matrix::t(object[[assay]][slot])
  annotated_table = object[[assay]]@meta.data
  indices =  GetTissueCoordinates(object)[c("x", "y")]
  mzs =  annotated_table["raw_mz"]
  annotated_json  = '['
  for (i in 1:nrow(mzs)) {
    annotated_json =
      paste0(
        annotated_json,
        '[',
        jsonlite::toJSON(annotated_table$raw_mz[i]),
        ',',
        jsonlite::toJSON(annotated_table$mz_names[i]),
        ',',
        jsonlite::toJSON(gsub(
          "\\]",
          "_",
          gsub("\\[", "_", annotated_table$all_IsomerNames[i])
        )),
        ',',
        jsonlite::toJSON(annotated_table$all_Isomers[i]),
        '],',
        collapse = ""
      )
  }
  annotated_json  = paste0(annotated_json,
                           ']')


  # Histogram information, mz, mean_intensity, max_intensity

  histogram_data_to_be_added  = ""
  mean_intensity = Matrix::colSums(mass_matrix) / nrow(mass_matrix)
  for (i in 1:nrow(mzs)) {
    histogram_data_to_be_added = paste0(histogram_data_to_be_added,
                                        "{ mz:",
                                        mzs[i,],
                                        ", mean:",
                                        mean_intensity[i],
                                        "},")
  }

  #KDE curve
  data_points <-
    rep(as.numeric(sub("mz-", "", mzs$raw_mz)), as.integer(mean_intensity / 5))

  # Generate KDE
  kde <- density(data_points, bw = 0.05)
  kde$y <-  kde$y * max(mean_intensity) / max(kde$y)
  kde_json = paste0("[", jsonlite::toJSON(kde$x), ",", jsonlite::toJSON(kde$y), "]")

  # Initialize an empty character vector to store the elements - the density map required item
  elements <- character(ncol(mass_matrix))
  # Loop through each column of mass_matrix
  print("Parsing mass matrix information")
  pb = txtProgressBar(
    min = 0,
    max = ncol(mass_matrix),
    initial = 0,
    style  =  3
  )
  for (i in 1:ncol(mass_matrix)) {
    # Combine row and column indices
    element <- paste("[",
                     indices[which(mass_matrix[, i] != 0), 1],
                     ",",
                     indices[which(mass_matrix[, i] != 0), 2],
                     "]",
                     sep = "" ,
                     collapse = ",")

    # Store the element in the vector
    elements[i] <- paste("[", element, "]", sep = "")

    setTxtProgressBar(pb, i)
  }
  close(pb)
  # Combine all elements into a JSON array
  json_array <-
    paste("[", paste(elements, collapse = ","), "]", sep = "")

  # Correct double commas
  json_array <- sub("],]", "]]", json_array)

  # Combine all elements into a JSON array
  json_array <-
    paste("[", paste(elements, collapse = ","), "]", sep = "")

  # Correct double commas
  json_array <- sub("],]", "]]", json_array)

  #writeLines(json_array, paste0(file.path(folder),name,"/reconstuct_mz.json"))


  javascript = paste0(
    "
<!DOCTYPE html>
  <div id='Data'></div>
    <html lang='en'>

      <head>
      <meta charset='UTF-8'>
        <meta name='viewport' content='width=device-width, initial-scale=1.0'>
          <title>Interactive Bar and Dot Plot</title>
          <script src='https://cdn.plot.ly/plotly-latest.min.js'></script>
            <script src='https://cdnjs.cloudflare.com/ajax/libs/PapaParse/4.1.2/papaparse.min.js'></script>
              <script src='https://d3js.org/d3.v7.min.js'></script>
                  <script type='text/javascript' src='https://cdn.jsdelivr.net/gh/stdlib-js/stats-kde2d@umd/browser.js'></script>
                    <script type='text/javascript' src='https://cdn.jsdelivr.net/gh/stdlib-js/ndarray@umd/browser.js'></script>
                      </head>

                      <body>
                      <script type='text/javascript'>
                        (function () {
                          window.kde2d;
                        })();
                      </script>
                        <script type='text/javascript'>
                          (function () {
                            window.ndarray;
                          })();
                        </script>

    <style>
      #container {
          display: flex;
          flex-direction: column;
      }
      #plots {
          display: flex;
      }
      #barPlot, #secPlot {
          width: 100%;
          border: 0px solid black; /* Optional: to visualize the divs */
          box-sizing: border-box; /* Ensure borders are included in the dimensions */
          height: 450px;
          padding: 10px;
          margin-bottom: 0px;
          z-index: 1;
      }
      #buttonContainer {
        display: flex;
        flex-direction: column;
        position: absolute;
        left: 55%; /* Adjust to move horizontally within secPlot */
        top: 5px; /* Adjust to move vertically within secPlot */
      }
      #tablePlot {
          width: 100%;
          border: 2px solid black; /* Optional: to visualize the div */
          box-sizing: border-box; /* Ensure borders are included in the dimensions */
          height: 300px;
          padding: 30px;
      }
      #message {
        display: flex;
        flex-direction: column;
        position: absolute;
        left: 55%; /* Adjust to move horizontally within secPlot */
        top: 50px; /* Adjust to move vertically within secPlot */
      }
  </style>
      <div id='container'>
        <div id='plots'>
            <div id='barPlot'></div>
            <div id='secPlot'>
              <div id = 'warning'></div>
              <div id='buttonContainer'>
              <div id='message'>
                  <strong>Instructions:</strong> Entry a number to describe the number of equally spaced points at which the density is to be estimated.
                  <input type='number' id='numberEntry' min='0' max='",min(max(indices[,1],
                                                                               indices[,2])),"' z-index ='10' width = '30' value = '30'/>
              </div>
              </div>
            </div>
        </div>
        <div id='tablePlot'>
        </div>
    </div>
  <script>
    document.getElementById('numberEntry').addEventListener('input', function() {
    var value = parseInt(this.value);
    if (value > 100) {
        this.value = 100; // Cap the value at 100
    }
});
          const data = [",
    histogram_data_to_be_added,
    "];
          // Add event listener for unhover on bar plot
          var json_array =",
    json_array,
    "
          var annotated_table =",
    annotated_json,
    "
          var kde =",
    kde_json,
    "
          var trace = {
              x: kde[0],
              y: kde[1],
              mode: 'lines',
              name: 'Estimated average density spectrum',
              line: {shape: 'spline'},
              type: 'scatter',
              width: 0.3,
              hoverinfo: null
              //hoverinfo: 'Estimated density' + kde[1]
            };

            const xValues = data.map(d => d.mz);
            const yValues = data.map(d => d.mean);

            // Create trace for the bar plot
            var barTrace = {
              x: xValues,
              y: yValues,
              type: 'bar',
              width: 0.3,
              name: 'Peak bins',
              hoverinfo: 'm/z value' + xValues
            };

            // Create trace for the dot plot


            // Set layout for bar plot
            var barLayout = {
              hovermode: 'x unified',
              title: 'Bar Plot',
              xaxis: {
                title: 'm/z values'
              },
              yaxis: {
                title: 'intensity'
              },
              plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
              paper_bgcolor: 'rgba(0,0,0,0)',
              hoverdistance: 50
            };

            // Set layout for second plot

            // Create bar plot
            Plotly.newPlot('barPlot', [trace,barTrace], barLayout);

            function findArrayByNumber(array, number) {
              let searchString = 'mz-' + number;
              for (let i = 0; i < array.length; i++) {
                let element = array[i][0];
                if (element === searchString) {
                  return i;
                }
              }
              return -1;
            }
            document.querySelectorAll('.plotly .scatterlayer .trace').forEach(function(el, index) {
            if(index === 0) {
            el.style.pointerEvents = 'none';  // Disable interactivity for this trace
            }
            });
            // Add event listener for hover on bar plot
            var ind = 0
            function hideMessage() {
              warning.textContent = ''; // Clear the message text
              warning.style.display = 'none'; // Hide the message box
        }
            document.getElementById('barPlot').on('plotly_click', function (data_event){
              hideMessage()
              if(data_event.points.length == 2){
                var pointNumber = data_event.points[1].pointNumber;
              }else if (data_event.points[0].fullData.type == 'bar'){
                var pointNumber = data_event.points[0].pointNumber;
              }else{
                Plotly.purge('secPlot');
                warning.textContent = 'No explicit peak bin selected, consider closing the density spetrum and zoom in to select explicit m/z peak'; // Set the message text
                warning.style.display = 'block'; // Make the message box visible
              }

              ind = pointNumber
              console.log(pointNumber)
              var selected_mz = xValues[pointNumber]
              Plotly.purge('secPlot');
              displaydensity(pointNumber)
              var table_index = findArrayByNumber(annotated_table.map(coord => coord[1]), selected_mz);
              if(table_index != -1){
                displayTable(table_index);
              }else{
                var table_update = [[selected_mz],['Not annotated by given adduct/db'],
                                    ['Not annotated by given adduct/db'],
                                    ['Not annotated by given adduct/db']]
                var tableData = [
                  {
                    type: 'table',
                    header: {
                      values: ['mz', 'mz_name', 'annotated metabolites', 'entry of the metabolites'],
                      align: 'center',
                      line: {width: 1, color: 'black'},
                      fill: {color: 'grey'},
                      font: {family: 'Arial', size: 12, color: 'white'}
                    },
                    cells: {
                      values: table_update,
                      align: 'center',
                      line: {color: 'black', width: 1},
                      fill: {color: ['white', 'lightgrey']},
                      font: {family: 'Arial', size: 11, color: ['black']}
                    }
                  }
                ];
                var layout_table = {
                  title: 'Annotation table for m/z: ' + data[table_index].mz,
                  autosize: true,
                  plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
                  paper_bgcolor: 'rgba(0,0,0,0)',
                };
                Plotly.newPlot('tablePlot', tableData,layout_table);
              }
            });
            // Get the density plot done
            // If the element exists, hide the loading message
            var space  = 30
            document.getElementById('numberEntry').addEventListener('blur', function() {
            space = Number(event.target.value)
            displaydensity(ind)
            });
            function displaydensity(index) {
              var temp_x = json_array[index].map(coord => coord[0]);
              var temp_y = json_array[index].map(coord => coord[1]);
              var n = space;
              var shape = [n, 2];
              var strides = [1, n];
              var offset = 0;
              var order = 'column-major';

              var out = kde2d(temp_x, temp_y, {
                'n': n,
                'buffer': temp_x.concat(temp_y),
                'order': 'column-major',
                'offset': offset,
                'strides': strides
              });

              let twoDArray = [];
              for (let i = 0; i < out.z._buffer.length; i += n) {
                twoDArray.push(out.z._buffer.slice(i, i + n));
              }
              var arr = twoDArray;
              console.log(twoDArray)
              var data_u = [{
                z: twoDArray,
                type: 'surface'
              }];
              var layout_U = {
                title: 'kde2d plot for m/z: ' + data[index].mz,
                autosize: true,
                plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
                paper_bgcolor: 'rgba(0,0,0,0)',
              };
              //loadingMessageElement.innerText = '';
              // Create dot plot
              Plotly.newPlot('secPlot', data_u, layout_U);
            }
            displaydensity(0)



            //table plot
            function displayTable(index) {
              var rowData = annotated_table[index];
              var tableData = [
                {
                  type: 'table',
                  header: {
                    values: ['mz', 'mz_name', 'annotated metabolites', 'entry of the metabolites'],
                    align: 'center',
                    line: {width: 1, color: 'black'},
                    fill: {color: 'grey'},
                    font: {family: 'Arial', size: 12, color: 'white'}
                  },
                  cells: {
                    values: rowData,
                    align: 'center',
                    line: {color: 'black', width: 1},
                    fill: {color: ['white', 'lightgrey']},
                    font: {family: 'Arial', size: 11, color: ['black']}
                  }
                }
              ];
              var layout_table = {
                title: 'Annotation table for m/z: ' + data[index].mz,
                autosize: true,
                plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background for the plot area
                paper_bgcolor: 'rgba(0,0,0,0)',
              };

              Plotly.newPlot('tablePlot', tableData,layout_table);
            }
            displayTable(0);
            </script>
              </body>
              </html>"
  )

  if (!file.exists(folder)) {
    dir.create(paste0(file.path(folder), name))
    writeLines(javascript, paste0(file.path(folder), "/mzs_density_map.html"))
  } else{
    writeLines(javascript, paste0(file.path(folder), "/mzs_density_map.html"))
  }
}






#' Interactive Spatial Plot for visualising different m/z bin sizes
#'
#' This function implements a shinyApp for interactive viusalisation of the samples mean mass intensity spectra.
#' Users can adjust the resolution/bin size and observe changes in spatial expression respectively.
#'
#' @param obj SpaMTP object containing the intensity data to plot.
#' @param assay Character defining the Seurat object assay to extract the intensity data from (default = "Spatial").
#' @param slot Character defining the Seurat object assay slot to extract the intensity data from (default = "counts").
#' @param image Character defining the the name of the image to use for tissue coordinates and spatial plotting (default = "slice1").
#'
#' @import shiny
#' @importFrom shinyjs useShinyjs reset
#' @importFrom plotly plotlyOutput renderPlotly plot_ly add_bars add_lines layout
#' @importFrom ggplot2 theme_void
#'
#' @export
#'
#' @examples
#' #InteractiveSpatialPlot(spamtp)
InteractiveSpatialPlot <- function(obj, assay = "Spatial", slot = "counts", image = "slice1"){

  means <- Matrix::as.matrix(Matrix::rowMeans(obj[[assay]][slot]))
  df <- data.frame(means)
  df$x <- as.numeric(gsub("mz-", "", rownames(df)))
  df$y <- df$means
  rownames(df) <- NULL

  # Shiny app UI
  ui <- fluidPage(
    titlePanel("Interactive Spatial Intensity Plot"),
    sidebarLayout(
      sidebarPanel(
        numericInput("curve_center", "m/z Value:", value = 450, min = min(df$means), max = max(df$means), step = 1),
        radioButtons("bin_mode", "Binning Mode:",
                     choices = c("Absolute (m/z)" = "mz", "Relative (ppm)" = "ppm"),
                     selected = "mz"),
        fluidRow(
          column(6,
                 sliderInput("curve_width_slider", "Bin Width:", min = 0, max = 100, value = 10, step = 0.1)
          ),
          column(6,
                 numericInput("curve_width","", value = 10, min = 0, max = 1000)
          )
        ),
        numericInput("spot_size", "Plot Spot Size:", value = 10, min = 0, max = 100),  # Added spot size input

        fluidRow(
          column(6, numericInput("x_min", "X-Axis Minimum:", value = min(df$means), min = min(df$means), max = max(df$means), step = 1)),
          column(6, numericInput("x_max", "X-Axis Maximum:", value = max(df$means), min = min(df$means), max = max(df$means), step = 1))
        )
      ),
      mainPanel(
        plotlyOutput("plot"),
        verbatimTextOutput("selected_indices"),
        plotOutput("spatial_plot")
      )
    )
  )

  server <- function(input, output, session) {

    observeEvent(input$curve_width_slider, {
      updateNumericInput(session, "curve_width", value = input$curve_width_slider)
    })

    # Keep slider in sync with numeric input
    observeEvent(input$curve_width, {
      updateSliderInput(session, "curve_width_slider", value = input$curve_width)
    })


    output$plot <- renderPlotly({
      curve_center <- input$curve_center
      bin_mode <- input$bin_mode

      if (bin_mode == "ppm") {
        tolerance <- input$curve_width * 1e-6
        lower_bound <- curve_center - (curve_center * tolerance)
        upper_bound <- curve_center + (curve_center * tolerance)
        curve_width <- (upper_bound - lower_bound) / 2
      } else {
        curve_width <- input$curve_width
      }


      x_min <- input$x_min
      x_max <- input$x_max

      # Calculate normal distribution values based on 'means'
      curve_vals <- dnorm(df$means, mean = curve_center, sd = curve_width)

      # Scale curve to match max y value at curve_center
      max_y_at_center <- df$y[which.min(abs(df$means - curve_center))]

      # Generate new x-values for the curve
      new_x_vals <- seq(from = curve_center - curve_width, to = curve_center + curve_width, length.out = 100)
      new_curve_vals <- dnorm(new_x_vals, mean = curve_center, sd = curve_width)
      new_curve_vals <- (new_curve_vals - min(new_curve_vals)) / (max(new_curve_vals) - min(new_curve_vals))
      new_curve_vals <- new_curve_vals * max_y_at_center

      # Create histogram with curve overlay
      plot_ly() %>%
        add_bars(x = df$x, y = df$y, name = "Histogram") %>%
        add_lines(x = new_x_vals, y = new_curve_vals, name = "Curve", line = list(color = 'blue')) %>%
        layout(
          xaxis = list(title = "Means (m/z values)", range = c(x_min, x_max)),
          yaxis = list(title = "Frequency"),
          showlegend = TRUE
        )
    })

    selected_points <- reactive({
      curve_center <- input$curve_center
      curve_width <- input$curve_width
      bin_mode <- input$bin_mode

      if (bin_mode == "ppm") {
        tolerance <- curve_width * 1e-6
        lower_bound <- curve_center - (curve_center * tolerance)
        upper_bound <- curve_center + (curve_center * tolerance)
        curve_width <- (upper_bound - lower_bound) / 2
      } else {
        lower_bound <- curve_center - curve_width
        upper_bound <- curve_center + curve_width
      }

      df$x[df$x >= lower_bound & df$x <= upper_bound]
    })


    output$selected_indices <- renderText({
      indices <- selected_points()

      if (length(indices) > 0) {
        paste("Selected feature:", paste(paste0("mz-", indices), collapse = ", "))
      } else {
        "No rows under the curve."
      }
    })

    output$spatial_plot <- renderPlot({
      indices <- selected_points()
      spot_size <- input$spot_size

      if (length(indices) > 0) {

        curve_center <- input$curve_center
        curve_width <- input$curve_width
        bin_mode <- input$bin_mode

        if (bin_mode == "ppm") {
          tolerance <- curve_width * 1e-6
          lower_bound <- curve_center - (curve_center * tolerance)
          upper_bound <- curve_center + (curve_center * tolerance)
          curve_width <- (upper_bound - lower_bound) / 2
        }

        features_to_plot <- indices  # Example gene names based on index
        features_to_plot <- paste0("mz-",features_to_plot)
        feature_name <- sprintf("mz-%d \u00B1 %.1f", input$curve_center, input$curve_width)
        binned_obj <- BinMetabolites(obj, mzs = features_to_plot, assay = "Spatial", slot = "counts", bin_name = feature_name)
        Seurat::SpatialFeaturePlot(binned_obj, images = image, features = feature_name,  pt.size.factor = spot_size) &theme_void()
      } else {
        plot(1, type = "n", xlab = "", ylab = "", main = "No selected features")
        text(1, 1, "No data available")
      }
    })



  }

  shinyApp(ui, server)
}


########################################################################################################################################################################################################################


