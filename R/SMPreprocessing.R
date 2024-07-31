library(Seurat)
library(ggplot2)
library(ggridges)
library(Matrix)
library(edgeR)
library(stats)
library(dplyr)
library(tidyr)


#' Normalizes m/z intensity data stored in a Seurat Object
#'
#' @param data Seurat Object to be normalized.
#' @param normalisation.type Character string defining the normalization method to run. Options are either c("TIC", "LogNormalize", "RC") which represent Total Ion Current (TIC) normalization, Log Normalization or counts per million (RC), respectively (default = "TIC").
#' @param scale.factor Numeric value that sets the scale factor for pixel/spot level normalization. Following normalization the total intensity value across each pixel will equal this value. If scale.factor = NULL, TIC normalization will use a scale factor = number of m/z and Log Normalisation will use a scale factor = 10000 (default = NULL).
#' @param assay Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").
#' @param slot  Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A Seurat Object with intensity values normalized. Normalized data is stored in the $data slot of the specified assay
#' @export
#'
#' @examples
#' # normalised_data <- NormalizeSMData(SeuratObject)
NormalizeSMData <- function(data, normalisation.type = 'TIC', scale.factor = NULL, assay = "Spatial", slot = "counts", verbose = TRUE) {

  if (is.null(normalisation.type)) {
    stop("Error: no normalisation.type is select. Please enter either 'LogNormalize' or 'TIC'")
  }

  if (!(normalisation.type == "TIC" | normalisation.type == "LogNormalize"| normalisation.type == "RC")) {
    stop("Error: incorrect normalisation.type is select. Please enter either 'LogNormalize', 'TIC' or 'RC'")
  }

  if (is.null(scale.factor)) {
    scale.factor <- 10000
  }


  if (normalisation.type == 'TIC') {
    scale.factor <- length(rownames(data))
    normalisation.type <- "RC"
  }

  normalised.data <- Seurat::NormalizeData(data, normalization.method = normalisation.type, scale.factor = scale.factor, verbose = verbose)

  return(normalised.data)

}

#' Performs TMM normalization between categories based on a specified ident
#'
#' @param combined.obj Seurat object that contains groups being normalized.
#' @param ident Character string defining the column name or ident group to normalize between.
#' @param refIdent Character string specifying one class/group type to use as a reference for TMM normalisation.
#' @param normalisation.type Character string defining the normalization method to run. Options are either c("CPM", "TIC", "LogNormalize") which represent counts per million (CPM), Total Ion Current (TIC) normalization or Log Normalization, respectively (default = "CPM").
#' @param CPM.scale.factor Numeric value that sets the scale factor for pixel/spot level normalization. Following normalization the total intensity value across each pixel will equal this value (default = 1e6).
#' @param assay Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").
#' @param slot  Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).
#'
#' @return Seurat object with count values normalised and corrected for between categories
#' @export
#'
#' @examples
#' # norm.data <- TMMNormalize(SeuratObj, ident = "samples", refIdent = "sample1", normalisation.type = "CPM")
TMMNormalize <- function(combined.obj, ident, refIdent, normalisation.type = "CPM", CPM.scale.factor = 1e6, assay = "Spatial", slot = "counts", verbose = FALSE) {

  data_list <- list()
  Seurat::Idents(combined.obj) <- ident
  if (length(unique(Seurat::Idents(combined.obj))) <= 1){
    stop("Specified ident has 0 or 1 catagory in seurat object. Length of Idents must be > 1 for TMM (between sample) normalisation factors to be calculated")

  }


  if (!(refIdent %in% unique(Seurat::Idents(combined.obj)))){
    stop("The refIdent supplied is not present in the ident column. Please specify a group that is found within the ident column specififed")
  }


  if (normalisation.type == "CPM") {
    normalisation.type <- "RC"
  }

  if (normalisation.type == "TIC") {
    normalisation.type <- "RC"
    CPM.scale.factor <- length(rownames(combined.obj))
  }


  for (name in unique(Idents(combined.obj))){

        sub <- subset_SPM(combined.obj, idents = name, verbose = verbose)

    data_list[[name]] <- sub
  }

  df <- data.frame(mz = rownames(data_list[[1]]))
  rownames(df) <- df$mz

  for (dataset in names(data_list)){
    rowsum <- Matrix::rowSums(data_list[[dataset]][[assay]][slot])
    df[dataset] <- rowsum
  }
  mtx <- Matrix::as.matrix(df[names(data_list)])

  factors <- edgeR::calcNormFactors(mtx, method = "TMM", refColumn = refIdent)

  norm_data_list <- list()
  for (name in names(data_list)){
    norm.data <- NormalizeSMData(data_list[[name]], normalisation.type = normalisation.type, scale.factor = (CPM.scale.factor / factors[[name]]), assay = assay, slot = slot)
    norm_data_list[[name]] <- norm.data
  }

  merged.data <- SeuratObject::JoinLayers(merge(norm_data_list[[1]], y = norm_data_list[2: length(names(norm_data_list))]), merge.data = TRUE)
return(merged.data)
}



########### Pre-processing plots #############


#' Helper function for QC plots by generating intensity count data
#'
#' @param seurat.obj Seruat object containing the intensity data.
#' @param group.by Character string specifying the meta.data column to group by (default = NULL).
#' @param assay Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").
#' @param slot  Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").
#' @param bottom.cutoff Numeric value defining the percent of data to exclude for the lower end of the distribution. A bottom.cutoff = 0.05 will remove the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for the upper end of the distribution. A top.cutoff = 0.05 will remove the top 5% of data point (default = NULL).
#' @param log.data Boolean value indicating whether to log transform the y-axis values (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).

#'
#' @return A data.frame containing the relative transformed and sum counts required for various QC plots
#' @export
#'
#' @examples
#' # df <- statPlot(SeuratObj, group.by = "sample", bottom.cutoff = 0.05, top.cutoff = 0.05, log.data = TRUE)
statPlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", bottom.cutoff = NULL, top.cutoff = NULL, log.data = FALSE, verbose = FALSE){

  data_list <- list()
  if (!(is.null(group.by))){
    Seurat::Idents(seurat.obj) <- group.by
    for (ident in unique(Seurat::Idents(seurat.obj))){
          sub <- subset_SPM(seurat.obj, idents = ident, verbose = verbose)

      data_list[[ident]] <- sub
    }
  } else {
    data_list[["data"]] <- seurat.obj
  }

  df <- data.frame(mz = rownames(data_list[[1]]))
  rownames(df) <- df$mz

  for (dataset in names(data_list)){
    rowsum <- Matrix::rowSums(data_list[[dataset]][[assay]][slot])
    df[dataset] <- rowsum
  }

  df2 <- tidyr::pivot_longer(df, cols =  names(data_list), names_to = "var", values_to = "x")

  df2 <- df2 %>% dplyr::arrange(x)

  if (!(is.null(bottom.cutoff))){
    df2 <- df2 %>%
      dplyr::group_by(var) %>%
      dplyr::mutate(bottom_cutoff = stats::quantile(x, bottom.cutoff))

    df2 <- df2 %>%  dplyr::group_by(var) %>%dplyr::filter(x >= bottom_cutoff)

    verbose_message(message_text = paste0("Removing bottom ", bottom.cutoff*100, "% of datapoints"), verbose = verbose)
  }

  if (!(is.null(top.cutoff))){
    df2 <- df2 %>%
      dplyr::group_by(var) %>%
      dplyr::mutate(top_cutoff = quantile(x, 1- top.cutoff))

    df2 <- df2 %>%  dplyr::group_by(var) %>% dplyr::filter(x <= top_cutoff)

    verbose_message(message_text = paste0("Removing top ", top.cutoff*100, "% of datapoints"), verbose = verbose)

  }

  if (log.data) {
    df2$x <- log10(df2$x + 1)
  }

  return(df2)
}






#' Generates a ridge plot of Spatial Metabolomic intensity data
#'
#' @param seurat.obj Seurat object containing the metabolomic intensity data.
#' @param group.by Character string specifying the meta.data column to group by (default = NULL).
#' @param assay Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").
#' @param slot  Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").
#' @param title Character string of the plot title (default = "RidgePlot").
#' @param x.lab Character string of the x-axis label (default = "var").
#' @param y.lab Character string of the y-axis label (default = "intensity").
#' @param bottom.cutoff Numeric value defining the percent of data to exclude for the lower end of the distribution. A bottom.cutoff = 0.05 will remove the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for the upper end of the distribution. A top.cutoff = 0.05 will remove the top 5% of data point (default = NULL).
#' @param bins number of bins to group
#' @param log.data Boolean value indicating whether to log transform the y-axis values (default = FALSE).
#' @param cols Vector of strings defining the colours to use for plotting. This vector should match the length of unique groups (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).
#'
#' @export
#'
#' @examples
#' # MZRidgePlot(SeuratObj, group.by = "sample")
MZRidgePlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", title = "RidgePlot", x.lab = "var", y.lab = "intensity", bottom.cutoff = NULL, top.cutoff = NULL, bins = 1000,log.data = FALSE, cols = NULL, verbose = FALSE){
  data <- statPlot(seurat.obj = seurat.obj,
                   group.by = group.by,
                   assay = assay,
                   slot = slot,
                   bottom.cutoff = bottom.cutoff,
                   top.cutoff = top.cutoff,
                   log.data = log.data,
                   verbose = verbose)

  ridge_plot <- ggplot2::ggplot(data, ggplot2::aes(y=var, x=x,  fill=var)) +
    ggridges::geom_density_ridges(alpha=0.6, stat="binline", bins=bins) +
    ggridges::theme_ridges() +  ggplot2::labs(title = title, x = x.lab, y = y.lab)

  if (!(is.null(cols))){
    ridge_plot <- ridge_plot + ggplot2::scale_fill_manual(values = cols)
  }

  return(ridge_plot)
}



#' Generates a Violin plot of Spatial Metabolomic intensity data
#'
#' @param seurat.obj Seurat object containing the metabolomic intensity data.
#' @param group.by Character string specifying the meta.data column to group by (default = NULL).
#' @param assay Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").
#' @param slot  Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").
#' @param title Character string of the plot title (default = "VlnPlot").
#' @param x.lab Character string of the x-axis label (default = "var").
#' @param y.lab Character string of the y-axis label (default = "intensity").
#' @param show.points Boolean value describing whether to show each individual data point (default = TRUE).
#' @param bottom.cutoff Numeric value defining the percent of data to exclude for the lower end of the distribution. A bottom.cutoff = 0.05 will remove the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for the upper end of the distribution. A top.cutoff = 0.05 will remove the top 5% of data point (default = NULL).
#' @param log.data Boolean value indicating whether to log transform the y-axis values (default = FALSE).
#' @param cols Vector of strings defining the colours to use for plotting. This vector should match the length of unique groups (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).
#'
#' @export
#'
#' @examples
#' # MZVlnPlot(SeuratObj, group.by = "sample",  bottom.cutoff = 0.05)
MZVlnPlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", title = "VlnPlot", x.lab = "var", y.lab = "intensity", show.points = TRUE, bottom.cutoff = NULL, top.cutoff = NULL,log.data = FALSE, cols = NULL, verbose = FALSE){

  data <- statPlot(seurat.obj = seurat.obj,
                   group.by = group.by,
                   assay = assay,
                   slot = slot,
                   bottom.cutoff = bottom.cutoff,
                   top.cutoff = top.cutoff,
                   log.data = log.data,
                   verbose = verbose)

  violin_plot <- ggplot2::ggplot(data, ggplot2::aes(x = var , y = x, fill = var)) +
    ggplot2::geom_violin() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title, x = x.lab, y = y.lab)

  if (show.points){
    violin_plot <- violin_plot + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.2), alpha = 0.5)
  }

  if (!(is.null(cols))){
    violin_plot <- violin_plot + ggplot2::scale_fill_manual(values = cols)
  }

  return(violin_plot)
}


#' Generates a Boxplot of Spatial Metabolomic intensity data
#'
#' @param seurat.obj Seurat object containing the metabolomic intensity data.
#' @param group.by Character string specifying the meta.data column to group by (default = NULL).
#' @param assay Character string defining the name of the Seurat Object assay to pull the corresponding intensity data from (default = "Spatial").
#' @param slot  Character string defining the name of the slot within the Seurat Object assay to pull the corresponding intensity data from (default = "counts").
#' @param title Character string of the plot title (default = "BoxPlot").
#' @param x.lab Character string of the x-axis label (default = "var").
#' @param y.lab Character string of the y-axis label (default = "intensity").
#' @param show.points Boolean value describing whether to show each individual data point (default = TRUE).
#' @param bottom.cutoff Numeric value defining the percent of data to exclude for the lower end of the distribution. A bottom.cutoff = 0.05 will remove the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for the upper end of the distribution. A top.cutoff = 0.05 will remove the top 5% of data point (default = NULL).
#' @param log.data Boolean value indicating whether to log transform the y-axis values (default = FALSE).
#' @param cols Vector of strings defining the colours to use for plotting. This vector should match the length of unique groups (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = FALSE).
#'
#' @export
#'
#' @examples
#' # MZBoxPlot(SeuratObj, group.by = "sample",  bottom.cutoff = 0.05)
MZBoxPlot <- function (seurat.obj, group.by = NULL, assay = "Spatial", slot = "counts", title = "BoxPlot", x.lab = "var", y.lab = "intensity", show.points = TRUE, bottom.cutoff = NULL, top.cutoff = NULL,log.data = FALSE, cols = NULL, verbose = FALSE){
  data <- statPlot(seurat.obj = seurat.obj,
                   group.by = group.by,
                   assay = assay,
                   slot = slot,
                   bottom.cutoff = bottom.cutoff,
                   top.cutoff = top.cutoff,
                   log.data = log.data,
                   verbose = verbose)

  box_plot <- ggplot2::ggplot(data, ggplot2::aes(x = var , y = x, fill = var)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title, x = x.lab, y = y.lab)

  if (show.points){
    box_plot <- box_plot + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.2), alpha = 0.5)
  }

  if (!(is.null(cols))){
    box_plot <- box_plot + ggplot2::scale_fill_manual(values = cols)
  }

  return(box_plot)
}













