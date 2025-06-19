#' Visualise Significant Pathways
#'
#' Displays the pathway analysis results form running the 'FishersPathwayAnalysis()' function
#'
#' @param SpaMTP SpaMTP Seurat object used to run FishersPathwayAnalysis function.
#' @param pathway_df Dataframe containing the pathway enrichment results (output from SpaMTP::FishersPathwayAnalysis function).
#' @param assay Character string defining the SpaMTP assay that contains m/z values (default = "SPM").
#' @param slot Character string defining the assay slot contain the intensity values (default = "counts").
#' @param min_n Integer value specifying the minimum number of analytes required to be present in a pathway (default = 3).
#' @param p_val_threshold The p-val cutoff to keep the pathways generated from fisher exact test (default = "0.1").
#' @param method Character string defining the statistical method used to calculate hclust (default = "ward.D2").
#' @param verbose Boolean indicating whether to show informative messages. If FALSE these messages will be suppressed (default = TRUE).
#' @param ... The arguments pass to stats::hclust
#'
#' @return A combined gg, ggplot object with pathway and dendrogram
#' @export
#'
#' @import grid
#'
#' @examples
#' #SpaMTP:::VisualisePathways(SpaMTP =seurat,pathway_df = pathway_df,p_val_threshold = 0.1,assay = "Spatial",slot = "counts")
VisualisePathways = function(SpaMTP,
                             pathway_df,
                             assay = "SPM",
                             slot = "counts",
                             min_n = 3,
                             p_val_threshold = 0.1,
                             method = "ward.D2",
                             verbose = TRUE,
                             ...) {
  if (is.null(min_n)){
    stop("Incorrect minimum analyte number! `min_n` must be set to a value > 1. Please adjust this value accordingly ...")
  }
  pathway_df = pathway_df[which(pathway_df$analytes_in_pathways>=min_n),]
  pathway_df$duplicate_pathways <- NA
  verbose_message(message_text = "Reducing synonymous pathways", verbose = verbose)
  index = c(1:nrow(pathway_df))
  merged_pathways = data.frame()
  pb = txtProgressBar(
    min = 0,
    max = length(index),
    initial = 0,
    style = 3
  )
  while (length(index) != 0) {
    pattern = stringr::str_extract(pathway_df$pathway_name[index[1]], pattern = "[A-Z][A-Z]\\([a-z0-9]")
    if (length(pattern) != 0 & !is.na(pattern)) {
      pattern = sub("\\(", "\\\\(", pattern)
      name = strsplit(pathway_df$pathway_name[index[1]], pattern)[[1]][1]
      full_name = paste0(name, pattern)
      frst_ind = which(grepl(
        pathway_df$pathway_name,
        pattern = full_name
      ))
      all_pathways = pathway_df[frst_ind, ]
      second_ind = which(duplicated(all_pathways$p_val))
      index = index[-which(index %in% frst_ind)]
      if (length(second_ind)>0){
        name = stringr::str_trim(name, side = "right")

        unique_pathways <- all_pathways[-second_ind, ]
        for (i in seq_along(unique_pathways$p_val)) {
          # Get indices of duplicates for each unique p_val
          duplicates <- which(all_pathways$p_val == unique_pathways$p_val[i])
          duplicate_ids <- all_pathways$pathway_id[duplicates]
          unique_pathways$duplicate_pathways[i] <- paste(duplicate_ids, collapse = ", ")
        }

        unique_pathways$pathway_name <- name
        merged_pathways = rbind(merged_pathways, unique_pathways)
      } else{
        merged_pathways = rbind(merged_pathways, all_pathways[frst_ind, ])
      }

    } else{
      merged_pathways = rbind(merged_pathways, pathway_df[index[1], ])
      index = index[-1]
    }
    setTxtProgressBar(pb, nrow(pathway_df) - length(index))
  }
  close(pb)
  merged_pathways = merged_pathways %>% filter(p_val <= p_val_threshold) %>% mutate(signif_at_005level =   ifelse(p_val <= 0.05, "Significant", "Non-significant"))

  retain_ind = 1:nrow(merged_pathways)

  gg_bar1 = with(
    merged_pathways,
    ggplot() +  geom_bar(
      aes(
        x = paste0(pathway_name, "(", pathway_id, ")"),
        y = total_in_pathways,
        color = "Total analytes in pathway"
      ),
      stat = "identity",
      fill = NA,
      position = "dodge"
    ) +
      geom_bar(
        aes(
          x = paste0(pathway_name, "(", pathway_id, ")"),
          y = analytes_in_pathways,
          colour = "Analytes detected in pathway"
        ),
        stat = "identity",
        fill = "blue",
        position = "dodge"
      ) +
      geom_point(
        aes(
          x = paste0(pathway_name, "(", pathway_id, ")"),
          y = 1.05 * max(total_in_pathways),
          size = p_val,
          fill = signif_at_005level
        ),
        color = "black",
        shape = 21,
        position = position_nudge(y = 0.5)
      ) +
      scale_size_area(max_size = 5, name = "p value") +
      scale_fill_manual(
        values = c(
          "Significant" = "green",
          "Non-significant" = "red"
        ),
        name = ""
      ) + scale_colour_manual(
        values = c(
          "Analytes detected in pathway" = "blue",
          "Total analytes in pathway" = "grey"
        ),
        name = ""
      ) +
      theme_minimal() + coord_flip() +
      theme(
        legend.position = "left",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )+labs(y= "Number of analytes in pathway", x = "Pathway names")
  )

  # Get the raster illustration for each pathway
  image_raster = list()
  verbose_message(message_text = "Generating images for visualising pathway enrichment across the sample ... ", verbose = verbose)
  pb = txtProgressBar(
    min = 0,
    max = nrow(merged_pathways),
    initial = 0,
    style = 3
  )
  mass_matrix = Matrix::t(SpaMTP[[assay]]@layers[[slot]])

  for(z in 1:nrow(merged_pathways)){
    mzs = paste0("mz-",
                 stringr::str_extract_all(merged_pathways$adduct_info[z], "\\d+\\.\\d+")[[1]])
    if ("mz-" %in% mzs){
      image_raster[[z]] <- NULL
    } else{
      mat_ind = which(row.names(SpaMTP[[assay]]@features) %in% mzs)
      pca_result <- prcomp(mass_matrix[, mat_ind])
      nc = 3
      pca_re_df <- pca_result[["x"]][, 1:nc]
      pca_df_normalized <- as.data.frame(apply(
        pca_re_df,
        MARGIN = 2 ,
        FUN =  function(x)
          (x - min(x)) / (max(x) - min(x))
      ))
      coords <- GetTissueCoordinates(SpaMTP)[c("x", "y")]
      # Convert the normalized UMAP result to an image matrix
      # 3 channels side by side

      rgb_m = array(dim = c(max(coords[, 1]), max(coords[, 2]), nc))
      for (j in 1:nc) {
        rgb_m[, , j] = matrix(pca_df_normalized[, j], nrow = max(coords[, 1]))
      }
      # Convert the image matrix to a raster object
      image_raster[[z]] <- as.raster(rgb_m)
      # # Plot the RGB image using ggplot2
      # ggplot() + annotation_custom(rasterGrob(image_raster, width = unit(1, "npc"), height = unit(1, "npc"))) +
      #   theme_void()
      setTxtProgressBar(pb, z)
    }
  }
  close(pb)
  if (length(image_raster)>0){
    for (k in 1:length(image_raster)) {
      if (!is.null(image_raster[[k]])){
        gg_bar1 =  gg_bar1 + ggplot2::annotation_custom(
          grid::rasterGrob(
            image_raster[[k]],
            width = unit(1, "npc"),
            height = unit(1, "npc")
          ),
          xmin = k - 0.5,
          xmax = k + 0.5,
          ymin = -3.5 ,
          ymax = -0.3
        )
      }
    }
  }
  gg_bar1 =  gg_bar1 + ylim(-2, max(merged_pathways$total_in_pathways) + 5)
  #gg_bar1

  # Adding the dendrogram

  # data(pathway, package = "SpaMTP")
  # data(analytehaspathway, package = "SpaMTP")
  jaccard_matrix = matrix(nrow = nrow(merged_pathways),
                          ncol = nrow(merged_pathways))

  for (i in 1:(nrow(merged_pathways) - 1)) {
    pathway_id_i = merged_pathways$pathway_id[i]
    pathway_content_i = unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_i)])])
    for (j in (i + 1):nrow(merged_pathways)) {
      pathway_id_j = merged_pathways$pathway_id[j]
      pathway_content_j = unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_j)])])

      jc_simi = length(intersect(pathway_content_i, pathway_content_j)) / length(union(pathway_content_i, pathway_content_j))
      jaccard_matrix[i, j] = jaccard_matrix[j, i] = jc_simi
    }
  }
  diag(jaccard_matrix) = 1


  # Generate a dendrogram
  hc <- as.dendrogram(hclust(as.dist(jaccard_matrix), method = method, ...))
  segment_hc <- with(ggdendro::segment(ggdendro::dendro_data(hc)),
                     data.frame(
                       x = y,
                       y = x,
                       xend = yend,
                       yend = xend
                     ))
  pos_table <- with(ggdendro::dendro_data(hc)$labels,
                    data.frame(
                      y_center = x,
                      gene = as.character(label),
                      height = 1
                    ))
  axis_limits <- with(pos_table, c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1)

  plt_dendr <- ggplot(segment_hc) +
    geom_segment(aes(
      x = sqrt(x),
      y = y,
      xend = sqrt(xend),
      yend = yend
    )) +
    scale_x_continuous(expand = c(0, 0.5),
                       limits = c(0,max(segment_hc$xend))) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = pos_table$y_center,
      labels = pos_table$gene,
      limits = axis_limits,
      position = "right"
    ) +
    labs(
      x = "Jacard distance",
      y = "",
      colour = "",
      size = ""
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

  # # Combine the dendrogram and bar plot
  # combined_plot <- grid.arrange(gg_bar1,ggplotify::as.ggplot(dendro), widths = c(5, 1), nrow =1)

  combined_plot = cowplot::plot_grid(gg_bar1,
                                     plt_dendr,
                                     align = 'h',
                                     rel_widths = c(6, 1))
  verbose_message(message_text = "\nDone", verbose = verbose)
  # Show the plot
  return(combined_plot)
}



#' Plot significantly enriched pathways per region
#'
#' Visualisation of Set Enrichment Analysis Results from `SpaMTP::FindRegionalPathways()`.
#'
#' @param regpathway A dataframe generated by the `SpaMTP::FindRegionalPathways()` function, containing identified regional pathways.
#' @param ident.column A character string specifying the column name of the `regpathway` dataframe containing the idents or clusters to compare (default = "Cluster_id").
#' @param selected_pathways A character vector specifying the names or IDs of pathways to be included in the analysis (e.g., `c("Amino acid metabolism", "WP1902", "Aspartate and asparagine metabolism")`). This argument is not case-sensitive (default = NULL).
#' @param sig_cutoff A numeric value defining the p-value cutoff for classifying significant pathways. If NULL the p-value cutoff will be = 0.05 (default = NULL).
#' @param num_display An integer specifying the number of pathways to display in the plot. If set to null will plot the smaller of either 10 pathways or the number of unique pathways provided in 'regpathway' (default = NULL).
#' @param text_size A numeric value controlling the size of the text elements in the plot. If NULL the default text size is 12 (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A `ggplot` object representing the set enrichment analysis results.
#' @export
#'
#' @examples
#' # PlotRegionalPathways(SpaMTP, ident = "clusters", regpathway = pathway_df)
PlotRegionalPathways <- function(regpathway,
                                 ident.column = "Cluster_id",
                                 selected_pathways = NULL,
                                 sig_cutoff = NULL,
                                 num_display = NULL,
                                 text_size = NULL,
                                 verbose = TRUE) {

  ## Checks for ident column
  if (!is.null(regpathway)){
    if (! "pathwayName" %in% colnames(regpathway)){
      stop("`pathwayName` column not found in `regpathway`! Please provide a data.frame containing the column named 'pathwayName', which stores the name of each pathway. `regpathway` can be generated by running `FindRegionalPathways()` ...")
    }
    if (!ident.column %in% colnames(regpathway)) {
      stop(
        "Column '",ident, "' not found in the provided `regpathway` dataframe! Make sure the ident.column is present ...")
    }
  } else {
    stop("Incorrect DE data.frame! A input for `regpathway` has not been provided. Please provide a data.frame containing significant pathways per ident. `regpathway` can be generated by running `FindRegionalPathways()` ...")
  }

  sig_pval <- sig_cutoff  %||% 0.05

  regpathway = regpathway %>% dplyr::group_by(pathwayName) %>% dplyr::mutate(
    Significance = ifelse(
      pval <= sig_pval,
      paste0("Significant at ",sig_pval," significance level"),
      "Not statistically significant"
    )
  ) %>% dplyr::mutate(group_importance = sum(abs(NES)))

  if (is.null(selected_pathways)) {
    regpathway = regpathway %>%  dplyr::filter(group_importance %in% sort(unique(regpathway$group_importance), decreasing = T)[1:(num_display %||% min(10, length(unique(
      regpathway$pathwayName
    ))))])
  } else{
    regpathway = regpathway %>% dplyr::filter((pathwayName %in% selected_pathways) |
                                                (sourceId %in% selected_pathways))
    regpathway = regpathway %>%  dplyr::filter(group_importance %in% sort(unique(regpathway$group_importance), decreasing = T)[1:(num_display %||% min(10, length(unique(
      regpathway$pathwayName
    ))))])
  }

  regpathway = regpathway %>% dplyr::mutate(pathnameid = paste0(pathwayName, "(", sourceId, ")"))
  pathwaynames = unique(regpathway$pathnameid)
  n = length(unique(regpathway$pathnameid))
  jaccard_matrix = matrix(nrow = n, ncol = n)
  colnames(jaccard_matrix) = rownames(jaccard_matrix) = pathwaynames
  verbose_message(message_text = "computing jaccard distance between pathways" , verbose = verbose)
  for (i in 1:(n - 1)) {
    pathway_id_i = sub(".*\\(([^)]+)\\).*", "\\1", pathwaynames[i])
    pathway_content_i = unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_i)])])
    for (j in (i + 1):n) {
      pathway_id_j = sub(".*\\(([^)]+)\\).*", "\\1", pathwaynames[j])
      pathway_content_j = unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_j)])])
      jc_simi = length(intersect(pathway_content_i, pathway_content_j)) / length(union(pathway_content_i, pathway_content_j))
      jaccard_matrix[i, j] = jaccard_matrix[j, i] = jc_simi
    }
  }
  diag(jaccard_matrix) = 1

  # Generate a dendrogram
  hc <- as.dendrogram(hclust(as.dist(jaccard_matrix)))
  segment_hc <- with(ggdendro::segment(ggdendro::dendro_data(hc)),
                     data.frame(
                       x = y,
                       y = x,
                       xend = yend,
                       yend = xend
                     ))
  pos_table <- with(ggdendro::dendro_data(hc)$labels,
                    data.frame(
                      y_center = x,
                      path = as.character(label),
                      height = 1
                    ))

  axis_limits <- with(pos_table, c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1)
  plt_dendr <- ggplot(segment_hc) +
    geom_segment(aes(
      x = sqrt(x),
      y = y,
      xend = sqrt(xend),
      yend = yend
    )) +
    scale_x_continuous(expand = c(0, 0.1), limits = c(0, sqrt(max(segment_hc$xend)))) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = pos_table$y_center,
      labels = stringr::str_wrap(pos_table$path, width = 40),
      limits = axis_limits,
      position = "right"
    ) +
    labs(
      x = "Jacard distance",
      y = "",
      colour = "",
      size = ""
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = ((text_size) %||% 12) - 3),
          axis.text.y = element_text(color= "black"))


  #### generate dot plot
  suppressWarnings({
    gg_dot = ggplot(data = regpathway, aes(
      x = factor(!!sym(ident.column), levels = sort(unique(!!sym(ident.column)))),
      y = factor(pathnameid, levels = ggdendro::dendro_data(hc)$labels$label)
    )) +
      geom_point(aes(colour = as.numeric(NES), size = as.numeric(sqrt(size)) * 5.1), stroke= 0) +
      scale_colour_gradient2(
        name = "Normalised enrichment score",
        low = "blue",
        # Original low color
        high = "red",
        # Original high color
        limits = c(min(regpathway$NES), max(regpathway$NES))
      ) +
      scale_size_continuous(name = "Number of altered metabolites \n in the pathway") +
      ggnewscale::new_scale_colour() +
      geom_point(shape = 1,
                 aes(colour = Significance, size = as.numeric(sqrt(size)) * 5.1 +
                       0.5, stroke = 1)) +
      scale_color_manual(
        values = setNames(
          c("black", "lightgrey"),  # Colors
          c(paste0("Significant at ", sig_pval, " significance level"), "Not statistically significant")  # Labels
        )
      ) +
      #labs(title = "Comparason of pathways expression between different cluster", y = "Pathways", x = "Clusters") +
      theme(
        title = element_text(size = text_size %||% 12, face = 'bold'),
        axis.text.x = element_text(
          size = text_size %||% 12,
          angle = 90,
          vjust = 0.5,
          hjust = 1,
        ),
        #20
        axis.title = element_text(size = text_size %||% 12),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = ((text_size) %||% 12) - 3),
        legend.text = element_text(size = ((text_size) %||% 12) - 3)
      ) +
      ggnewscale::new_scale_colour() + theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "lightgrey", size = 0.4)
      ) +   theme(
        legend.position = "left",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "horizontal"
      ) + xlab("Clusters")
  })
  u = gg_dot|plt_dendr
  return(u)
}






#' Plots the expression profile of a feature set corresponding to specified pathways onto a 2D scatter plot based on a dimensionality reduction technique.
#'
#' This function has been adapted from fgsea package. For more detail about function inputs please visit their [documentation](https://github.com/alserglab/fgsea/blob/master/R/geseca-plot.R#L266C1-L266C31).
#'
#' @param pathways Character vector of pathway names to plot. These pathways must be those stored in the `chempathway` object. To check call `unique(chempathway$pathwayName)`.
#' @param object SpaMTP Seurat object containing the dimensionality reduction to use for plotting.
#' @param title Optional title for the plot. If set to `NULL` the title will be the pathway name (default = NULL).
#' @param assay Character string defining the name of assay to use for data extraction. This `assay` should contain RAMP_IDs as feature names (default = `SeuratObject::DefaultAssay(object)`).
#' @param slot Character string defining the name of slot to extract data from (default = "scale.data").
#' @param reduction Character string defining which dimensionality reduction to use. If not specified, first searches for 'umap', then 'tsne', then 'pca' (default = NULL).
#' @param colors Vector of colors for the gradient (default = c("darkblue", "lightgrey", "darkred")).
#' @param guide Character string stating the type of legend to display (default = "colourbar").
#' @param ... Additional inputs taken by `Seurat::FeaturePlot()`. Check the relative documentation for more infomation.
#'
#' @return A ggplot object visualizing the score of each pathway across the relative reduction.
#' @export
#'
#' @examples
#' #PlotPathways(c("Glycolysis", "Acylcarnitine 3-Butenylcarnitine", "ABC transporters"), spamtp_obj, reduction = "umap"))
PlotPathways <- function(pathways, object, title=NULL,
                         assay=DefaultAssay(object),
                         slot = "scale.data",
                         reduction=NULL,
                         colors=c("darkblue", "lightgrey", "darkred"),
                         guide="colourbar",
                         ...) {

  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")
  pathway_db <- split(chempathway$rampId, chempathway$pathwayName)
  pathway_db <- pathway_db[!duplicated(tolower(names(pathway_db)))]

  pathway_db <- pathway_db[pathways]

  if (length(pathway_db) < 1) {
    stop("None of the pathway names provided were found in the pathway database. Please check the pathway names provided!")
  }

  if (is.null(title)) {
    titles <- names(pathway_db)
  } else if (length(title) != length(pathway_db)) {
    stop("Length of the specified titles does not match count of pathways")
  } else {
    titles <- title
  }

  # Use lapply to iterate over pathways and generate plots
  pathway_plots <- lapply(seq_along(pathway_db), function(i) {
    PlotSinglePathway(pathway=pathway_db[[i]], object=object,
                      title=titles[i], assay=assay, slot = slot,
                      reduction=reduction, colors=colors,
                      guide=guide, ...)
  })

  names(pathway_plots) <- pathways
  return(pathway_plots)
}









#' Plots the expression profile of a feature set corresponding to specified pathways onto a 2D scatter plot based on a dimensionality reduction technique.
#'
#' NOTE: This is a helper function for `PlotPathways`. This function has been adapted from fgsea package. For more detail about function inputs please visit their [documentation](https://github.com/alserglab/fgsea/blob/master/R/geseca-plot.R#L266C1-L266C31).
#'
#' @param pathway Character vector of pathway names to plot. These pathways must be those stored in the `chempathway` object. To check call `unique(chempathway$pathwayName)`.
#' @param object SpaMTP Seurat object containing the spatial data and coordinates for plotting.
#' @param title Optional title for the plot. If set to `NULL` the title will be the pathway name (default = NULL).
#' @param assay Character string defining the name of assay to use for data extraction. This `assay` should contain RAMP_IDs as feature names (default = `SeuratObject::DefaultAssay(object)`).
#' @param slot Character string defining the name of slot to extract data from (default = "scale.data").
#' @param reduction Character string defining which dimensionality reduction to use. If not specified, first searches for 'umap', then 'tsne', then 'pca' (default = NULL).
#' @param colors Vector of colors for the gradient (default = c("darkblue", "lightgrey", "darkred")).
#' @param guide Character string stating the type of legend to display (default = "colourbar").
#' @param ... Additional inputs taken by `Seurat::FeaturePlot()`. Check the relative documentation for more infomation.
#'
#' @return A ggplot object visualizing the pathway score spatially.
#' @export
#'
#' @examples
#' #glycolysis_list <- list("Glycolysis" = c("RAMP_C_000218730","RAMP_G_000012583","RAMP_G_000001171","RAMP_G_000007564","RAMP_C_000218226","RAMP_C_000040403","RAMP_C_000001115","RAMP_G_000008859"))
#' #PlotSinglePathway(glycolysis_list, spamtp_obj, reduction = "umap")
PlotSinglePathway <- function(pathway, object, title=NULL,
                              assay=DefaultAssay(object),
                              slot = "scale.data",
                              reduction=NULL,
                              colors=c("darkblue", "lightgrey", "darkred"),
                              guide="colourbar",
                              ...) {
  obj2 <- addGesecaScores(list(pathway=pathway), object, assay=assay, slot = slot, scale=TRUE)

  p <- Seurat::FeaturePlot(obj2, features = "pathway",
                           combine = FALSE, reduction=reduction, ...)

  p2 <- lapply(1:length(p), function(x){

    p2 <- p[[x]] + ggplot2::coord_fixed()
    p2$scales$scales[p2$scales$find("color")] <- NULL

    # Modify color scale
    suppressMessages(
      p2 <- p2 + ggplot2::scale_color_gradientn(limits=c(-3, 3), breaks=c(-3, 0, 3),
                                                colors=colors, oob=scales::squish,
                                                guide=guide, name="z-score"))

    if (!is.null(title)) {
      p2 <- p2 & ggplot2::ggtitle(title)
    }
    p2
  })


  if(length(p2)> 1){
    p2 <- cowplot::plot_grid(plotlist = p2, ncol = length(p2))
  } else{
    p2 <- p2[[1]]
  }

  return(p2)
}



#' Plots the spatial expression profile of a feature set corresponding to specified pathways
#'
#' This function has been adapted from fgsea package. For more detail about function inputs please visit their [documentation](https://github.com/alserglab/fgsea/blob/master/R/geseca-plot.R#L266C1-L266C31).
#'
#' @param pathways Character vector of pathway names to plot. These pathways must be those stored in the `chempathway` object. To check call `unique(chempathway$pathwayName)`.
#' @param object SpaMTP Seurat object containing the spatial data and coordinates for plotting.
#' @param images Character vector specifying which images to use for plotting from the SpaMTP Seurat Object.
#' @param title Optional title for the plot. If set to `NULL` the title will be the pathway name (default = NULL).
#' @param image.alpha Numeric value between 0 and 1 controlling the transparency of the underlying image (default = 1).
#' @param assay Character string defining the name of assay to use for data extraction. This `assay` should contain RAMP_IDs as feature names (default = `SeuratObject::DefaultAssay(object)`).
#' @param slot Character string defining the name of slot to extract data from (default = "scale.data").
#' @param colors Vector of colors for the gradient (default = c("darkblue", "lightgrey", "darkred")).
#' @param guide Character string stating the type of legend to display (default = "colourbar").
#' @param crop Boolean logical indicating whether to crop images (default = TRUE).
#' @param min.cutoff Numeric value defining the minimum cutoff value for color scale. If `NA` this value will be determined based on the minimum value of the dataset (default = NA).
#' @param max.cutoff Numeric value defining the maximum cutoff value for color scale. If `NA` this value will be determined based on the maximum value of the dataset (default = NA).
#' @param ncol Numeric value defining Number of columns for arranging plots. Note this is useful if multiple images have been supplied/containing within the one SpaMTP Seurat Object (default = NULL).
#' @param pt.size.factor Numeric value defining the spot point size for plotting (default = 1.6).
#' @param alpha Numeric vector controlling the transparency of points (default = c(1,1).
#' @param image.scale Character string defining the image resolution to use. This can be either "hires" or "lowres" (default = "lowres").
#' @param shape Integer value defining the shape of points. If shape = 21, circles will be plotted (default = 21).
#' @param stroke Numeric value stating the stroke width for points (default = NA).
#' @param interactive Boolean logical indicating whether to create interactive plots. NOTE: this functionality is implemented through `Seurat::SpatialFeaturePlot` (default = FALSE).
#' @param information An optional data.frame or matrix of extra information to be displayed on hover. NOTE: this functionality is implemented through `Seurat::SpatialFeaturePlot` (default = NULL).
#' @param image.labels Character vector specifying optional labels for multiple images (default = NULL).
#'
#' @return A ggplot object visualizing the pathway score spatially.
#' @export
#'
#' @examples
#' #PlotPathwaysSpatially(c("Glycolysis", "Acylcarnitine 3-Butenylcarnitine", "ABC transporters"), spamtp_obj)
PlotPathwaysSpatially <- function(pathways, object, images, title=NULL,image.alpha = 1,
                                  assay=SeuratObject::DefaultAssay(object), slot="scale.data",
                                  colors=c("darkblue", "lightgrey", "darkred"),
                                  guide="colourbar",
                                  crop = TRUE,
                                  min.cutoff = NA,
                                  max.cutoff = NA,
                                  ncol = NULL,
                                  pt.size.factor = 1.6,
                                  alpha = c(1, 1),
                                  image.scale = "lowres",
                                  shape = 21,
                                  stroke = NA,
                                  interactive = FALSE,
                                  information = NULL,
                                  image.labels = NULL
) {

  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")
  pathway_db <- split(chempathway$rampId, chempathway$pathwayName)
  pathway_db <- pathway_db[!duplicated(tolower(names(pathway_db)))]

  pathway_db <- pathway_db[pathways]

  if (length(pathway_db) < 1) {
    stop("None of the pathway names provided were found in the pathway database. Please check the pathway names provided!")
  }

  if (is.null(title)) {
    titles <- names(pathway_db)
  } else if (length(title) != length(pathway_db)) {
    stop("Length of the specified titles does not match count of pathways")
  } else {
    titles <- title
  }

  # Use lapply to iterate over pathways and generate plots
  ps <- lapply(seq_along(pathway_db), function(i) {
    PlotSinglePathwaySpatially(pathway_db[[i]], object, images, image.alpha = image.alpha,
                               title=titles[i], assay=assay, slot=slot,
                               colors=colors, guide=guide, crop = crop,
                               ncol = ncol,
                               min.cutoff = min.cutoff,
                               max.cutoff = max.cutoff,

                               pt.size.factor = pt.size.factor,
                               alpha = alpha,
                               image.scale = image.scale,
                               shape = shape,
                               stroke = stroke,
                               interactive = interactive,
                               information = information,
                               image.labels = image.labels)
  })

  names(ps) <- pathways
  return(ps)
}




#' Plot expression profile of a single RAMP pathway spatially
#'
#' NOTE: This is a helper function for `PlotPathwaysSpatially`. This function has been adapted from fgsea package. For more detail about function inputs please visit their [documentation](https://github.com/alserglab/fgsea/blob/master/R/geseca-plot.R#L266C1-L266C31).
#'
#' @param pathway Character vector of analyte IDs (genes or metabolites) in the pathway.
#' @param object SpaMTP Seurat object containing the spatial data and coordinates for plotting.
#' @param images Character vector specifying which images to use for plotting from the SpaMTP Seurat Object.
#' @param title Optional title for the plot. If set to `NULL` the title will be the pathway name (default = NULL).
#' @param image.alpha Numeric value between 0 and 1 controlling the transparency of the underlying image (default = 1).
#' @param assay Character string defining the name of assay to use for data extraction. This `assay` should contain RAMP_IDs as feature names (default = `SeuratObject::DefaultAssay(object)`).
#' @param slot Character string defining the name of slot to extract data from (default = "scale.data").
#' @param colors Vector of colors for the gradient (default = c("darkblue", "lightgrey", "darkred")).
#' @param guide Character string stating the type of legend to display (default = "colourbar").
#' @param crop Boolean logical indicating whether to crop images (default = TRUE).
#' @param min.cutoff Numeric value defining the minimum cutoff value for color scale. If `NA` this value will be determined based on the minimum value of the dataset (default = NA).
#' @param max.cutoff Numeric value defining the maximum cutoff value for color scale. If `NA` this value will be determined based on the maximum value of the dataset (default = NA).
#' @param ncol Numeric value defining Number of columns for arranging plots. Note this is useful if multiple images have been supplied/containing within the one SpaMTP Seurat Object (default = NULL).
#' @param pt.size.factor Numeric value defining the spot point size for plotting (default = 1.6).
#' @param alpha Numeric vector controlling the transparency of points (default = c(1,1).
#' @param image.scale Character string defining the image resolution to use. This can be either "hires" or "lowres" (default = "lowres").
#' @param shape Integer value defining the shape of points. If shape = 21, circles will be plotted (default = 21).
#' @param stroke Numeric value stating the stroke width for points (default = NA).
#' @param interactive Boolean logical indicating whether to create interactive plots. NOTE: this functionality is implemented through `Seurat::SpatialFeaturePlot` (default = FALSE).
#' @param information An optional data.frame or matrix of extra information to be displayed on hover. NOTE: this functionality is implemented through `Seurat::SpatialFeaturePlot` (default = NULL).
#' @param image.labels Character vector specifying optional labels for multiple images (default = NULL).
#'
#' @return A ggplot object visualizing the pathway score spatially.
#' @export
#'
#' @examples
#' #glycolysis_list <- list("Glycolysis" = c("RAMP_C_000218730","RAMP_G_000012583","RAMP_G_000001171","RAMP_G_000007564","RAMP_C_000218226","RAMP_C_000040403","RAMP_C_000001115","RAMP_G_000008859"))
#' #PlotSinglePathwaySpatially(glycolysis_list, spamtp_obj)
PlotSinglePathwaySpatially <- function(pathway, object, images, title=NULL, image.alpha= 1,
                                       assay=SeuratObject::DefaultAssay(object), slot="scale.data",
                                       colors=c("darkblue", "lightgrey", "darkred"),
                                       guide="colourbar",
                                       crop = TRUE,
                                       min.cutoff = NA,
                                       max.cutoff = NA,
                                       ncol = NULL,
                                       pt.size.factor = 1.6,
                                       alpha = c(1, 1),
                                       image.scale = "lowres",
                                       shape = 21,
                                       stroke = NA,
                                       interactive = FALSE,
                                       information = NULL,
                                       image.labels = NULL
) {

  obj2 <- addGesecaScores(list(pathway_x=pathway), object, assay=assay, slot=slot, scale=TRUE)


  if (class(obj2@images[[images]]) == "FOV"){
    p <- Seurat::ImageFeaturePlot(
      obj2,
      features="pathway_x",
      fov=images,
      boundaries = NULL,
      cols = colors[2:3],
      size = pt.size.factor,
      min.cutoff = min.cutoff,
      max.cutoff = max.cutoff,
      alpha = alpha[1],
      dark.background = "white",
      crop = crop
    )

  }else{
    p <- Seurat::SpatialFeaturePlot(obj2, features="pathway_x", images=images,
                                    combine=FALSE,
                                    image.alpha=image.alpha,
                                    crop = crop,
                                    min.cutoff = min.cutoff,
                                    max.cutoff = max.cutoff,
                                    pt.size.factor = pt.size.factor,
                                    alpha = alpha,
                                    image.scale = image.scale,
                                    shape = shape,
                                    stroke = stroke,
                                    interactive = interactive,
                                    information = information)

  }

  if(length(p)> 1){
    if(!is.null(image.labels)){
      title <- paste0(title, " (", image.labels, ")")
    } else{
      title <- rep(title, times = length(p))
    }
  }

  p2 <- lapply(1:length(p), function(x){

    p[[x]]$scales$scales[p[[x]]$scales$find("fill")] <- NULL

    suppressMessages({
      p2 <- p[[x]] + ggplot2::scale_fill_gradientn(limits=c(-3, 3), breaks=c(-3, 0, 3),
                                                   oob=scales::squish, colors=colors,
                                                   guide=guide, name="z-score") +
        ggplot2::theme(legend.position = ggplot2::theme_get()$legend.position)
    })

    if (!is.null(title)) {
      p2 <- p2 & ggplot2::ggtitle(title[x])
    }
    p2
  })



  if(length(p2)> 1){
    if(!is.null(ncol)){
      p2 <- cowplot::plot_grid(plotlist = p2, ncol = ncol)
    } else {
      p2 <- cowplot::plot_grid(plotlist = p2, ncol = length(images))
    }
  } else{
    p2 <- p2[[1]]
  }


  #p$scales$scales[p$scales$find("fill")] <- NULL

  #suppressMessages({
  #    p2 <- p + ggplot2::scale_fill_gradientn(limits=c(-3, 3), breaks=c(-3, 0, 3),
  #                                   oob=scales::squish, colors=colors,
  #                                   guide=guide, name="z-score") +
  #        ggplot2::theme(legend.position = ggplot2::theme_get()$legend.position)
  #})

  #if (!is.null(title)) {
  #    p2 <- p2 #+ ggplot2::ggtitle(title)
  #}
  return(p2)
}




#' Add GESECA Scores to SpaMTP Object
#'
#' Adds Feature Set Enrichment Score for a specific pathway to a SpaMTP Seurat object's metadata. This is a helper function for `PlotPathwaysSpatially`.
#'
#' @param pathways List of gene/metabolite sets where the list names (pathway name) become metadata column names.
#' @param object SpaMTP Seurat object to add the pathway scores to.
#' @param assay Character string defining the name of assay to use for data extraction (default  = `SeuratObject::DefaultAssay(object)).
#' @param slot Character string stating the name of slot to extract data from (default = "scale.data").
#' @param prefix Character string to append before each pathway name in metadata columns (default = "").
#' @param scale Boolean logical value indicating whether to scale the GESECA score (default = FALSE).
#'
#' @return SpaMTP Seurat object with added GESECA scores in metadata.
#'
#' @examples
#' ### HELPER FUNCTION
#' #glycolysis_list <- list("Glycolysis" = c("RAMP_C_000218730","RAMP_G_000012583","RAMP_G_000001171","RAMP_G_000007564","RAMP_C_000218226","RAMP_C_000040403","RAMP_C_000001115","RAMP_G_000008859"))
#' #spamtp_obj <- addGesecaScores(pathways = glycolysis_list, object = spamtp_obj ,assay = "RNA", scale = TRUE)
addGesecaScores <- function(pathways,
                            object,
                            assay=SeuratObject::DefaultAssay(object),
                            slot = "scale.data",
                            prefix="",
                            scale=FALSE) {
  x <- Seurat::GetAssay(object, assay)
  E <- x[slot]

  res <- object


  for (i in seq_along(pathways)) {
    pathway <- pathways[[i]]
    pathway <- intersect(unique(pathway), rownames(E))
    score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
    score <- scale(score, center=TRUE, scale=scale)
    res@meta.data[[paste0(prefix, names(pathways)[i])]] <- score
  }


  return(res)
}








