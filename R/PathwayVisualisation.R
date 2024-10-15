#' Displays the pathway analysis results form running the 'FishersPathwayAnalysis()' function
#'
#' @param SpaMTP SpaMTP Seurat object used to run FishersPathwayAnalysis function.
#' @param pathway_df Dataframe containing the pathway enrichment results (output from SpaMTP::FishersPathwayAnalysis function).
#' @param assay Character string defining the SpaMTP assay that contains m/z values (default = "SPM").
#' @param slot Character string defining the assay slot contatin ght intesity values (default = "counts").
#' @param p_val_threshold The p-val cuttoff to keep the pathways generated from fisher exact test (default = "0.1").
#' @param method Character string defining the statistical method used to calculate hclust (default = "ward.D2").
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
                             p_val_threshold = 0.1,
                             method = "ward.D2",
                             verbose = TRUE,
                             ...) {
  pathway_df = pathway_df[which(pathway_df$analytes_in_pathways>=3),]
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
      frst_ind = which(grepl(
        pathway_df$pathway_name,
        pattern = sub("\\(", "\\\\(", pattern)
      ))
      all_pathways = pathway_df[frst_ind, ]
      second_ind = which(duplicated(
        paste0(
          all_pathways$gene_name_list,
          all_pathways$met_name_list
        )
      ))
      index = index[-which(index %in% frst_ind)]
      merged_pathways = rbind(merged_pathways, all_pathways[-second_ind, ])
    } else{
      merged_pathways = rbind(merged_pathways, pathway_df[index[1], ])
      index = index[-1]
    }
    setTxtProgressBar(pb, nrow(pathway_df) - length(index))
  }
  close(pb)
  merged_pathways = merged_pathways %>% filter(p_val <= p_val_threshold) %>% mutate(signif_at_005level =   ifelse(p_val <= 0.05, "Significant", "Non-significant"))

  retain_ind = 1:nrow(merged_pathways)
  for(z in 1:nrow(merged_pathways)){
    mzs = paste0("mz-",
                 stringr::str_extract_all(merged_pathways$adduct_info[z], "\\d+\\.\\d+")[[1]])
    mat_ind = which(row.names(SpaMTP[[assay]]@features) %in% mzs)
    if(length(mat_ind)<=2){
      retain_ind =   retain_ind[-which(retain_ind == z)]
    }
  }
  merged_pathways  =  merged_pathways[retain_ind,]
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
  verbose_message(message_text = "\nRunning PCA for dimension reduction", verbose = verbose)
  pb = txtProgressBar(
    min = 0,
    max = nrow(merged_pathways),
    initial = 0,
    style = 3
  )
  mass_matrix = Matrix::t(SpaMTP[[assay]]@layers[[slot]])



  for (i in 1:nrow(merged_pathways)) {
    mzs = paste0("mz-",
                 stringr::str_extract_all(merged_pathways$adduct_info[i], "\\d+\\.\\d+")[[1]])
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
    image_raster[[i]] <- as.raster(rgb_m)
    # # Plot the RGB image using ggplot2
    # ggplot() + annotation_custom(rasterGrob(image_raster, width = unit(1, "npc"), height = unit(1, "npc"))) +
    #   theme_void()
    setTxtProgressBar(pb, i)
  }
  close(pb)
  for (k in 1:length(image_raster)) {
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
  # dendro <- ggtree(as.phylo(hc), layout = "rectangular")+scale_x_reverse()
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



#' Visualisation of Set Enrichment Analysis Results from `SpaMTP::FindRegionalPathways()`.
#'
#' @param SpaMTP A `SpaMTP` Seurat object containing spatial metabolomics (SM) and/or spatial transcriptomics (ST) data. If SM data is included, it must be annotated using the `SpaMTP::AnnotateSM()` function.
#' @param ident A character string specifying the cluster identifier used to group regions, corresponding to a column name in the `SpaMTP@meta.data` slot.
#' @param DE.list A list containing differential expression results from the `FindAllMarkers()` function, with items matching the order of the `analyte_types` argument.
#' @param regpathway A dataframe generated by the `SpaMTP::FindRegionalPathways()` function, containing identified regional pathways.
#' @param selected_pathways A character vector specifying the names or IDs of pathways to be included in the analysis (e.g., `c("Amino acid metabolism", "WP1902", "Aspartate and asparagine metabolism")`). This argument is not case-sensitive.
#' @param pval_cutoff_pathway A numeric value between 0 and 1, defining the cutoff for the adjusted p-value from the permutation test used to compute significant pathways.
#' @param num_display An integer specifying the number of pathways to display in the plot.
#' @param text_size A numeric value controlling the size of the text elements in the plot.
#'
#' @return A `ggplot` object representing the set enrichment analysis results.
#' @export
#'
#' @examples
#' # PlotRegionalPathways(SpaMTP, ident = "clusters", regpathway = pathway_df)
PlotRegionalPathways <- function(SpaMTP,
                   ident,
                   regpathway,
                   selected_pathways = NULL,
                   pval_cutoff_pathway = NULL,
                   num_display = NULL,
                   text_size = NULL) {
  ## Checks for ident in SpaMTP Object
  if (!(ident %in% colnames(SpaMTP@meta.data))) {
    stop(
      "Ident: ",
      ident,
      " not found in SpaMTP object's @meta.data slot ... Make sure the ident column is in your @metadata and is a factor!"
    )
  }

  clus = as.factor(gsub("\\,.*", "", SpaMTP@meta.data[[ident]]))
  coordnate = sapply(
    SpaMTP@meta.data[, which(grepl(
      colnames(SpaMTP@meta.data),
      pattern = "coord",
      ignore.case = T
    ))],
    FUN = function(x) {
      as.numeric(gsub("\\,.*", "", x))
    }
  )
  palette <- suppressWarnings({
    brewer.pal(length(levels(clus)), "Set3")
  })
  if (length(levels(clus)) > length(palette)) {
    palette <- colorRampPalette(brewer.pal(9, "Set3"))(length(levels(clus)))
  }
  names(palette) = levels(clus)
  # SPM coordiante
  uid = unique(regpathway$Cluster_id)
  clu_names = data.frame(cluster = levels(clus),
                         clu_name = paste0("Cluster", levels(clus)))
  colnames(regpathway)[1] = "pathwayName"
  regpathway = regpathway %>% dplyr::group_by(pathwayName) %>% dplyr::mutate(
    Significance = ifelse(
      pval <= 0.05,
      "Significant at 5% significance level",
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
  # dendro <- ggtree(as.phylo(hc), layout = "rectangular")+scale_x_reverse()
  segment_hc <- with(ggdendro::segment(dendro_data(hc)),
                     data.frame(
                       x = y,
                       y = x,
                       xend = yend,
                       yend = xend
                     ))
  pos_table <- with(dendro_data(hc)$labels,
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
      labels = str_wrap(pos_table$path, width = 40),
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
          axis.text = element_text(size = ((text_size) %||% 12) - 3))
  #### ggplot
  suppressWarnings({
    gg_dot = ggplot(data = regpathway, aes(
      x = factor(Cluster_id, levels = sort(unique(Cluster_id))),
      y = factor(pathnameid, levels = dendro_data(hc)$labels$label)
    )) +
      geom_point(aes(colour = as.numeric(NES), size = as.numeric(sqrt(size)) * 5.1)) +
      scale_colour_gradient2(
        name = "Normalised enrichment score",
        low = "blue",
        # Original low color
        high = "red",
        # Original high color
        limits = c(min(regpathway$NES), max(regpathway$NES))
      ) +
      scale_size_continuous(name = "Number of altered metabolites \n in the pathway") +
      new_scale_colour() +
      geom_point(shape = 1,
                 aes(colour = Significance, size = as.numeric(size) * 2.1 +
                       0.1)) +
      scale_color_manual(
        values = c(
          "Significant at 5% significance level" = "green",
          "Not statistically significant" = "black"
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
      new_scale_colour() + theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey")
      ) +   theme(
        legend.position = "left",
        axis.text.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "horizontal"
      ) + xlab("Clusters")
  })
  u = ggarrange(gg_dot,
                plt_dendr,
                ncol = 2,
                widths = c(1, 0.6))
  print(u)
  svglite::svglite(file = "~/figure.svg",
                   width = 12,
                   height = 7)
  print(u)
  dev.off()
  return(u)
}


