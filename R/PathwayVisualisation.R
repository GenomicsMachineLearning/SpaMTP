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



#' This the function used to compute the exact fisher test for over-representation based pathway analysis
#'
#' @param SpaMTP A seurat object contains spatial metabolomics/transcriptomics features or both.
#' @param selected_pathways A character vector contains the candidate pathway names that need to be plotted (.
#' @param pval_cutoff_pathway A numerical value between 0 and 1 describe the cutoff adjusted p value for the permutation test used to compute output pathways (default = NULL).
#' @param num_display Number of pathways that to be output in the figures (default = NULL).
#' @param text_size A numerical value controls the text size of the plot (default = NULL).
#'
#' @return A combined gg, ggplot object displaying differentially expressed pathway and a cluster dendrogram
#' @export
#'
#' @import grid
#'
#' @examples
#' # PlotRegionalPathways(SpaMTP_obj)
PlotRegionalPathways = function(SpaMTP,
                   selected_pathways = NULL,
                   pval_cutoff_pathway = NULL,
                   num_display = NULL,
                   text_size = NULL) {
  #dendrogram based on jaccard distance
  enriched_names = names(SpaMTP@misc)[which(grepl(names(SpaMTP@misc), pattern = "set_enriched"))]
  repeat {
    # Prompt user for input
    cat(
      paste0(
        "Please enter the index of one the following (i.e. 1,2...) to specify the set enrichment result (in the format set_enriched+ The cluster used + background indices used) to use: \n",
        paste0(1:length(enriched_names), ".", enriched_names , collapse = " \n")
      )
    )
    user_input = readline(prompt = "Please enter:")
    # Check if user wants to exit
    user_input = as.numeric(user_input)
    if (user_input %in% 1:length(enriched_names)) {
      gsea_all_cluster_sig = SpaMTP@misc[[which(names(SpaMTP@misc) == enriched_names[user_input])]]
      break
    } else{
      cat(paste0(
        "\n Please enter correct one of followings:",
        paste0(enriched_names, collapse = "\n"),
        collapse = ";"
      ))
    }
  }
  gsea_all_cluster_sig = gsea_all_cluster_sig %>% filter(Cluster_id!="Cluster9")

  cluster_indent = sub("set_enriched_", "", sub("_[^_]*$", "", enriched_names[user_input]))
  clus = as.factor(gsub("\\,.*", "", SpaMTP@meta.data[[cluster_indent]]))
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
  palette <- brewer.pal(length(levels(clus)), "Set3")
  if (length(levels(clus)) > length(palette)) {
    palette <- colorRampPalette(brewer.pal(9, "Set3"))(length(levels(clus)))
  }
  names(palette) = levels(clus)
  # SPM coordiante
  uid = unique(gsea_all_cluster_sig$Cluster_id)
  clu_names = data.frame(cluster = levels(clus),
                         clu_name = paste0("Cluster", levels(clus)))
  colnames(gsea_all_cluster_sig)[1] = "pathwayName"
  gsea_all_cluster_sig = gsea_all_cluster_sig %>% dplyr::group_by(pathwayName) %>% dplyr::mutate(
    Significance = ifelse(
      pval <= 0.05,
      "Significant at 5% significance level",
      "Not statistically significant"
    )
  ) %>% dplyr::mutate(group_importance = sum(abs(NES)))
  if(is.null(selected_pathways )){
    gsea_all_cluster_sig = gsea_all_cluster_sig %>%  dplyr::filter(group_importance %in% sort(
      unique(gsea_all_cluster_sig$group_importance),
      decreasing = T
    )[1:(num_display %||% min(10,length(unique(gsea_all_cluster_sig$pathwayName))))])
  }else{
    gsea_all_cluster_sig = gsea_all_cluster_sig %>% dplyr::filter(pathwayName %in% selected_pathways)
    gsea_all_cluster_sig = gsea_all_cluster_sig %>%  dplyr::filter(group_importance %in% sort(
      unique(gsea_all_cluster_sig$group_importance),
      decreasing = T
    )[1:(num_display %||% min(10,length(unique(gsea_all_cluster_sig$pathwayName))))])
  }

  gsea_all_cluster_sig = gsea_all_cluster_sig %>% dplyr::mutate(pathnameid = paste0(pathwayName, "(", sourceId, ")"))
  pathwaynames = unique(gsea_all_cluster_sig$pathnameid)
  n = length(unique(gsea_all_cluster_sig$pathnameid))
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
    scale_x_continuous(expand = c(0, 0.1), limits = c(0, max(segment_hc$xend))) +
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
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = ((text_size) %||% 11)-3))
  #### ggplot
  suppressWarnings({
    gg_dot = ggplot(data = gsea_all_cluster_sig, aes(
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
        limits = c(
          min(gsea_all_cluster_sig$NES),
          max(gsea_all_cluster_sig$NES)
        )
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
        legend.title = element_text(size = ((text_size) %||% 8)-3),
        legend.text = element_text(size = ((text_size) %||% 8)-3)
      ) +
      new_scale_colour() + theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey")
      ) +   theme(
        legend.position = "left",
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "horizontal"
      )
  })

  gg_dotnl = gg_dot + theme(legend.position = "none")


  uid = c("Cluster9", "Cluster2")
  raster_image = list()
  for (j in 1:length(uid)) {
    temp_matrix = data.frame(coordnate) %>% mutate(colour = "red")
    if (sub("cluster", "", tolower(uid[j])) %in% tolower(clus)) {
      temp_matrix$colour = "red"
      temp_matrix$colour[which(tolower(clus) != sub("cluster", "", tolower(uid[j])))] = "#A9A9A9"
    } else{
      temp_matrix$colour = "red"
      temp_matrix$colour[which(tolower(clus) %in% sub("cluster", "", tolower(uid)))] = "#A9A9A9"
    }
    colnames(temp_matrix)[1] = "x"
    colnames(temp_matrix)[2] = "y"
    temp_matrix = na.omit(unique(temp_matrix))
    matrix_result <- as.matrix(reshape2::dcast(temp_matrix, y ~ x, value.var = "colour"))[, -1]
    raster_image[[j]] = as.raster(matrix_result)
  }

  ggano = ggplot(data = gsea_all_cluster_sig, aes(x = factor(Cluster_id, levels = sort(
    unique(Cluster_id)
  )), y = 1)) + geom_point() + ylim(-0.5, 0.5) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_blank()
  )
  gganot  = ggano + theme(axis.text.x = element_text(size = text_size %||% 12))

  g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x)
      x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }

  legend <- as.ggplot(g_legend(gg_dot))+ theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_blank()
  )

  for (k in 1:length(uid)) {
    gganot <- gganot + annotation_custom(
      rasterGrob(raster_image[[k]], interpolate = TRUE),
      xmin = k - 0.4,
      xmax = k + 0.4,
    )
  }

  suppressWarnings({left = ggarrange(gg_dotnl, gganot, nrow = 2, heights = c(3, 1))})
  suppressWarnings({right = ggarrange(plt_dendr,
                                      legend,
                                      nrow = 2,
                                      heights = c(3.65, 1))})
  ggarrange(left,right,ncol = 2)
  u = plot_grid(left, right, align = 'h', rel_widths = c(1, 1.5))

  print(u)
  return(u)
}
