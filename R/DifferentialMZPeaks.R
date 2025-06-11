
#### SpaMTP Differential Metabolite Analysis Functions ########################################################################################################################################################################################



#' Pools SpaMTP Seurat object into random pools for pseudo-bulking.
#'
#' Runs pooling of a SpaMTP dataset to generate pseudo-replicates for each unique identity provided.
#' This function is used by `FindAllDEMs()`.
#'
#' @param data.filt A Seurat Object containing count values for pooling.
#' @param idents A character string defining the idents column to pool the data against.
#' @param n An integer defining the amount of pseudo-replicates to generate for each sample (default = 3).
#' @param assay Character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param slot Character string defining the assay storage slot to pull the relative mz intensity values from (default = "counts").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A SinglCellExpereiment object which contains pooled (n)-pseudo-replicate counts data based on the Seurat Object input
#' @export
#'
#' @examples
#' # run_pooling <- list(seuratObj, idents = "sample", n = 3, assay = "Spatial", slot = "counts")
run_pooling <- function(data.filt, idents, n, assay, slot, verbose = TRUE) {

  cell_metadata <- data.filt@meta.data
  samples <- unique(cell_metadata[[idents]])

  verbose_message(message_text = paste0("Pooling one sample into ", n ," replicates..."), verbose = verbose)

  nrg <- n
  for(i in c(1:length(samples))){
    set.seed(i)
    wo<-which(cell_metadata[[idents]]== samples[i])
    cell_metadata[wo,'orig.ident2']<-paste(samples[i],sample(c(1:n),length(wo)
                                                             ,replace=T,prob=rep(1/nrg,nrg)),sep='_')
  }
  gene_data <- row.names(data.filt)
  filtered.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data.filt[[assay]][slot]),
                                       colData = cell_metadata)


  tempf=strsplit(filtered.sce@colData[["orig.ident2"]],'_')
  pid=NULL
  for(i in 1:length(tempf)){
    pidone=tempf[[i]]
    if(length(pidone)!=3){
      pidone=c(pidone[1],'yes',pidone[2])
    }
    pid=rbind(pid,pidone)
  }

  filtered.sce@colData$type=pid[,2]

  summed <- scater::aggregateAcrossCells(filtered.sce,
                                 id=SingleCellExperiment::colData(filtered.sce)[,'orig.ident2'])

  return(summed)
}




#' Runs EdgeR analysis for pooled data
#'
#' Worker function for calculating differentially abundant metabolites per pooling group.
#' This function is used by by `FindAllDEMs()`.
#'
#' @param pooled_data A SingleCellExperiment object which contains the pooled pseudo-replicate data.
#' @param seurat_data A Seurat object containing the merged Xenium data being analysed (this is subset).
#' @param ident A character string defining the ident column to perform differential expression analysis against.
#' @param output_dir A character string defining the ident column to perform differential expression analysis against.
#' @param run_name A character string defining the title of this DE analysis (will be used when saving DEMs to .csv file).
#' @param n An integer that defines the number of pseudo-replicates per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param annotation.column Character string defining the column where annotation information is stored in the assay metadata. This requires AnnotateSeuratMALDI() to be run where the default column to store annotations is "all_IsomerNames" (default = "None").
#' @param assay A character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#' @param return.individual Boolean value defining whether to return a list of individual edgeR objects for each designated ident. If FALSE, one merged edgeR object will be returned (default = FALSE).
#'
#' @returns A modified edgeR object which contains the relative pseudo-bulking analysis outputs, including a DEMs data.frame with a list of differential expressed m/z metabolites
#' @export
#'
#' @examples
#' # pooled_obj <- run_pooling(SeuratObj, "sample", n = 3)
#' # run_DE(pooled_obj, SeuratObj, "sample", "~/Documents/DE_output/", "run_1", n = 3, logFC_threshold = 1.2, annotation.column = "all_IsomerNames", assay = "Spatial")
run_DE <- function(pooled_data, seurat_data, ident, output_dir, run_name, n, logFC_threshold, annotation.column, assay, return.individual = FALSE, verbose = TRUE){

  verbose_message(message_text = paste("Running edgeR DE Analysis for ", run_name, " -> with samples [", paste(unique(unlist(seurat_data@meta.data[[ident]])), collapse = ", "), "]"), verbose = verbose)

  annotation_result <- list()

  for (condition in unique(seurat_data@meta.data[[ident]])) {

    # Create groups
    groups <- SingleCellExperiment::colData(pooled_data)[[ident]]
    groups <- gsub(condition, "Comp_A", groups)
    groups <- ifelse(groups != "Comp_A", "Comp_B", groups)

    # Extract continuous expression data (e.g., intensity matrix)
    expression_data <- SingleCellExperiment::counts(pooled_data)  # Or the assay holding your continuous data

    # Optional: If your data is raw intensities, log-transform it here (add small offset if needed)
    expression_data <- log2(expression_data + 1)

    #keep <- rowMeans(expression_data) > some_threshold  # Define threshold based on your data
    #expression_data <- expression_data[keep, ]

    # Create design matrix
    design <- model.matrix(~groups)
    design[, 2] <- 1 - design[, 2]  # To match your original contrast logic

    # Fit linear model
    fit <- lmFit(expression_data, design)

    # Empirical Bayes moderation
    fit <- eBayes(fit, robust = TRUE)

    # Use treat for log fold change threshold testing if desired
    res <- treat(fit, lfc = log2(logFC_threshold), robust = TRUE)

    all_decisions <- decideTests(res)[, ncol(decideTests(res))]

    res_table <- topTreat(fit, coef = ncol(fit$design), n = nrow(expression_data), lfc = log2(logFC_threshold))

    res_table$regulate <- dplyr::recode(
      as.character(all_decisions[rownames(res_table),]),
      "0" = "Normal",
      "1" = "Up",
      "-1" = "Down"
    )

    # Order by p-value or FDR
    de_group_limma <- res_table[order(res_table$adj.P.Val), ]

    # Rename adj.P.Val to FDR
    colnames(de_group_limma) <- ifelse(colnames(de_group_limma) == "adj.P.Val", "FDR", colnames(de_group_limma))

    # Add gene/metabolite names
    de_group_limma$gene <- rownames(de_group_limma)

    # Add annotations if requested
    if (!is.null(annotation.column)) {
      annotation.data <- seurat_data[[assay]]@meta.data
      if (!(annotation.column %in% colnames(annotation.data))) {
        stop("Warning: The annotation column provided does not exist in seurat_data[[assay]]@meta.data")
      } else {
        rownames(annotation.data) <- annotation.data$mz_names
        annotation.data_subset <- annotation.data[rownames(de_group_limma), ]
        de_group_limma$annotations <- annotation.data_subset[[annotation.column]]
      }
    }

    # Write CSV output if directory specified
    if (!is.null(output_dir)) {
      utils::write.csv(de_group_limma, file.path(output_dir, paste0(condition, "_", run_name, ".csv")))
    }

    # Store results
    y$DEMs <- de_group_limma
    annotation_result[[condition]] <- y

    verbose_message(message_text = paste("Analysis complete for condition:", condition), verbose = verbose)

  }



  if (return.individual){
    annotation_result <- lapply(names(annotation_result), function(x){
      annotation_result[[x]]$DEMs$cluster <- x
      annotation_result[[x]]
    })

    return(annotation_result)
  } else {

    edger <- edgeR::DGEList(
      counts = annotation_result[[1]]$counts,
      samples = annotation_result[[1]]$samples
    )
    edger$samples$group <- edger$samples$ident
    edger$samples$condition <- NULL

    dems <- lapply(names(annotation_result), function(x){
      annotation_result[[x]]$DEMs$cluster <- x
      rownames(annotation_result[[x]]$DEMs) <- NULL
      annotation_result[[x]]$DEMs
    })

    combined_dems <- do.call(rbind, dems)
    rownames(combined_dems) <- 1:length(combined_dems$cluster)

    edger$DEMs <- combined_dems
    return(edger)
  }

}


#' Finds differentially expressed m/z values/metabolites between all comparison groups.
#'
#' @param data A Seurat object containing mz values for differential expression analysis.
#' @param ident A character string defining the metadata column or groups to compare mz values between.
#' @param n An integer that defines the number of pseudo-replicates (pools) per sample (default = 3).
#' @param logFC_threshold A numeric value indicating the logFC threshold to use for defining significant genes (default = 1.2).
#' @param DE_output_dir A character string defining the directory path for all output files to be stored. This path must a new directory. Else, set to NULL as default.
#' @param run_name A character string defining the title of this DE analysis that will be used when saving DEMs to .csv file (default = 'FindAllDEMs').
#' @param annotation.column Character string defining the column where annotation information is stored in the assay metadata. This requires AnnotateSeuratMALDI() to be run where the default column to store annotations is "all_IsomerNames" (default = "None").
#' @param assay A character string defining the assay where the mz count data and annotations are stored (default = "Spatial").
#' @param slot Character string defining the assay storage slot to pull the relative mz intensity values from. Note: EdgeR requires raw counts, all values must be positive (default = "counts").
#' @param return.individual Boolean value defining whether to return a list of individual edgeR objects for each designated ident. If FALSE, one merged edgeR object will be returned (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns Returns an list() contains the EdgeR DE results. Pseudo-bulk counts are stored in $counts and DEMs are in $DEMs.
#' @export
#'
#' @examples
#' # FindAllDEMs(SeuratObj, "sample",DE_output_dir = "~/Documents/DE_output/", annotations = TRUE)
FindAllDEMs <- function(data, ident, n = 3, logFC_threshold = 1.2, DE_output_dir = NULL, run_name = "FindAllDEMs", annotation.column = NULL, assay = "Spatial", slot = "counts", return.individual = FALSE, verbose = TRUE){

  if (!(is.null(DE_output_dir))){
    if (dir.exists(DE_output_dir)){
      warning("Please supply a directory path that doesn't already exist")
      stop("dir.exists(DE_output_dir) = TRUE")
    } else{
      dir.create(DE_output_dir)
    }
  }

  #Step 1: Run Pooling to split each unique ident into 'n' number of pseudo-replicate pools
  pooled_data <- run_pooling(data,ident, n = n, assay = assay, slot = slot, verbose = verbose)

  #Step 2: Run EdgeR to calculate differentially expressed m/z peaks
  DEM_results <- run_DE(pooled_data, data, ident = ident, output_dir = DE_output_dir, run_name = run_name, n=n, logFC_threshold=logFC_threshold, annotation.column = annotation.column, assay = assay, verbose = verbose, return.individual = return.individual)

  # Returns an EDGEr object which contains the pseudo-bulk counts in $counts and DEMs in $DEMs
  return(DEM_results)

}




#' Heatmap of Differentially Expressed Metabolites
#'
#' Generates a heatmap of DEMs generated from edgeR analysis run using `FindAllDEMs()`.
#' This function uses `pheatmap` to plot data.
#'
#' @param edgeR_output A list containing outputs from edgeR analysis (from FindAllDEMs()). This includes pseudo-bulked counts and DEMs.
#' @param n A numeric integer that defines the number of UP and DOWN regulated peaks to plot (default = 25).
#' @param only.pos Boolean indicating if only positive markers should be returned (default = FALSE).
#' @param FDR.threshold Numeric value that defines the FDR threshold to use for defining most significant results (default = 0.05).
#' @param logfc.threshold Numeric value that defines the logFC threshold to use for filtering significant results (default = 0.5).
#' @param order.by Character string defining which parameter to order markers by, options are either 'FDR' or 'logFC' (default = "FDR").
#' @param scale A character string indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param color A vector of colors used in heatmap (default = grDevices::colorRampPalette(c("navy", "white", "red"))(50)).
#' @param cluster_cols Boolean value determining if columns should be clustered or hclust object (default = F).
#' @param cluster_rows Boolean value determining if rows should be clustered or hclust object (default = T).
#' @param fontsize_row A numeric value defining the fontsize of rownames (default = 15).
#' @param fontsize_col A numeric value defining the fontsize of colnames (default = 15).
#' @param cutree_cols A numeric value defining the number of clusters the columns are divided into, based on the hierarchical clustering(using cutree), if cols are not clustered, the argument is ignored (default = 9).
#' @param silent Boolean value indicating if the plot should not be draw (default = TRUE).
#' @param plot_annotations_column Character string indicating the column name that contains the metabolite annotations to plot. Annotations = TRUE must be used in FindAllDEMs() for edgeR output to include annotations. If plot_annotations_column = NULL, m/z vaues will be plotted (default = NULL).
#' @param save_to_path Character string defining the full filepath and name of the plot to be saved as.
#' @param plot.save.width Integer value representing the width of the saved pdf plot (default = 20).
#' @param plot.save.height Integer value representing the height of the saved pdf plot (default = 20).
#' @param nlabels.to.show Numeric value defining the number of annotations to show per m/z (default = NULL).
#' @param annotation_colors List for specifying annotation_row and annotation_col track colors manually. Check pheatmap R-Package documentation for details. If set to 'NA', default coloring will be used (default = NA).
#'
#' @returns A heatmap plot of significantly differentially expressed metabolites defined in the edgeR ouput object.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' # DEMs <- FindAllDEMs(SeuratObj, "sample")
#'
#' # DEMsHeatmap(DEMs)
DEMsHeatmap <- function(edgeR_output,
                         n = 5,
                         only.pos = FALSE,
                         FDR.threshold = 0.05,
                         logfc.threshold = 0.5,
                         order.by = "FDR",
                         scale ="row",
                         color = grDevices::colorRampPalette(c("navy", "white", "red"))(50),
                         cluster_cols = F,
                         cluster_rows = T,
                         fontsize_row = 15,
                         fontsize_col = 15,
                         cutree_cols = 9,
                         silent = TRUE,
                         plot_annotations_column = NULL,
                         save_to_path = NULL,
                         plot.save.width = 20,
                         plot.save.height = 20,
                         nlabels.to.show = NULL,
                         annotation_colors = NULL){


  degs <- edgeR_output$DEMs
  degs <- subset(degs, FDR < FDR.threshold)

  if (order.by == "FDR"){

    grouped_pos<- degs %>%
      group_by(cluster) %>%
      filter( logFC > logfc.threshold) %>%
      arrange(desc(regulate)) %>%
      slice_head(n = n)


    if (only.pos) {
      grouped_neg <- NULL

    } else {
      grouped_neg <- degs %>%
        group_by(cluster) %>%
        filter(logFC < - logfc.threshold) %>%
        arrange(regulate) %>%
        slice_head(n = n)
    }
    df <- do.call(rbind, list(grouped_pos,grouped_neg))
    df <- df[order(df$cluster, dplyr::desc(df$regulate)), ]

  } else {
    if ( order.by != "logFC"){
      warning("order.by has invalid argument. Must be either 'FDR' or 'logFC'. Heatmap defaulting to order by logFC")
    }

    grouped_pos<- degs %>%
      group_by(cluster) %>%
      filter(logFC > logfc.threshold) %>%
      arrange(-logFC) %>%
      slice_head(n = n)


    if (only.pos) {
      grouped_neg <- NULL
    } else {
      grouped_neg <- degs %>%
        group_by(cluster) %>%
        filter(logFC < - logfc.threshold) %>%
        arrange(logFC) %>%
        slice_head(n = n)
    }
    df <- do.call(rbind, list(grouped_pos,grouped_neg))
    df <- df[order(df$cluster, -df$logFC), ]
  }



  col_annot <- data.frame(sample = edgeR_output$samples$ident)
  row.names(col_annot) <- colnames(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE)))

  if (!is.null(annotation_colors)){
    annotation_colors <- list(sample = unlist(annotation_colors))
  } else {
    annotation_colors <- NA
  }

  mtx <- as.matrix(as.data.frame(edgeR::cpm(edgeR_output,log=TRUE))[unique(df$gene),])
  if (!(is.null(plot_annotations_column))){
    if (is.null(edgeR_output$DEMs[[plot_annotations_column]])){
      warning("There are no annotations present in the edgeR_output object. Run 'annotate.SeuratMALDI()' prior to 'FindAllDEMs' and set annotations = TRUE .....\n Heatmap will plot default m/z values ... ")
    } else{
      if (!is.null(nlabels.to.show)){
        df[[plot_annotations_column]] <- labels_to_show(df[[plot_annotations_column]], n = nlabels.to.show)
      }
      rownames(mtx) <- unique(df[[plot_annotations_column]])
    }
  }

  p <- pheatmap::pheatmap(mtx,scale=scale,color=color,cluster_cols = cluster_cols, annotation_col=col_annot, cluster_rows = cluster_rows,
                          fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_cols = cutree_cols, silent = silent, annotation_colors = annotation_colors)

   if (!(is.null(save_to_path))){
     save_pheatmap_as_pdf(pheatmap = p, filename = save_to_path, width = plot.save.width, height = plot.save.height)
   }

  return(p)
}


#' Saves a DEMsHeatmap as a PDF
#'
#' @param pheatmap A pheatmap plot object that is being saved.
#' @param filename Character string defining the full filepath and name of the plot to be saved as.
#' @param width Integer value representing the width of the saved pdf plot (default = 20).
#' @param height Integer value representing the height of the saved pdf plot (default = 20).
#'
#' @export
#'
#' @examples
#' # save_pheatmap_as_pdf(pheatmap, filename = "/Documents/plots/pheatmap1")
save_pheatmap_as_pdf <- function(pheatmap, filename, width=20, height=20){

  pdf(paste0(filename,".pdf"), width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(pheatmap$gtable)
  dev.off()
}





########################################################################################################################################################################################################################
