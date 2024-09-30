library(RColorBrewer)
library(ggnewscale)
library(ggplotify)
library(ggpubr)
library(ggdendro)
library(cowplot)

#' This the function is used to check the normalisation/scaling/clustering/dimension reduction status
#'
#' @param SpaMTP A seurat object contains spatial metabolomics/transcriptomics features or both.
#' @param assay Character string defining the SpaMTP assay to extract intensity values from (default = "SPM").
#' @param resolution Goes to seurat::FindClusters Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param npcs A numerical integer, representing the number of principle components that the PCA will reduce to.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python
#' @param cluster_name cluster_name
#' @param ... The other parameters goes to Seurat methods
#'
#' @return A normalised/scaled/clustered and dimensionally reduced SpaMTP seurat object
#' @export
#'

check_spamtp = function(SpaMTP,
                        assay,
                        resolution = 0.5,
                        npcs = 30,
                        algorithm = 1,
                        verbose = T,
                        cluster_name = NULL,
                        ...) {
  DefaultAssay(SpaMTP) <- assay
  verbose_message(message_text = paste0("Checking ", assay, "...") ,
                  verbose = verbose)
  
  if (any(grepl(names(SpaMTP@commands), pattern = "NormalizeData") &
          grepl(names(SpaMTP@commands), pattern = assay))) {
    verbose_message(message_text = "Data normlaised, continue... " , verbose = verbose)
  } else{
    verbose_message(message_text = "Normalising data, continue... " , verbose = verbose)
    SpaMTP = NormalizeData(SpaMTP , verbose = FALSE, ...)
  }
  
  
  if (any(grepl(names(SpaMTP@commands), pattern = "ScaleData_") &
          grepl(names(SpaMTP@commands), pattern = assay))) {
    verbose_message(message_text = "Data scaled, continue... " , verbose = verbose)
  } else{
    verbose_message(message_text = "Scaling data, continue... " , verbose = verbose)
    SpaMTP = ScaleData(SpaMTP, verbose = FALSE, ...)
  }
  
  
  if (any(grepl(names(SpaMTP@commands), pattern = "FindVariableFeatures") &
          grepl(names(SpaMTP@commands), pattern = assay))) {
    verbose_message(message_text = "Data variable feature found, continue... " , verbose = verbose)
  } else{
    verbose_message(message_text = "Finding variable feature of data, continue... " , verbose = verbose)
    SpaMTP =  FindVariableFeatures(SpaMTP, verbose = FALSE, ...)
  }
  
  if (any(grepl(names(SpaMTP@commands), pattern = "RunPCA") &
          grepl(names(SpaMTP@commands), pattern = assay))) {
    verbose_message(message_text = "PCA found, continue..." , verbose = verbose)
  } else{
    verbose_message(
      message_text = paste0(
        "Runing PCA to reduced to ",
        npcs,
        " principle components, continue... "
      ) ,
      verbose = verbose
    )
    SpaMTP = RunPCA(
      SpaMTP,
      npcs = npcs,
      reduction.name = paste0(assay, "pca"),
      verbose = FALSE,
      ...
    )
  }
  # Creating neighbouts
  if (any(grepl(names(SpaMTP@commands), pattern = "FindNeighbors") &
          grepl(names(SpaMTP@commands), pattern = assay))) {
    verbose_message(message_text = "Nearest neighbour found, continue..." , verbose = verbose)
  } else{
    verbose_message(
      message_text = paste0("Founding Nearest neighbours, continue... ") ,
      verbose = verbose
    )
    #SpaMTP =  FindNeighbors(SpaMTP, dims = 1:npcs,graph.name = paste0(assay,".nn"), reduction = paste0(assay,"pca"), verbose = FALSE)
    SpaMTP =  FindNeighbors(
      SpaMTP,
      dims = 1:npcs,
      reduction = paste0(assay, "pca"),
      verbose = FALSE,
      ...
    )
  }
  
  if (any(grepl(names(SpaMTP@commands), pattern = "RunUMAP") &
          grepl(names(SpaMTP@commands), pattern = assay))) {
    verbose_message(message_text = "UMAP done, continue..." , verbose = verbose)
  } else{
    verbose_message(message_text = paste0("Computing UMAP, continue...") ,
                    verbose = verbose)
    SpaMTP =   RunUMAP(
      SpaMTP,
      dims = 1:npcs,
      reduction = paste0(assay, "pca"),
      reduction.name = paste0(assay, "umap"),
      verbose = FALSE,
      ...
    )
  }
  
  
  if (any(grepl(names(SpaMTP@commands), pattern = "FindClusters") &
          grepl(names(SpaMTP@commands), pattern = assay))) {
    verbose_message(message_text = "Clustering done, continue..." , verbose = verbose)
  } else{
    verbose_message(message_text = paste0("Clustering, continue... ") ,
                    verbose = verbose)
    SpaMTP = FindClusters(
      SpaMTP,
      resolution = resolution,
      algorithm = algorithm,
      cluster.name = cluster_name,
      verbose = FALSE,
      ...
    )
    #SpaMTP =   FindClusters(combined.data, resolution = resolution,algorithm =algorithm, graph.name = paste0(assay,".nn"),cluster.name = cluster_name, verbose = FALSE)
  }
  return(SpaMTP)
}

#' This the function used to compute the exact fisher test for over-representation based pathway analysis
#'
#' @param SpaMTP A seurat object contains spatial metabolomics/transcriptomics features or both.
#' @param polarity The polarity of the spatial metabolomics experiment. Inputs must be either NULL, 'positive' or 'negative'. If NULL, pathway analysis will run in neutral mode (default = NULL).
#' @param adduct A vector subset of SpaMTP::adduct_file$adduct_name, e.g. c("M+K","M+H "), default will take the full list of SpaMTP::adduct_file$adduct_name
#' @param analyte_types A subset of c("gene", "metabolites"), can be c("gene"), c("metabolites") or both
#' @param spatial_metabolomic_assay A Character string defining descrbing slot name for spatial metabolomics data in SpaMTP to extract intensity values from, default is "SPM"
#' @param spatial_transcriptomic_assay A Character string defining descrbing slot name for spatial transcriptomics data in SpaMTP to extract RNA count values from(default = "SPT").
#' @param resolution Goes to seurat::FindClusters Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python
#' @param slot The slot name in each assay to access the matrix data of each omics, default is "count"
#' @param npcs A numerical integer, representing the number of principle components that the PCA will reduce to.
#' @param st_cluster_name A character to describe the clustered result for spatial transcriptomics if the cluster haven't been done
#' @param sm_cluster_name A character to describe the clustered result for spatial metabolomics  if the cluster haven't been done"
#' @param max_path_size The max number of metabolites in a specific pathway (default = 500).
#' @param min_path_size The min number of metabolites in a specific pathway (default = 5).
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion `[ion mass,m/z]`/`[Full width at half height]` (default = 30000).
#' @param ppm_error A numerical value, standing for the parts-per-million error tolerance of matching m/z value with potential metabolites (default = NULL, calculated by )
#' @param cluster_vector A factor vector where each levels indicates the different clustered regions in tissue
#' @param background_cluster A vector consist of 0 and 1, where 1 indicates the intended background region
#' @param pval_cutoff_pathway A numerical value between 0 and 1 describe the cutoff adjusted p value for the permutation test used to compute output pathways
#' @param pval_cutoff_mets A numerical value between 0 and 1 describe the cutoff adjusted p value for the differential expression analysis for metabolites
#' @param pval_cutoff_genes A numerical value between 0 and 1 describe the cutoff adjusted p value for the differential expression analysis for RNAs
#' @param verbose A boolean value indicates whether verbose is shown
#' @param ... The other parameters goes to Seurat methods

#'
#' @return A SpaMTP object with set enrichment on given analyte types.
#' @export
#'
#' @examples
#' # SpaMTP = region_sea(SpaMTP, polarity = "positive")



region_sea = function(SpaMTP,
                      polarity,
                      adduct = NULL,
                      analyte_types = c("gene", "metabolites"),
                      spatial_metabolomic_assay = "SPM",
                      spatial_transcriptomic_assay = "SPT",
                      slot = "counts",
                      npcs = 30,
                      resolution = 0.5,
                      algorithm = 1,
                      test_use = "wilcox",
                      st_cluster_name = "ST_clusters",
                      sm_cluster_name = "SM_clusters",
                      tof_resolution = 30000,
                      min_path_size = 5,
                      max_path_size = 500,
                      cluster_vector = NULL,
                      background_cluster = NULL,
                      ppm_error = NULL,
                      pval_cutoff_pathway = NULL,
                      pval_cutoff_mets = NULL,
                      pval_cutoff_genes = NULL,
                      verbose = T,
                      ...) {
  if ("gene" %in% analyte_types) {
    if (is.null(SpaMTP@assays[[spatial_transcriptomic_assay]])) {
      stop(
        paste0(
          "If you are use genetice data with 'gene' in 'analyte_types' (e.g spatial transcriptolomics), please ensure that '",
          spatial_transcriptomic_assay,
          "' is inside the input seurat object"
        )
      )
    } else{
      SpaMTP = check_spamtp(
        SpaMTP,
        spatial_transcriptomic_assay,
        resolution = resolution,
        npcs = npcs,
        algorithm = algorithm,
        verbose = verbose,
        cluster_name = st_cluster_name,
        ...
      )
      gene_matrix = Matrix::t(SpaMTP[[spatial_transcriptomic_assay]]@layers[[slot]])
    }
  }
  if ("metabolites" %in% analyte_types) {
    if (is.null(SpaMTP@assays[[spatial_metabolomic_assay]])) {
      stop(
        paste0(
          "If you are use metabolic data with 'metabolites' in 'analyte_types' (e.g spatial metabolomics), please ensure that '",
          spatial_metabolomic_assay,
          "' is inside the input seurat object"
        )
      )
    } else{
      SpaMTP = check_spamtp(
        SpaMTP,
        spatial_metabolomic_assay,
        resolution = resolution,
        npcs = npcs,
        algorithm = algorithm,
        verbose = verbose,
        cluster_name = sm_cluster_name,
        ...
      )
      mass_matrix = Matrix::t(SpaMTP[[spatial_metabolomic_assay]]@layers[[slot]])
    }
  }
  empty_theme = theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
  if (!is.null(cluster_vector)) {
    cluster_vector = as.factor(cluster_vector)
    
    if ((("metabolites" %in% analyte_types)&length(cluster_vector)!= nrow(mass_matrix)) |
        (("gene" %in% analyte_types)&length(cluster_vector)!= nrow(gene_matrix)) ) {
      stop(
        "Please make sure the input cluster is a vector has same length as the number of elements and is a factor"
      )
    } else{
      
      assignment = cluster_vector
      cluster = levels(cluster_vector)
      sscs = "Custome_clustering"
      user_input = 1
      sscs[user_input]
      SpaMTP@meta.data[[sscs[user_input]]] = assignment
    }
  } else{
    verbose_message(message_text = "Checking for whether spatial Shrunken centroid or clustering is finisher" , verbose = verbose)
    sscs = names(SpaMTP@meta.data)[which(
      grepl(names(SpaMTP@meta.data), pattern = "ssc") |
        grepl(names(SpaMTP@meta.data), pattern = st_cluster_name) |
        grepl(names(SpaMTP@meta.data), pattern = sm_cluster_name)
    )]
    if (length(sscs) == 0) {
      stop(
        "No spatial shrunken centroid done, please either run spatial shrunken centroid or provide a clustering index (cluster)"
      )
    } else{

      gg_cluster = lapply(sscs, function(x) {
        assignment = as.factor(gsub("\\,.*", "", SpaMTP@meta.data[which(names(SpaMTP@meta.data) == x)][, 1]))
        assignment[which(is.na(assignment))] = sample(levels(assignment), size =
                                                        1)
        cluster = levels(assignment)
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
        palette <- brewer.pal(length(cluster), "Set3")
        if (length(cluster) > length(palette)) {
          palette <- colorRampPalette(brewer.pal(9, "Set3"))(length(cluster))
        }
        names(palette) = cluster
        
        gg = ggplot() + geom_raster(aes(
          x = coordnate[, 1],
          y = coordnate[, 2],
          fill = assignment
        )) + scale_fill_manual(values = palette) + empty_theme + ggtitle(x)
        return(gg)
      })
      suppressWarnings({
        print(do.call(grid.arrange, c(gg_cluster, nrow = 2)))
      })
      repeat {
        # Prompt user for input
        user_input = readline(
          prompt = paste0(
            "Please enter the index of one the following (i.e. 1,2...) to specify the spatial neighborhood radius of nearby pixels to consider: \n",
            paste0(1:length(sscs), ".", sscs, collapse = " \n")
          )
        )
        # Check if user wants to exit
        user_input = as.numeric(user_input)
        if (user_input %in% 1:length(sscs)) {
          cat(paste0("Selected ", sscs[user_input]))
          assignment = as.factor(gsub("\\,.*", "", SpaMTP@meta.data[which(names(SpaMTP@meta.data) == sscs[user_input])][, 1]))
          assignment[which(is.na(assignment))] = sample(levels(assignment), size =
                                                          1)
          cluster = levels(assignment)
          break
        } else{
          cat(paste0(
            "\n Please enter correct one of followings:",
            paste0(sscs, collapse = "\n"),
            collapse = ";"
          ))
        }
      }
    }
    assignment = as.factor(gsub("\\,.*", "", SpaMTP@meta.data[which(names(SpaMTP@meta.data) == sscs[user_input])][, 1]))
  }
  
  
  
  #(1) Get the rank entry for each cluster

  
  
  if (is.null(background_cluster)) {
    #assignment = as.factor(gsub("\\,.*", "", SpaMTP@meta.data[which(names(SpaMTP@meta.data) == sscs[user_input])][, 1]))
    assignment[which(is.na(assignment))] = sample(levels(assignment), size =
                                                    1)
    cluster = levels(assignment)
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
    palette <- brewer.pal(length(cluster), "Set3")
    if (length(cluster) > length(palette)) {
      palette <- colorRampPalette(brewer.pal(9, "Set3"))(length(cluster))
    }
    names(palette) = cluster
    
    gg = ggplot() + geom_raster(aes(x = coordnate[, 1], y = coordnate[, 2], fill = assignment)) + scale_fill_manual(values = palette) +
      empty_theme + ggtitle(ifelse(exists("sscs"),sscs[user_input],"Custom cluster"))
    suppressWarnings({
      print(gg)
    })
    
    background_cluster = rep(0, times = length(assignment))
    sel_bacground = c("pseudo")
    background_clu = c()
    repeat {
      # Prompt user for input
      cc = c(cluster, "pseudo")
      verbose_message(
        message_text = paste0(
          "Select any regions can be considered as background/Baseline signals(if no need to selected, type 'none', if all background regions selected, type 'done'): \n",
          paste0(c(cc[-which(cc %in% sel_bacground)], "none", "done"), collapse =
                   " \n")
        ) ,
        verbose = verbose
      )
      user_input_background = readline(prompt = "")
      background_clu = c(background_clu, user_input_background)
      if (user_input_background  %in% c("none", "done")) {
        break
      } else if (user_input_background %in% cc) {
        sel_bacground = c(sel_bacground, user_input_background)
        background_cluster[which(assignment == user_input_background)] = 1
      } else{
        verbose_message(
          paste0(
            "\n Please enter correct one of followings:\n",
            paste0(c(cc[-which(cc %in% sel_bacground)], "none", "done"), collapse =
                     " \n"),
            collapse = ";"
          ),
          verbose = verbose
        )
      }
    }
    background_cluster = as.factor(background_cluster)
  } else{
    background_cluster = as.factor(background_cluster)
    if ((!all(levels(background_cluster) %in% c(0, 1))) |
        (length(background_cluster) != length(assignment))) {
      stop(
        "Please ensure that background_cluster is a vecter of 0 and 1 where 1 indicates the background region"
      )
    }
  }
  
  
  # (2) Annotation
  if (!is.null(SpaMTP@assays[[spatial_metabolomic_assay]])) {
    if (!all(c("pathway_db", "annotated_db") %in% names(SpaMTP@assays[[spatial_metabolomic_assay]]@misc))) {
      tryCatch({
        input_mz = data.frame(cbind(
          row_id = 1:ncol(mass_matrix),
          mz = as.numeric(
            stringr::str_extract(row.names(SpaMTP[[spatial_metabolomic_assay]]@features), pattern = "\\d+\\.?\\d*")
          )
        ))
      }, error = function(cond) {
        stop(
          "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
        )
      }, warning = function(cond) {
        stop(
          "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
        )
      })
      # Set the db that you want to search against
      db = rbind(HMDB_db, Chebi_db, Lipidmaps_db)
      # set which adducts you want to search for
      #load("data/adduct_file.rda")
      
      if (polarity == "positive") {
        test_add_pos = adduct_file$adduct_name[which(adduct_file$charge > 0)]
        test_add_pos = test_add_pos[which(test_add_pos %in% (adduct %||% adduct_file$adduct_name))]
        # Using Chris' pipeline for annotation
        # 1) Filter DB by adduct.
        db_1 = db_adduct_filter(db,
                                test_add_pos,
                                polarity = "pos",
                                verbose = verbose)
      } else if (polarity == "negative") {
        test_add_neg = adduct_file$adduct_name[which(adduct_file$charge < 0)]
        test_add_neg = test_add_neg[which(test_add_neg %in% (adduct %||% adduct_file$adduct_name))]
        # Using Chris' pipeline for annotation
        # 1) Filter DB by adduct.
        db_1 = db_adduct_filter(db,
                                test_add_neg,
                                polarity = "neg",
                                verbose = verbose)
      } else if (polarity == "neutral") {
        # Using Chris' pipeline for annotation
        # 1) Filter DB by adduct.
        db_1 = db %>% mutate("M" = `M-H ` + 1.007276)
      } else{
        stop("Please enter correct polarity from: 'positive', 'negative', 'neutral'")
      }
      db_2 = formula_filter(db_1)
      ppm_error = ppm_error %||% (1e6 / tof_resolution / sqrt(2 * log(2)))
      db_3 = proc_db(input_mz, db_2, ppm_error) %>% dplyr::mutate(entry = stringr::str_split(Isomers, pattern = "; "))
      verbose_message(message_text = "Query necessary data and establish pathway database" , verbose = verbose)
      input_id = lapply(db_3$entry, function(x) {
        x = unlist(x)
        index_hmdb = which(grepl(x, pattern = "HMDB"))
        x[index_hmdb] = paste0("hmdb:", x[index_hmdb])
        index_chebi = which(grepl(x, pattern = "CHEBI"))
        x[index_chebi] = tolower(x[index_chebi])
        return(x)
      })
      
      db_3 = db_3 %>% dplyr::mutate(inputid = input_id) %>%  dplyr::mutate(chem_source_id = input_id)
      rampid = c()
      chem_source_id = unique(chem_props$chem_source_id)
      
      verbose_message(message_text = "Query db for addtional matching" , verbose = verbose)
      
      db_3 = merge(chem_props, db_3, by = "chem_source_id")
      
      verbose_message(message_text = "Query finished!" , verbose = verbose)
      # get names for the ranks
      name_rank = lapply(input_mz$mz, function(x) {
        return(unique(na.omit(db_3[which(db_3$observed_mz == x), ])))
      })
      
      
      name_rank_ids = lapply(input_mz$mz, function(x) {
        return(unique(na.omit(db_3$rampid[which(db_3$observed_mz == x)])))
      })
      # Get pathway db
      verbose_message(message_text = "Constructing pathway database ..." , verbose = verbose)
      chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")
      
      pathway_db = split(chempathway$rampId, chempathway$pathwayName)
      pathway_db = pathway_db[which(!duplicated(names(pathway_db)))]
      pathway_db = pathway_db[lapply(pathway_db, length) >= min_path_size  &
                                lapply(pathway_db, length) <= max_path_size]
      
      SpaMTP@assays[[spatial_metabolomic_assay]]@misc = list(pathway_db = pathway_db, annotated_db =   db_3)
    } else{
      pathway_db = SpaMTP@assays[[spatial_metabolomic_assay]]@misc[["pathway_db"]]
      db_3 =   SpaMTP@assays[[spatial_metabolomic_assay]]@misc[["annotated_db"]]
    }
  }
  ###########################
  verbose_message(message_text = "Calculating metabolomics differentially expressed cluster ..." , verbose = verbose)
  indices = assignment
  levels(indices) = c(levels(indices), "Baseline")
  background_clu = background_clu[-which(background_clu == "done" |
                                           background_clu == "none")]
  indices[which(indices %in% background_clu)] = "Baseline"
  u = paste0(sscs[user_input], "_Baseline")
  SpaMTP@meta.data[[u]] = indices
  
  Idents(SpaMTP) = u
  DE_met <- FindAllMarkers(SpaMTP,
                           assay = "SPM",
                           only.pos = F,
                           test.use = test_use)
  db_3 = db_3 %>% mutate(mz_name = paste0("mz-", db_3$observed_mz))
  DE_met = DE_met %>% mutate(mz_name = gene)
  db_3 = merge(db_3 , DE_met, by = "mz_name")
  verbose_message(message_text = "Calculating transcriptomics differentially expressed cluster ..." , verbose = verbose)
  DE_rna <- FindAllMarkers(SpaMTP,
                           assay = "SPT",
                           only.pos = F,
                           test.use = test_use
  )
  DE_rna = DE_rna %>% mutate(commonName = toupper(gene))
  source_gene = merge(DE_rna, source_df[which(grepl(source_df$rampId,
                                                    pattern = "RAMP_G")),], by = "commonName")
  
  SpaMTP@misc[[paste0("DE_",
                      spatial_metabolomic_assay,
                      "_",
                      sscs[user_input],
                      "_",
                      paste0(c(background_clu, "Baseline"), collapse = "."))]] = db_3
  SpaMTP@misc[[paste0("DE_",
                      spatial_transcriptomic_assay,
                      "_",
                      sscs[user_input],
                      "_",
                      paste0(c(background_clu, "Baseline"), collapse = "."))]] = source_gene
  
  gc()
  
  gsea_all_cluster = data.frame()
  all_ranks = list()
  non_bac_cluster = na.omit(unique(indices[which(background_cluster == 0)]))
  pb3 = txtProgressBar(
    min = 0,
    max = ifelse(
      all(background_cluster == 0),
      length(non_bac_cluster),
      (length(non_bac_cluster) + 1)
    ),
    initial = 0,
    style = 3
  )
  for (i in 1:(length(non_bac_cluster) + 1)) {
    # Consider the background pathway
    ranks = c()
    if (i == (length(non_bac_cluster) + 1)) {
      if (all(background_cluster == 0)) {
        next
      } else{
        clu_wise = which(background_cluster == 1)
      }
      # Get coordinates for the elements in the cluster
      sub_db3 = db_3[which(db_3$cluster == "Baseline"), ] %>% filter(p_val_adj <= pval_cutoff_mets %||% 0.05) %>% filter(!duplicated(ramp_id))
      ranks = scale(sub_db3$avg_log2FC, center = 0)
      names(ranks) = sub_db3$ramp_id
      
      sub_de_gene = source_gene[which(source_gene$cluster == "Baseline"), ] %>% filter(p_val_adj <= pval_cutoff_genes %||% 0.05) %>% filter(!duplicated(rampId))
      ranks_gene_vector = scale(sub_de_gene$avg_log2FC, center = 0)
      names(ranks_gene_vector) = sub_de_gene$rampId
      # Genes and metabolites
      ranks = c(ranks, ranks_gene_vector)
      ####################################################### non background
    } else{
      sub_db3 = db_3[which(as.character(db_3$cluster) == as.character(non_bac_cluster[i])), ] %>% filter(p_val_adj <= pval_cutoff_mets %||% 0.05) %>% filter(!duplicated(ramp_id))
      ranks = sub_db3$avg_log2FC
      names(ranks) = sub_db3$ramp_id
      
      sub_de_gene = source_gene[which(as.character(source_gene$cluster) == as.character(non_bac_cluster[i])), ] %>% filter(p_val_adj <= pval_cutoff_genes %||% 0.05) %>% filter(!duplicated(rampId))
      ranks_gene_vector = scale(sub_de_gene$avg_log2FC, center = 0)
      names(ranks_gene_vector) = sub_de_gene$rampId
      # Genes and metabolites
      ranks = c(ranks, ranks_gene_vector)
    }
    ranks = ranks[which(!duplicated(names(ranks)))]
    all_ranks[[i]] = ranks[is.finite(ranks)]
    names(all_ranks)[i] = paste0("cluster", i)
    suppressWarnings({
      gsea_result = fgsea::fgsea(
        pathways =  pathway_db,
        stats = ranks,
        minSize = 5,
        maxSize = 500
      )  %>%  dplyr::mutate(Cluster_id = ifelse(
        i <= length(non_bac_cluster),
        paste0("Cluster", non_bac_cluster[i]),
        "Baseline"
      ))
    })
    gsea_result = na.omit(gsea_result)
    short_source = source_df[which((source_df$rampId %in% names(ranks)) &
                                     !duplicated(source_df$rampId)), ]
    addtional_entry = do.call(rbind, lapply(1:nrow(gsea_result), function(x) {
      temp = unique(unlist(gsea_result$leadingEdge[x]))
      temp_ref = db_3[which(db_3$ramp_id %in% temp), ] %>% dplyr::mutate(adduct_info = paste0(observed_mz, "[", Adduct, "]")) %>% dplyr::filter(!duplicated(adduct_info))
      temp_rna = short_source[which((short_source$rampId %in% temp) &
                                      (grepl(short_source$rampId, pattern = "RAMP_G"))), ]
      return(
        data.frame(
          adduct_info = paste0(temp_ref$adduct_info, collapse = ";"),
          leadingEdge_metabolites = paste0(temp_ref$IsomerNames, collapse = ";"),
          leadingEdge_genes = paste0(temp_rna$commonName, collapse = ";"),
          regulation = paste0(ifelse(ranks[which(names(ranks) %in% temp)] >= 0, "↑", "↓"), collapse = ";")
        )
      )
    }))
    gsea_result = cbind(gsea_result , addtional_entry)
    gsea_all_cluster = rbind(gsea_all_cluster, gsea_result)
    setTxtProgressBar(pb3, i)
  }
  close(pb3)
  gsea_all_cluster_sig = gsea_all_cluster %>% dplyr::group_by(pathway) %>% dplyr::filter(any(as.numeric(pval) <= (pval_cutoff_pathway %||% 0.05))) %>% dplyr::mutate(
    Significance = ifelse(
      pval <= 0.05,
      "Significant at 5% significance level",
      "Not statistically significant"
    )
  ) %>% dplyr::mutate(group_importance = sum(abs(NES)))
  gsea_all_cluster_return =   gsea_all_cluster %>%  dplyr::mutate(
    Significance = ifelse(
      pval <= 0.05,
      "Significant at 5% significance level",
      "Not statistically significant"
    )
  )
  colnames(gsea_all_cluster_return)[1] = "pathwayName"
  gsea_all_cluster_return = merge(gsea_all_cluster_return, pathway, by = "pathwayName")
  ########################################################
  #Plot#
  SpaMTP@misc[[paste0("set_enriched_", sscs[user_input], "_", paste0(c(background_clu, "baseline"), collapse =
                                                                       "."))]] = gsea_all_cluster_return
  return(SpaMTP)
}


#' This the function used to compute the exact fisher test for over-representation based pathway analysis
#'
#' @param SpaMTP A seurat object contains spatial metabolomics/transcriptomics features or both.
#' @param selected_pathways A character vector contains the candidate pathway names that need to be plotted.
#' @param pval_cutoff_pathway A numerical value between 0 and 1 describe the cutoff adjusted p value for the permutation test used to compute output pathways
#' @param num_display Number of pathways that to be output in the figures
#' @param text_size A numerical value contols the text size of the plot


plotsea = function(SpaMTP,
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

