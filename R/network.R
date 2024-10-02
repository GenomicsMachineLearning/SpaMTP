#' This the function used to compute the exact fisher test for over-representation based pathway analysis
#'
#' @param SpaMTP A seurat object contains spatial metabolomics/transcriptomics features or both, after processing set enrichment via SpaMTP::region_sea
#' @param candidate_pathways A char vector includes the pathway names that are used to construct the networks, e.g c("Amino acid metabolism","Aspartate and asparagine metabolism"), not case sensitive.
#' @param path The directory to write the output, default is current working directory
#' @param spatial_metabolomic_assay A Character string defining descrbing slot name for spatial metabolomics data in SpaMTP to extract intensity values from, default is "SPM"
#' @param spatial_transcriptomic_assay A Character string defining descrbing slot name for spatial transcriptomics data in SpaMTP to extract RNA count values from(default = "SPT")
#' @param slot The slot name in each assay to access the matrix data of each omics, default is "count"
#' @param analyte_types A subset of c("gene", "metabolites"), can be c("gene"), c("metabolites") or both


#'
#' @return Am interactive html with network structure of given candidate pathways
#' @export
#'
#' @examples
#' spa_network(SpaMTP, c("3-Phosphoglycerate dehydrogenase deficiency","Cystathionine Beta-Synthase Deficiency","DNA Repair"))

spa_network  = function(SpaMTP,
                        candidate_pathways,
                        path = getwd(),
                        slot = "counts",
                        colour_palette = NULL,
                        spatial_metabolomic_assay = "SPM",
                        spatial_transcriptomic_assay = "SPT",
                        analyte_types = c("gene", "metabolites")) {
  if ("gene" %in% analyte_types) {
    if (is.null(SpaMTP@assays[[spatial_transcriptomic_assay]])) {
      stop(
        paste0(
          "If you are use genetice data with 'gene' in 'analyte_types' (e.g spatial transcriptolomics), please ensure that '",
          spatial_transcriptomic_assay,
          "' is inside the input seurat object"
        )
      )
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
    }
  }
  if (("metabolites" %in% analyte_types) &
      ("gene" %in% analyte_types)) {
    if (nrow(mass_matrix) != nrow(gene_matrix)) {
      stop("Please align the spatial data before processing")
    }
  }
  
  repeat {
    # Prompt user for input
    
    cat(paste0(
      "Please check the cluster ident used for this function: \n",
      paste0(1:length(names(
        SpaMTP@meta.data
      )), ".", names(SpaMTP@meta.data), collapse = " \n"),
      "\n"
    ))
    user_input = readline(prompt = paste0("Enter a number to select one of above: \n"))
    # Check if user wants to exit
    user_input = as.numeric(user_input)
    if (user_input %in% 1:length(names(SpaMTP@meta.data))) {
      cat(paste0("Selected ident:", names(SpaMTP@meta.data)[user_input]))
      assignment = as.factor(gsub("\\,.*", "", SpaMTP@meta.data[which(names(SpaMTP@meta.data) == names(SpaMTP@meta.data)[user_input])][, 1]))
      assignment[which(is.na(assignment))] = sample(levels(assignment), size =
                                                      1)
      cluster = levels(assignment)
      break
    } else{
      cat(" \n")
      cat(paste0(
        "\n Please enter correct one of followings: \n",
        paste0(
          1:length(names(SpaMTP@meta.data)),
          ".",
          names(SpaMTP@meta.data),
          collapse = " \n"
        )
      ))
    }
  }
  
  
  if (any(grepl(
    names(SpaMTP@misc),
    pattern = names(SpaMTP@meta.data)[user_input],
    ignore.case = T
  ))) {
    if ("metabolites" %in% analyte_types) {
      met_can_de_list = names(SpaMTP@misc)[which((grepl(
        names(SpaMTP@misc),
        pattern = names(SpaMTP@meta.data)[user_input],
        ignore.case = T
      )) &  (
        grepl(
          names(SpaMTP@misc),
          pattern = spatial_metabolomic_assay,
          ignore.case = T
        )
      ))]
      if (length(met_can_de_list) == 0) {
        stop(
          'Please ensure that differential expression analysis for metabolomics has been done or remove "metabolites" from analyte_types'
        )
      }
      repeat {
        # Prompt user for input
        cat(
          paste0(
            "Please check the name of differential expression for metabolomics data: \n",
            paste0(
              1:length(met_can_de_list),
              ".",
              met_can_de_list,
              collapse = " \n"
            ),
            "\n"
          )
        )
        user_input2 = readline(prompt = paste0("Enter a number to select one of above: \n"))
        # Check if user wants to exit
        user_input2 = as.numeric(user_input2)
        if (user_input2 %in% 1:length(met_can_de_list)) {
          de_mets = SpaMTP@misc[[which(names(SpaMTP@misc) == met_can_de_list[user_input2])]]
          break
        } else{
          cat(" \n")
          cat(paste0("\nPlease enter correct one of followings: \n"))
        }
      }
    }
    if ("gene" %in% analyte_types) {
      rna_can_de_list = names(SpaMTP@misc)[which((grepl(
        names(SpaMTP@misc),
        pattern = names(SpaMTP@meta.data)[user_input],
        ignore.case = T
      )) &  (
        grepl(
          names(SpaMTP@misc),
          pattern = spatial_transcriptomic_assay,
          ignore.case = T
        )
      ))]
      if (length(rna_can_de_list) == 0) {
        stop(
          'Please ensure that differential expression analysis for transcriptolomics has been done or remove "gene" from analyte_types'
        )
      }
      repeat {
        # Prompt user for input
        cat(
          paste0(
            "Please check the name of differential expression for transcriptolomics data: \n",
            paste0(
              1:length(rna_can_de_list),
              ".",
              rna_can_de_list ,
              collapse = " \n"
            ),
            "\n"
          )
        )
        user_input3 = readline(prompt = paste0("Enter a number to select one of above: \n"))
        # Check if user wants to exit
        user_input3 = as.numeric(user_input3)
        if (user_input3 %in% 1:length(rna_can_de_list)) {
          de_rna = SpaMTP@misc[[which(names(SpaMTP@misc) == rna_can_de_list[user_input3])]]
          break
        } else{
          cat(" \n")
          cat(paste0("\nPlease enter correct one of followings: \n"))
        }
      }
    }
    
  } else{
    stop(
      "Please check whether DE analysis accross given ident has been done prior to this function"
    )
  }
  
  
  # get enrichment dataframe
  
  candidate_enriched = names(SpaMTP@misc)[which(grepl(
    names(SpaMTP@misc),
    pattern = "set_enriched",
    ignore.case = T
  ))]
  if (length(candidate_enriched) == 0) {
    stop("Please run region_sea prior to this function")
  }
  repeat {
    # Prompt user for input
    cat(
      paste0(
        "Please select which dataframe to be used for network construction: \n",
        paste0(
          1:(length(candidate_enriched) + 1),
          ".",
          c(candidate_enriched, "default") ,
          collapse = " \n"
        ),
        "\n"
      )
    )
    user_input4 = readline(prompt = paste0("Enter a number to select one of above: \n"))
    # Check if user wants to exit
    user_input4 = as.numeric(user_input4)
    if (user_input4 %in% 1:(length(candidate_enriched) + 1)) {
      if (user_input4 == (length(candidate_enriched) + 1)) {
        default_enrich = candidate_enriched[which(grepl(
          candidate_enriched,
          pattern = names(SpaMTP@meta.data)[user_input],
          ignore.case = T
        ))]
        verbose_message(
          message_text = paste0("Selected default dataframe:", default_enrich),
          verbose = verbose
        )
        enriched_df = SpaMTP@misc[[which(names(SpaMTP@misc) == default_enrich)]] %>% dplyr::filter(sub("Cluster", "", Cluster_id) %in% as.character(cluster))
        returnname = default_enrich
      } else{
        enriched_df = SpaMTP@misc[[which(names(SpaMTP@misc) == candidate_enriched[user_input4])]] %>% dplyr::filter(sub("Cluster", "", enriched_df$Cluster_id) %in% as.character(cluster))
        returnname = candidate_enriched[user_input4]
      }
      break
    } else{
      cat(" \n")
      cat(paste0("\nPlease enter correct one of followings: \n"))
    }
  }
  
  SpatialColors <- grDevices::colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
  colour_palette = colour_palette %||% SpatialColors(100)
  
  
  node_ind = paste0("const networks = [")
  sub_enriched = enriched_df[which(tolower(enriched_df$pathwayName) %in% tolower(candidate_pathways)), ]
  
  
  type = sub_enriched$type
  
  check_topology_existence = function(type, pathway_of_interest) {
    temp_db = c()
    if (type == "wiki") {
      index = which(tolower(names(humanwiki_ramp)) == tolower(pathway_of_interest))
      if (length(index) == 0) {
        temp_db = humanwiki_ramp
      }
    } else if (type == "reactome") {
      index = which(tolower(names(humanReactome_ramp)) == tolower(pathway_of_interest))
      if (length(index) != 0) {
        temp_db = humanReactome_ramp
      }
    } else if (type == "kegg") {
      index = which(tolower(names(humankegg_ramp)) == tolower(pathway_of_interest))
      if (length(index) != 0) {
        temp_db = humankegg_ramp
      }
    } else if (type == "hmdb") {
      index = which(tolower(names(humansmp_ramp)) == tolower(pathway_of_interest))
      if (length(index) != 0) {
        temp_db = humansmp_ramp
      }
    } else{
      all_list = c(humanwiki_ramp,
                   humanReactome_ramp,
                   humankegg_ramp,
                   humansmp_ramp)
      all_names = names(all_list)
      index = which(tolower(all_names) == tolower(pathway_of_interest))
      if (length(index) == 0) {
        return(NULL)
      } else{
        return(all_list[index])
      }
    }
    if (exists("temp_db")) {
      return(temp_db[[index]])
    } else{
      return(NULL)
    }
  }
  
  topodb = apply(
    sub_enriched[!duplicated(sub_enriched$pathwayName), ],
    MARGIN = 1,
    FUN = function(x) {
      check_topology_existence(x$type, x$pathwayName)
    }
  )
  dblengths = unlist(lapply(topodb, length))
  names(topodb) = sub_enriched$pathwayName[!duplicated(sub_enriched$pathwayName)]
  if (any(dblengths == 0)) {
    verbose_message(
      message_text = paste0(
        "Following pathway(s) cannot find matched topological structure: \n",
        paste0(sub_enriched[!duplicated(sub_enriched$pathwayName), ]$pathwayName[which(dblengths == 0)], collapse = ";\n")
      ),
      verbose = verbose
    )
  }
  topodb = topodb[which(dblengths != 0)]
  
  
  
  sub_enriched = sub_enriched[which(tolower(sub_enriched$pathwayName) %in% tolower(names(topodb))), ]
  pathway_names = unique(sub_enriched$pathwayName)
  
  network = paste0('const networks = [')
  matrix_ids = c()
  ucid = naturalsort::naturalsort(unique(sub_enriched$Cluster_id))
  
  simplified_content = c()
  for (i in 1:length(ucid)) {
    sub_cluster = sub_enriched[which(sub_enriched$Cluster_id == ucid[i]), ]
    # cluster wise metabolites
    sub_expr_met = de_mets[which(tolower(de_mets$cluster) ==  str_extract(tolower(ucid [i]), "[0-9]+")), ]
    # cluster wise rna
    sub_expr_rna = de_rna[which(tolower(de_rna$cluster) == str_extract(tolower(ucid[i]), "[0-9]+")), ]
    cluster = paste0("[")
    
    temp_analytes = unique(unlist(lapply(sub_cluster$leadingEdge, function(x) {
      return(unlist(x))
    })))
    simplified_content = c(simplified_content, temp_analytes)
    for (j in 1:length(pathway_names)) {
      # if pathway result present in given cluster
      if (tolower(pathway_names[j]) %in% tolower(sub_cluster$pathwayName)) {
        temp_topo = topodb[[which(unique(tolower(names(topodb))) == tolower(pathway_names[j]))]]
        
        
        
        path_df = na.omit(
          merge(unique(
            rbind(temp_topo[["mixedEdges"]], temp_topo[["metabolEdges"]], temp_topo[["metabolPropEdges"]], temp_topo[["protPropEdges"]])
          ), as.data.frame(reaction_type), by = "reaction_type") %>% rowwise() %>% mutate(detected = ifelse((src %in% temp_analytes) |
                                                                                                              (dest %in% temp_analytes), T, F
          ))
        )
        path_df = path_df %>%   add_count(src, name = "src_n") %>%   add_count(dest, name = "dest_n")
        pathway_analyte = unique(c(path_df$src, path_df$dest))
        
        #simplified
        simplified_pathdf =  path_df %>% dplyr::filter((src %in% temp_analytes) |
                                                         (dest %in% temp_analytes))
        simplified_pathway_analyte = unique(c(simplified_pathdf$src, simplified_pathdf$dest))
        
        # Add plottable candidtable analytes
        matrix_ids = c(matrix_ids, temp_analytes[which(temp_analytes %in% pathway_analyte)])
        # 1.1 Complete node sets
        temp_nodes = paste0("{nodes: [")
        for (k in 1:length(pathway_analyte)) {
          if (pathway_analyte[k] %in% temp_analytes) {
            met_expr = sub_expr_met[which(sub_expr_met$ramp_id == pathway_analyte[k]), ]
            met_expr = met_expr[which.min(met_expr$p_val_adj), ]
            rna_expr = sub_expr_rna[which(sub_expr_rna$rampId == pathway_analyte[k]), ]
            rna_expr = rna_expr[which.min(rna_expr$p_val_adj), ]
            temp_nodes  = paste0(
              temp_nodes,
              '{ id:"',
              pathway_analyte[k],
              '", group:"',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                'rna',
                'mets'
              ),
              '",',
              'expr:',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                rna_expr$avg_log2FC,
                met_expr$avg_log2FC
              ),
              ', shape: "',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                'rna',
                'mets'
              ),
              '", name: "',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                rna_expr$gene,
                met_expr$common_name
              ),
              '", size:',
              sqrt(abs(
                ifelse(
                  grepl(pathway_analyte[k], pattern = "RAMP_G"),
                  rna_expr$avg_log2FC,
                  met_expr$avg_log2FC
                )
              )),
              ', display:"',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                rna_expr$gene,
                paste0(
                  met_expr$mz_name,
                  '[',
                  met_expr$Adduct,
                  ']-',
                  met_expr$common_name
                )
              ),
              '"},'
            )
          } else{
            temp_nodes  = paste0(
              temp_nodes,
              '{ id:"',
              pathway_analyte[k],
              '", group:"',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                'rna',
                'mets'
              ),
              '",',
              'expr: "Not"',
              ', shape: "',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                'rna',
                'mets'
              ),
              '", name: "',
              source_df$commonName[which(source_df$rampId == pathway_analyte[k])][1],
              '", size:',
              0.7,
              "},"
            )
          }
        }
        temp_nodes = paste0(temp_nodes, "],")
        
        
        # 2.1 Complete link sets
        temp_links = paste0("links: [")
        
        for (z in 1:nrow(path_df)) {
          temp_path_df = path_df[z, ]
          temp_links = paste0(
            temp_links,
            '{source:',
            '"',
            temp_path_df$src,
            '", target:',
            '"',
            temp_path_df$dest,
            '"',
            ',type:"',
            temp_path_df$colour,
            '", weight:"',
            temp_path_df$dest_n + temp_path_df$src_n,
            '", style:"',
            temp_path_df$linetype,
            '" ,head:"',
            temp_path_df$arrowhead,
            '"},'
          )
        }
        temp_links  = paste0(temp_links, "]},")
        
        
        temp_this_path = paste0(temp_nodes, temp_links)
        cluster  = paste0(cluster, temp_this_path)
      } else{
        # if pathway result not present in given cluster
        temp_nodes  = paste0("{nodes: [],")
        temp_links  = paste0("links: []},")
        temp_this_path = paste0(temp_nodes, temp_links)
        cluster  = paste0(cluster, temp_this_path)
      }
    }
    cluster = paste0(cluster, "],")
    network  = paste0(network, cluster)
  }
  network = paste0(network, "];")
  
  
  tab_div = c()
  for (o in 1:length(pathway_names)) {
    tab_div = paste0(tab_div, '<div class="tab">', pathway_names[o], "</div>")
  }
  
  
  
  
  options = paste0('<option value="option1" selected>', ucid[1], "</option>")
  
  for (k in 2:length(ucid)) {
    options = paste0(options,
                     '<option value="option',
                     k,
                     '">',
                     ucid[k],
                     "</option>")
  }
  
  scale_legend = as.integer(sqrt(max(abs(
    c(de_mets$avg_log2FC, de_rna$avg_log2FC)
  ))))
  
  
  # Get coordinates
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
  non_na_ind = which((!is.na(coordnate[, 1])) &
                       (!is.na(coordnate[, 2])))
  coordnate = cbind(coordnate, assign = as.character(assignment))
  coordnate = na.omit(coordnate)
  max_x =  max(na.omit(as.numeric(coordnate[, 1])))
  max_y =  max(na.omit(as.numeric(coordnate[, 2])))
  
  
  coordi = paste0("const coordinates = [")
  
  for (t in 1:nrow(coordnate)) {
    coordi = paste0(
      coordi,
      '[',
      as.numeric(coordnate[t, 1]) * 180 / max_x ,
      ',',
      as.numeric(coordnate[t, 2]) * 200 / max_y,
      '],'
    )
  }
  
  coordi = paste0(coordi, '];')
  
  
  # USE matrix_ids  to obtain the m/z matrix/ genetric matrix to be displayed
  met_plot = paste0("const metplot = {")
  rna_plot = paste0("const rnaplot = {")
  for (a in 1:length(matrix_ids)) {
    # for metabolites
    if (grepl(matrix_ids[a], pattern = "RAMP_C")) {
      if ("metabolites" %in% analyte_types) {
        met_temp_sub  = de_mets[which(de_mets$ramp_id == matrix_ids[a]), ]
        mz_name = met_temp_sub$gene[which.min(met_temp_sub$p_val_adj)]
        met_mat_index = which(tolower(rownames(SpaMTP@assays[[spatial_metabolomic_assay]]@features)) ==  mz_name)
        mmat_return = mass_matrix[non_na_ind, met_mat_index]
        met_plot = paste0(met_plot,
                          '"',
                          matrix_ids[a],
                          '":[',
                          paste0(mmat_return, collapse = ","),
                          '],')
      }
    } else{
      if ("gene" %in% analyte_types) {
        rna_temp_sub = de_rna$commonName[which(de_rna$rampId == matrix_ids[a])]
        gene_mat_index = which(tolower(rownames(SpaMTP@assays[[spatial_transcriptomic_assay]]@features)) == tolower(unique(rna_temp_sub)[1]))
        gmat_return = gene_matrix[non_na_ind, gene_mat_index]
        rna_plot = paste0(rna_plot,
                          '"',
                          matrix_ids[a],
                          '":[',
                          paste0(gmat_return, collapse = ","),
                          '],')
      }
    }
  }
  rna_plot = paste0(rna_plot, "};")
  met_plot = paste0(met_plot, "};")
  
  
  # Cluster colour
  cluster_infor = paste0('const cluster_info = ["',
                         paste0(coordnate[, 3], collapse = '","'),
                         '"]')
  
  html = paste0(
    '
<!DOCTYPE html>
  <html lang="en">

    <head>
    <meta charset="UTF-8">
      <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Interactive Network Visualization</title>
        <style>
        body {
          display: flex;
          margin: 0;
          font-family: Arial, sans-serif;
        }

      #sidebar {
      width: 200px;
      background-color: #f4f4f4;
        padding: 10px;
      height: 92vh;
      overflow-y: auto;
      border-right: 1px solid #ddd;
      }

#network-container {
flex: 1;
position: relative;
height: 94vh;
padding: 0px;
}

#network-frame {
border: 2px solid #333;
padding: 0px;
height: 94vh;
background-color: #f9f9f9;
  }


.legend {
  top: 10px;
  right: 10px;
  background-color: #fff;
    border: 1px solid #ddd;
  padding: 10px;
  width: 200px;
  position: flex;
}

.legend-rectangle {
  width: 90%;
  height: 20px;
  background: linear-gradient(to right, blue,yellow,red);
  margin-bottom: 5px;
  margin: 0 15px;
}

.legend-labels {
  display: flex;
  justify-content: space-between;
  /* Space out labels */
    width: 100%;
  /* Match the width of the gradient */
    margin-top: 5px;
  /* Space above the labels */
}

.tab {
  padding: 10px;
  cursor: pointer;
  border-bottom: 1px solid #ddd;
}

.tab:hover {
  background-color: #eee;
}

svg {
  width: 100%;
  height: 100%;
  display: block;
  margin: 5px 0;
}

.link {
  stroke-width: 2px;
}

.node circle {
  stroke: #fff;
    stroke-width: 1.5px;
}

.node text {
  font: 12px sans-serif;
  pointer-events: none;
}

.arrowhead {
  fill: none;
  stroke-width: 2px;
}

.scale-bar {
  margin: 10px;
}

.scale-bar label {
  display: block;
  margin-bottom: 5px;
}

.scale-bar input[type="range"] {
  width: 100%;
}

.line-legend {
  display: flex;
  align-items: center;
  margin: 5px 0;
}


.line-text {
  font-size: 12px;
}

.tab {
  padding: 10px;
  cursor: pointer;
}
.rastercontainer {
            display: flex;
            align-items: flex-start;
        }
        .rasterCanvas {
            border: 1px solid black;
        }
        .legendCanvas {
            margin-left: 0px;
            padding: 0px;
            border: 1px solid transparent;
        }

.tab.active {
  background-color: #007BFF;
    color: white;
}

.line {
  width: 80px;
  height: 2px;
  background-color: black;
  margin-right: 10px;
}

.shape {
  width: 20px;
  height: 20px;
  margin-right: 10px;
  display: inline-block;
}

.shape.triangle {
  width: 0;
  height: 0;
  border-left: 10px solid transparent;
  border-right: 10px solid transparent;
  border-bottom: 20px solid black;
}

.shape.rectangle {
  background-color: black;
}

.shape.circle {
  border-radius: 50%;
  background-color: black;
}
p {
      display: inline-block;
      max-width: 200px;
      font-size: 12px;
}
.hex1::before {
  width: 0;
  height: 20px;
  border-left: 0px solid transparent;
  border-right: 8px solid transparent;
  border-bottom: 0px solid transparent;
  border-top: -5px solid transparent;
  content: "B22";
  color: black;
  font-size: 35px;
}

.slider {
  width: 100px;
  height: 40px;
  background-color: #FFFF00;
  border-radius: 20px;
  position: relative;
  cursor: pointer;
  display: flex;
  justify-content: center;
  align-items: center;
  transition: background-color 0.3s ease;
}

.slider.active {
  background-color: #4CAF50;
}

.shape-text {
  font-size: 12px;
}

.scrollable-select {
  width: 200px;
  height: 100px;
  /* Height determines the visible area */
    overflow-y: auto;
  /* Add vertical scroll */
}
</style>
  </head>
  <body>
  <div id="sidebar">',
tab_div,
'
  </div>
  <div id="network-container">
  <div id="network-frame">
  <svg id="network">
  </svg>
  </div>

  </div>

  <div class="legend">
  <h4>Legend: log2 fold change</h4>
  <p><span style="color: red;">●</span> Upregulated</p>
  <p><span style="color: blue;">●</span> Downregulated</p>
  <p><span style="color: gray;">●</span> Not detected</p>

  <div>
  <div class="legend-rectangle" id="colorLegend">
  </div>
  <div class="legend-labels">
  <span>',
-scale_legend ,
'</span>
    <span>',
-scale_legend / 2,
'</span>
  <span>0</span>
    <span>',
scale_legend / 2,
'</span>
    <span>',
scale_legend ,
'</span>
  </div>
  </div>
  <hr>


  <div class = "line-legend">
  <canvas id="arrowCanvas1" width="110" height="20"></canvas>
  <div class="line-text">Indirect</div>
  </div>
  <div class = "line-legend">
  <canvas id="arrowCanvas2" width="110" height="20"></canvas>
  <div class="line-text">Binding</div>
  </div>
  <div class = "line-legend">
  <canvas id="arrowCanvas10" width="110" height="20"></canvas>
  <div class="line-text">Dissociation</div>
  </div>

  <div class = "line-legend">
  <canvas id="arrowCanvas3" width="110" height="20"></canvas>
  <div class="line-text">Missing</div>
  </div>
  <div class = "line-legend">
  <canvas id="arrowCanvas4" width="110" height="20"></canvas>
  <div class="line-text">Inhibition</div>
  </div>
  <div class = "line-legend">
  <canvas id="arrowCanvas5" width="110" height="20"></canvas>
  <div class="line-text">Phosphorylationn</div>
  </div>
  <div class = "line-legend">
  <canvas id="arrowCanvas9" width="110" height="20"></canvas>
  <div class="line-text">Dephosphorylation</div>
  </div>


  <div class = "line-legend">
  <canvas id="arrowCanvas6" width="110" height="20"></canvas>
  <div class="line-text">Activation</div>
  </div>
  <div class = "line-legend">
  <canvas id="arrowCanvas7" width="110" height="20"></canvas>
  <div class="line-text">Expression</div>
  </div>
  <div class = "line-legend">
  <canvas id="arrowCanvas12" width="110" height="20"></canvas>
  <div class="line-text">Repression</div>
  </div>

  <div class = "line-legend">
  <canvas id="arrowCanvas8" width="110" height="20"></canvas>
  <div class="line-text">Ubiquitination</div>
  </div>


  <div class = "line-legend">
  <canvas id="arrowCanvas11" width="110" height="20"></canvas>
  <div class="line-text">Methylation</div>
  </div>



  <div class = "line-legend">
  <canvas id="arrowCanvas13" width="110" height="20"></canvas>
  <div class="line-text">Glycosylation</div>
  </div>

  <div class = "line-legend">
  <canvas id="arrowCanvas14" width="110" height="20"></canvas>
  <div class="line-text">State change</div>
  </div>

  <hr>
  <div class="line-legend">
  <div class="shape rectangle"></div>
  <div class="shape-text">Metabolites</div>
  </div>
  <div class="line-legend">
  <div class="shape circle"></div>
  <div class="shape-text">RNA</div>
  </div>
  <hr>
  </div>
  <div class="scale-bar">
  <label for="node-scale">Node Size:</label>
  <input type="range" id="node-scale" min="0.5" max="3" step="0.1" value="1">
    <label for="link-scale">Link Size:</label>
  <input type="range" id="link-scale" min="0.2" max="2" step="0.05" value="1">
  <div class="slider" id="slider" onclick="toggleState()">Simplified net</div>
  <label for="select1">Select cluster for network:</label>
  <select id="select1" class="scrollable-select" size="5">
',
options,
'
  </select>
  <div class="rastercontainer">
    <canvas id = "rc1" class="rasterCanvas" width="160" height="200"></canvas>
    <canvas id= "lca1" class = "legendCanvas" width="40" height="200"></canvas>
</div>
<p width ="200" display="inline-block"  id= "spatial_window">Display spatial data by clicking on vertices</p>
<div class="rastercontainer">
  <canvas id = "rc2" class="rasterCanvas" width="160" height="200"></canvas>
</div>
<p width ="200" display="inline-block"  id= "clu_window">Cluster 1</p>
  </div>
  <script src="https://d3js.org/d3.v7.min.js"></script>
  <script>
  function drawArrowLine(ctx, startX, startY, endX, endY, headLength,angle ) {
    const dx = endX - startX;
    const dy = endY - startY;
    // Draw shaft of the arrow (the main line)
    ctx.beginPath();
    ctx.moveTo(startX, startY);
    ctx.lineTo(endX, endY);
    ctx.stroke();

    // Draw arrowhead]
var arrowAngle2 = -angle;
const arrowX1 = endX - headLength * Math.cos(angle);
const arrowY1 = endY - headLength * Math.sin(angle);
const arrowX2 = endX - headLength * Math.cos(arrowAngle2);
const arrowY2 = endY - headLength * Math.sin(arrowAngle2);

// Draw the first side of the arrowhead
ctx.moveTo(endX, endY);
ctx.lineTo(arrowX1, arrowY1);

// Draw the second side of the arrowhead
ctx.moveTo(endX, endY);
ctx.lineTo(arrowX2, arrowY2);

ctx.stroke();
  }

// Set line width and color

// Draw arrow line with 3 segments: 1 shaft + 2 arrowhead lines
const canvas = document.getElementById("arrowCanvas1");
const ctx1 = canvas.getContext("2d");
ctx1.lineWidth = 2;
ctx1.strokeStyle = "black";
ctx1.setLineDash([8, 10]);
drawArrowLine(ctx1,0, 10, 180-80, 10, 10,Math.PI/4);

const canvas2 = document.getElementById("arrowCanvas2");
const ctx2 = canvas2.getContext("2d");
ctx2.lineWidth = 2;
ctx2.strokeStyle = "blue";
ctx2.setLineDash([]);
drawArrowLine(ctx2,0, 10, 175-80, 10, 10,Math.PI*1.3);


const canvas3 = document.getElementById("arrowCanvas3");
const ctx3 = canvas3.getContext("2d");
ctx3.lineWidth = 2;
ctx3.strokeStyle = "black";
ctx3.setLineDash([]);
drawArrowLine(ctx3,0, 10, 175-80, 10, 8,Math.PI*1.3);
drawArrowLine(ctx3,0, 10, 175-80, 10, 8,Math.PI*0.3);

const canvas4 = document.getElementById("arrowCanvas4");
const ctx4 = canvas4.getContext("2d");
ctx4.lineWidth = 2;
ctx4.strokeStyle = "black";
ctx4.setLineDash([]);
drawArrowLine(ctx4,0, 10, 175-80, 10, 10,Math.PI/2);

const canvas5 = document.getElementById("arrowCanvas5");
const ctx5 = canvas5.getContext("2d");
ctx5.lineWidth = 2;
ctx5.strokeStyle = "orange";
ctx5.setLineDash([]);
drawArrowLine(ctx5,0, 10, 175-80, 10, 10,Math.PI/4);

const canvas6 = document.getElementById("arrowCanvas6");
const ctx6 = canvas6.getContext("2d");
ctx6.lineWidth = 2;
ctx6.strokeStyle = "red";
ctx6.setLineDash([8, 10]);
drawArrowLine(ctx6,0, 10, 175-80, 10, 30,Math.PI/4);

const canvas7 = document.getElementById("arrowCanvas7");
const ctx7 = canvas7.getContext("2d");
ctx7.lineWidth = 2;
ctx7.strokeStyle = "green";
ctx7.setLineDash([8, 10]);
drawArrowLine(ctx7,0, 10, 175-80, 10, 30,Math.PI/4);


const canvas8 = document.getElementById("arrowCanvas8");
const ctx8 = canvas8.getContext("2d");
ctx8.lineWidth = 2;
ctx8.strokeStyle = "purple";
ctx8.setLineDash([8, 10]);
drawArrowLine(ctx8,0, 10, 175-80, 10, 10,Math.PI/4);


const canvas9 = document.getElementById("arrowCanvas9");
const ctx9 = canvas9.getContext("2d");
ctx9.lineWidth = 2;
ctx9.strokeStyle = "orange";
ctx9.setLineDash([8, 10]);
drawArrowLine(ctx9,0, 10, 175-80, 10, 10,Math.PI/4);

const canvas10 = document.getElementById("arrowCanvas10");
const ctx10 = canvas10.getContext("2d");
ctx10.lineWidth = 2;
ctx10.strokeStyle = "black";
ctx10.setLineDash([8, 10]);
drawArrowLine(ctx10,0, 10, 175-80, 10, 10,Math.PI*1.3);

const canvas11 = document.getElementById("arrowCanvas11");
const ctx11 = canvas11.getContext("2d");
ctx11.lineWidth = 2;
ctx11.strokeStyle = "red";
ctx11.setLineDash([]);
drawArrowLine(ctx11,0, 10, 175-80, 10, 10,Math.PI/4);


const canvas12 = document.getElementById("arrowCanvas12");
const ctx12 = canvas12.getContext("2d");
ctx12.lineWidth = 2;
ctx12.strokeStyle = "black";
ctx12.setLineDash([8, 10]);
drawArrowLine(ctx12,0, 10, 175-80, 10, 10,Math.PI/2);

const canvas13 = document.getElementById("arrowCanvas13");
const ctx13 = canvas13.getContext("2d");
ctx13.lineWidth = 2;
ctx13.strokeStyle = "black";
ctx13.setLineDash([]);
drawArrowLine(ctx13,0, 10, 169-80, 10, 10,Math.PI);
ctx13.beginPath();
ctx13.arc(95, 10, 5, 0, Math.PI * 2);
ctx13.stroke();


const canvas14 = document.getElementById("arrowCanvas14");
const ctx14 = canvas14.getContext("2d");
ctx14.lineWidth = 2;
ctx14.strokeStyle = "black";
ctx14.setLineDash([]);
drawArrowLine(ctx14,0, 10,170-70, 10, 10,Math.PI/4);
',
coordi,
met_plot,
'
',
rna_plot,
'
',
cluster_infor,
'
',
paste0('const simplified_elements= ["',
       paste0(simplified_content, collapse = '","'),
       '"]'),
'
',
paste0("const shrink_ratio_x =", 180 / max_x),
'
',
paste0("const shrink_ratio_y =", 200 / max_y),
'
',
paste0('const colour_platte = ["',
       paste0(colour_palette, collapse = '","'),
       '"]'),
'
// First layer cluster Second layer different pathways
',
network,
'

const colorScale = d3.scaleLinear()
.domain([-',
scale_legend,
', 0, ',
scale_legend,
'])  // Define the input domain
.range(["blue", "yellow", "red"]);
let selectedNetwork1 = 0

  var string_sel = document.getElementById("select1").selectedOptions[0].text
  var rasterCanvas = document.getElementById("rc2");
  var ctx = rasterCanvas.getContext("2d");
  ctx.clearRect(0, 0, rasterCanvas.width, rasterCanvas.height);
    coordinates.forEach((coord, index) => {
      var temp_clu = cluster_info[index];
      if(temp_clu == string_sel.slice(7)){
        var color = "red";
      }else{
        var color = "gray";
      }
      ctx.fillStyle = color;
      ctx.fillRect(coord[0], coord[1], shrink_ratio_x, shrink_ratio_y); // Draw a small square (10x10)
    });
    document.getElementById("clu_window").textContent = string_sel;

document.getElementById("select1").addEventListener("change", function (event) {
  selectedNetwork1 = event.target.selectedIndex;
  switchNetwork(network_ind, selectedNetwork1, upper = 1);
  console.log(selectedNetwork1);
  var rasterCanvas = document.getElementById("rc2");
  var ctx = rasterCanvas.getContext("2d");
  string_sel = document.getElementById("select1").selectedOptions[0].text
  ctx.clearRect(0, 0, rasterCanvas.width, rasterCanvas.height);
    coordinates.forEach((coord, index) => {
      var temp_clu = cluster_info[index];
      if(temp_clu == string_sel.slice(7)){
        var color = "red";
      }else{
        var color = "gray";
      }
      ctx.fillStyle = color;
      ctx.fillRect(coord[0], coord[1], shrink_ratio_x, shrink_ratio_y); // Draw a small square (10x10)
    });
    document.getElementById("clu_window").textContent = string_sel;
});
//var svg = d3.select("#network");
const networkContainer = document.getElementById("network-frame");
function checkBounds(d) {
  if (d.x < 0) d.x = 10;
  if (d.x > width) d.x = width - 30;
  if (d.y < 0) d.y = 30;
  if (d.y > height) d.y = height - 10;
}
// Initial dimensions
let width = networkContainer.clientWidth;
let height = networkContainer.clientHeight;


window.addEventListener("resize", () => {
  // Update dimensions on resize
  width = networkContainer.clientWidth;
  height = networkContainer.clientHeight;
  switchNetwork(network_ind, selectedNetwork1, 1)
  //console.log(`New size: ${width}px by ${height}px`);
});

var checked = 0

function updateNetwork(index, partition, upper,full) {
  const svg = d3.select("#network");
  svg.selectAll("*").remove(); // Clear the current network

  svg.attr("width", width)
  .attr("height", height);
  let network_full = networks[partition][index];
  let test = networks[0][0].links.filter(link => simplified_elements.includes(link.source.id) | simplified_elements.includes(link.target.id))
  console.log(test)
  if(full){
    var network = network_full;
  }else{
    if(test.length == 0){
      let filteredLinks = network_full.links.filter(link =>
     simplified_elements.includes(link.source) | simplified_elements.includes(link.target));
     console.log(network_full.links)
     console.log(simplified_elements)
     console.log(filteredLinks)
     let tobeadded = filteredLinks.flatMap(link => [link.source, link.target])
console.log(tobeadded)
     var filteredNodes = network_full.nodes.filter(b => tobeadded.includes(b.id));
    var network = {
  nodes: filteredNodes,
  links: filteredLinks};
    }else{
      let filteredLinks = network_full.links.filter(link =>
     simplified_elements.includes(link.source.id) | simplified_elements.includes(link.target.id));
     console.log(network_full.links)
     console.log(simplified_elements)
     console.log(filteredLinks)
     let tobeadded = filteredLinks.flatMap(link => [link.source.id, link.target.id])
console.log(tobeadded)
     var filteredNodes = network_full.nodes.filter(b => tobeadded.includes(b.id));
    var network = {
  nodes: filteredNodes,
  links: filteredLinks};
    }
  };
console.log(network)
  var container = svg.append("g")
  .attr("width", width)
  .attr("height", height);

  container.append("defs").selectAll("marker")
  .data(["end"])
  .enter().append("marker")
  .attr("id", "end")
  .attr("viewBox", "0 -5 10 10")
  .attr("refX", 10)
  .attr("refY", 0)
  .attr("markerWidth", 6)
  .attr("markerHeight", 6)
  .attr("orient", "auto")
  .append("path")
  .attr("d", "M0,-5L10,0L0,5")
  .attr("fill", "#000");



  var scale = 1;

  // Define arrowheads

  network.links.forEach((link, i) => {
 if(link.head == "arrow"){
    container.append("defs").append("marker")
    .attr("viewBox", "0 -5 10 10")
    .attr("id", `end-arrow${i}`)
    .attr("refX", 13)
    .attr("refY", 0)
    .attr("markerWidth", 5)
    .attr("markerHeight", 5)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M0,-5L10,0L0,5")
    .attr("fill", link.type);  // Set the arrowhead color to match the stroke
    }
    if(link.head == "crosshead"){
    container.append("defs").append("marker")
    .attr("id", `end-arrow${i}`)
    .attr("refX", 13)
    .attr("refY", 3.5)
    .attr("markerWidth", 7)
    .attr("markerHeight", 7)
    .attr("orient", "auto")
    var marker_sel = d3.select(`#end-arrow${i}`);
    marker_sel.append("line")
    .attr("class", "marker-line")
    .attr("x1","0")
    .attr("y1","0")
    .attr("x2","7")
    .attr("y2","7")
    .attr("stroke", link.type)
    .attr("stroke-width", "1")


    marker_sel.append("line")
    .attr("class", "marker-line")
    .attr("x1","7")
    .attr("y1","0")
    .attr("x2","0")
    .attr("y2","7")
    .attr("stroke", link.type)
    .attr("stroke-width", "1");  // Set the arrowhead color to match the stroke
    }

    if(link.head == "revarrow"){
    container.append("defs").append("marker")
    .attr("id", `end-arrow${i}`)
    .attr("refX", 14)
    .attr("refY", 3.5)
    .attr("markerWidth", 7)
    .attr("markerHeight", 7)
    .attr("orient", "auto")
    var marker_sel = d3.select(`#end-arrow${i}`);
    marker_sel.append("line")
    .attr("class", "marker-line")
    .attr("x1","3.5")
    .attr("y1","3.5")
    .attr("x2","7")
    .attr("y2","7")
    .attr("stroke", link.type)
    .attr("stroke-width", "1")
    //.attr("viewBox", "0 -5 10 10")
    marker_sel.append("line")
    .attr("class", "marker-line")
    .attr("x1","3.5")
    .attr("y1","3.5")
    .attr("x2","7")
    .attr("y2","0")
    .attr("stroke", link.type)
    .attr("stroke-width", "1");  // Set the arrowhead color to match the stroke
    }
    if(link.head == "thead"){
    container.append("defs").append("marker")
    .attr("id", `end-arrow${i}`)
    //.attr("viewBox", "0 -5 10 10")
    .attr("refX", 10)
    .attr("refY", 5)
    .attr("markerWidth", 3)
    .attr("markerHeight", 10)
    .attr("orient", "auto")
    .append("rect")
    .attr("class", "marker-rect")
    .attr("width","10")
    .attr("height","10")
    .attr("fill", link.type);
    }

    if(link.head == "round"){
    container.append("defs").append("marker")
    .attr("id", `end-arrow${i}`)
    //.attr("viewBox", "0 -5 10 10")
    .attr("refX", 9)
    .attr("refY", 5)
    .attr("markerWidth", 10)
    .attr("markerHeight", 10)
    .attr("orient", "auto")
    .append("circle")
    .attr("r", "2")
    .attr("cx","5")
    .attr("cy","5")
    .attr("stroke", link.type)
    .attr("fill", "transparent")
    .attr("stroke-width", "1");
    }
  });


  const simulation = d3.forceSimulation(network.nodes)
  .force("link", d3.forceLink(network.links).id(d => d.id).distance(d =>Math.log10(d.weight)*60))
  .force("charge", d3.forceManyBody().strength(-30))
  .force("center", d3.forceCenter(width / 2, height / 2))
  .force("collide", d3.forceCollide(30).strength(0.3))
  .on("tick", ticked);
  const drag = d3
  .drag()
  .on("start", dragstart)
  .on("drag", draggeddrag);

  function dragstart() {
    d3.select(this).classed("fixed", true);
  }

  function draggeddrag(event, d) {
    d.fx = clamp(event.x, 0, width);
    d.fy = clamp(event.y, 0, height);
    simulation.alpha(1).restart();
  }

  var link = container.append("g")
  .attr("class", "links")
  .selectAll("line")
  .data(network.links)
  .enter().append("line")
  .attr("stroke", d => d.type)
  .attr("stroke-dasharray", d => d.style === "dashed" ? "4 2" : null)
  .attr("marker-end", (d, i) => `url(#end-arrow${i})`) // Use the dynamic marker
  //.attr("marker-end", "url(#end)")
  .attr("stroke-width", 1);

  var node = container
  .selectAll("g.node")
  .data(network.nodes)
  .enter().append("svg:g")
  .attr("class", function (d) {
    if (d.shape === "mets") {
      return "rect node";
    } else if (d.shape === "protein") {
      return "hex node";
    } else if (d.shape === "rna") {
      return "circle node";
    }
  })
  .attr("transform", "translate(0, 0)")
  .attr("fill", d => typeof d.expr === "number" ? colorScale(d.expr) : "gray")
  .call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended));

  function click(event, d) {
    delete d.fx;
    delete d.fy;
    d3.select(this).classed("fixed", false);
    simulation.alpha(1).restart();
  }
  function clamp(x, lo, hi) {
    return x < lo ? lo : x > hi ? hi : x;
  }

  node.call(drag).on("click", click);

  d3.selectAll(".rect").append("rect")
  .data(network.nodes)
  .attr("width", d => 10 )
  .attr("height", d => 10 );

  d3.selectAll(".circle").append("circle")
  .data(network.nodes)
  .attr("r", d => 6);


  const hexagonPoints = (radius) => {
    const halfWidth = radius * Math.sqrt(3) / 2;
    return `
    0,${-radius}
    ${halfWidth},${-radius / 2}
    ${halfWidth},${radius / 2}
    0,${radius}
    ${-halfWidth},${radius / 2}
    ${-halfWidth},${-radius / 2}`;
  };
  d3.selectAll(".hex").append("polygon")
  .attr("points", function (d) {
    return hexagonPoints(6 )
  });


  node.append("title")
  .text(d => d.name);


  let label = container.append("g")
  .attr("class", "labels")
  .selectAll("text")
  .data(network.nodes)
  .enter().append("text")
  .attr("class", "label")
  .text(d => d.name).attr("font-size", "12px")
  .attr("fill", d => typeof d.expr === "number" ? colorScale(d.expr) : "gray");

  function dragstarted(event, d) {
    delete d.fx;
    delete d.fy;
    d3.select(this).classed("fixed", false);
    simulation.alpha(0.2).restart();
  }

  function dragged(event, d) {
    d3.select(this).attr("cx", d.fx).attr("cy", d.fy);
    d.fx = event.x;
    d.fy = event.y;
  }

  function dragended(event, d) {
    d.fx = null;
    d.fy = null;
  }

  function ticked() {

    link.attr("x1", function (d) {
      checkBounds(d.source);
      return d.source.x;
    })
    .attr("y1", function (d) {
      checkBounds(d.source);
      return d.source.y;
    })
    .attr("x2", function (d) {
      checkBounds(d.target);
      return d.target.x;
    })
    .attr("y2", function (d) {
      checkBounds(d.source);
      return d.target.y;
    }).attr("d", function (d) {

      var x1 = d.source.x,
      y1 = d.source.y,
      x2 = d.target.x,
      y2 = d.target.y,
      dx = x2 - x1,
      dy = y2 - y1,
      dr = Math.sqrt(dx * dx + dy * dy),
      drx = dr,
      dry = dr,
      xRotation = 0,
      largeArc = 0,
      sweep = 1;
      if (x1 === x2 && y1 === y2) {
        xRotation = -45;
        largeArc = 1;
        drx = 30;
        dry = 20;
        x2 = x2 + 1;
        y2 = y2 + 1;
      }
      return "M" + x1 + "," + y1 + "A" + drx + "," + dry + " " + xRotation + "," + largeArc + "," + sweep + " " + x2 + "," + y2;
    });

    node
    .attr("transform", function (d) {
      checkBounds(d);
      return "translate(" + d.x + ", " + d.y + ")";
    });


    label.attr("x", function (d) {
      checkBounds(d);
      return d.x;
    })
    .attr("y", function (d) {
      checkBounds(d);
      return d.y - 10;
    });
  }

  d3.select("#node-scale").on("input", function () {
    scale = +this.value;
    d3.select("#network").selectAll("*")
    .attr("r", 6 * scale)
    .attr("width", 10 * scale)
    .attr("height", 10 * scale)
    .attr("points", function (d) {
      return hexagonPoints(6*scale)
    });
    d3.select("#network2").selectAll("*")
    .attr("r", 6 * scale)
    .attr("width", 10 * scale)
    .attr("height", 10 * scale)
    .attr("points", function (d) {
      return hexagonPoints(6*scale)
    });
    d3.select("#network").selectAll("text")
    .attr("font-size", 12*scale + "px");
    d3.select("#network2").selectAll("text")
    .attr("font-size", 12*scale + "px");

  });


  d3.select("#link-scale").on("input", function () {
    scale = +this.value;
    d3.select("#network").selectAll("*")
    .attr("stroke-width", 1 * scale);
    d3.select("#network2").selectAll("*")
    .attr("stroke-width", 1 * scale);
  });
  node.on("click", function(event, d) {
    console.log(d.id)
    var rasterCanvas = document.getElementById("rc1");
    var ctx = rasterCanvas.getContext("2d");
    var legendCanvas = document.getElementById("lca1");
    var legendCtx = legendCanvas.getContext("2d");
    var legendHeight = legendCanvas.height-10;
    var numSteps = 100; // Number of steps in the gradient
    var numTicks = 6;
    var tickSpacing = legendHeight / (numTicks - 1);
    legendCtx.fillStyle = "black";
    legendCtx.font = "6px Arial";
    legendCtx.textAlign = "right";
    legendCtx.textBaseline = "right";
    ctx.clearRect(0, 0, rasterCanvas.width, rasterCanvas.height);
    legendCtx.clearRect(0, 0, legendCanvas.width, legendCanvas.height);
    console.log("Node clicked:", d);
    if(d.group === "rna"){
      var rna_plot_sub =  rnaplot[d.id]
      document.getElementById("spatial_window").textContent = d.display;
      // Normalize values to range 0-1
      var minValue = 0;
      var maxValue = getPercentileFromObject(rna_plot_sub ,1);

      coordinates.forEach((coord, index) => {
        const normalizedValue = normalize(rna_plot_sub[index],minValue,maxValue);
        const color =colour_platte[Math.floor(normalizedValue*(colour_platte.length))];
        ctx.fillStyle = color;
        ctx.fillRect(coord[0], coord[1], shrink_ratio_x, shrink_ratio_y); // Draw a small square (10x10)
      });
      // Draw color gradient on the legend
      for (let i = 0; i <= numSteps; i++) {
        const normalizedValue = i / numSteps;
        const color = colour_platte[Math.floor(i*(colour_platte.length)/numSteps)]
        legendCtx.fillStyle = color;
        legendCtx.fillRect(28, legendHeight - (i * (legendHeight / numSteps)),28 , (legendHeight / numSteps));
      }

      for (let i = 0; i < numTicks; i++) {
        const yPosition = legendHeight +5 - (i * tickSpacing);
        const valueAtTick = minValue + (i / (numTicks - 1)) * (maxValue - minValue);

        // Draw tick mark
        legendCtx.beginPath();
        legendCtx.moveTo(10, yPosition);
        legendCtx.lineTo(10, yPosition);
        legendCtx.stroke();

        // Draw label
        legendCtx.fillText(valueAtTick.toFixed(1),27.5, yPosition);
      }
    }else{
      var met_plot_sub =  metplot[d.id]
      document.getElementById("spatial_window").textContent = d.display;
      // Normalize values to range 0-1
      var minValue = 0;
      var maxValue = getPercentileFromObject(met_plot_sub,1);

      coordinates.forEach((coord, index) => {
        const normalizedValue = normalize(met_plot_sub[index],minValue,maxValue);
        const color = colour_platte[Math.floor(normalizedValue*(colour_platte.length))];
        ctx.fillStyle = color;
        ctx.fillRect(coord[0], coord[1], shrink_ratio_x, shrink_ratio_y); // Draw a small square (10x10)
      });

      // Draw color gradient on the legend
      for (let i = 0; i <= numSteps; i++) {
        const normalizedValue = i / numSteps;
        const color = colour_platte[Math.floor(i*(colour_platte.length)/numSteps)];
        legendCtx.fillStyle = color;
        legendCtx.fillRect(28, legendHeight - (i * (legendHeight / numSteps)),28 , (legendHeight / numSteps));
      }

      for (let i = 0; i < numTicks; i++) {
        const yPosition = legendHeight +5 - (i * tickSpacing);
        const valueAtTick = minValue + (i / (numTicks - 1)) * (maxValue - minValue);

        // Draw tick mark
        legendCtx.beginPath();
        legendCtx.moveTo(10, yPosition);
        legendCtx.lineTo(10, yPosition);
        legendCtx.stroke();

        // Draw label
        legendCtx.fillText(valueAtTick.toFixed(1),27.5, yPosition);
      }
    }
  });
  //d3.select("#edge-scale").on("input", function() {
    //    var edgeScale = +this.value;
    //    link
    //        .attr("stroke-width", 2 * edgeScale);
    //});
}
var network_ind = 0

function switchNetwork(index, partition, upper,full) {
  document.getElementById("node-scale").value = 1;
  updateNetwork(index, partition, upper,full);
  network_ind = index
}
const slider = document.getElementById("slider");
document.querySelectorAll(".tab").forEach(tab => {
  tab.addEventListener("click", function () {
    document.querySelectorAll(".tab").forEach(t => t.classList.remove("active"));
    this.classList.add("active");
    const index = Array.from(document.querySelectorAll(".tab")).indexOf(this);
    switchNetwork(index, selectedNetwork1, 0,slider.textContent==="Simplified net"?false:true);
    console.log(index)
  });
});
// Initialize the first network
switchNetwork(0, selectedNetwork1, 0, false);

function normalize(value, minValue, maxValue){
  return (value - minValue) / (maxValue - minValue);
}

// Function to map normalized values to gradient colors (from blue to red)
function getGradientColor(value) {
  const r = Math.floor(value * 255); // Red increases with value
  const g = 0;                       // No green for simplicity
  const b = Math.floor((1 - value) * 255); // Blue decreases with value
  return `rgb(${r},${g},${b})`;
}


function flattenArray(arr) {
  let result = [];

  arr.forEach(item => {
    if (Array.isArray(item)) {
      result = result.concat(flattenArray(item)); // Recursively flatten
    } else {
      result.push(item); // Add non-array values
    }
  });

  return result;
}

function getPercentileFromObject(obj,x) {
  let allValues = [];

  // Iterate through each key in the object
  Object.keys(obj).forEach(key => {
    allValues = allValues.concat(obj[key]);
  });

  // Filter out non-numeric or invalid values (optional)
  const validValues = allValues.filter(value => typeof value === "number" && !isNaN(value));

  validValues.sort((a, b) => a - b);
  const index = Math.ceil(x * validValues.length) - 1; // Subtract 1 for 0-based index
  return validValues[index];
}
function click_drag(event, d) {
  delete d.fx;
  delete d.fy;
  d3.select(this).classed("fixed", false);
  simulation.alpha(1).restart();
}
function toggleState() {
  const slider = document.getElementById("slider");

  slider.classList.toggle("active");

  if (slider.classList.contains("active")) {
    slider.innerHTML = "Full net";
    switchNetwork(network_ind, selectedNetwork1,0,true);
    checked = 1;
  } else {
    slider.innerHTML = "Simplified net";
    switchNetwork(network_ind, selectedNetwork1,0,false);
    checked = 1;
  }
}
</script>
</body>
</html>')
  
  full_path <- paste0(path, "/", returnname, ".html")
  if (file.access(path, mode = 2) != 0)
  {
    stop("Error: Cannot write to the specified path. Access denied.")
  }
  if (file.exists(full_path))
  {
    warning(paste(
      "Warning: File",
      full_path,
      "already exists and will be overwritten."
    ))
  }
  tryCatch({
    writeLines(html, full_path)
    message(paste("File successfully written to", full_path))
  }, error = function(e) {
    stop(paste("Error: Failed to write the file. Details:", e$message))
  })
}

