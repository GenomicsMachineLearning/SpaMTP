#' Constructs an interactive network for exploring spatial metabolomics and transcriptomics data.
#'
#' @param SpaMTP A `SpaMTP` Seurat object containing spatial metabolomics (SM) and/or spatial transcriptomics (ST) data. If SM data is included, it must be annotated using the `SpaMTP::AnnotateSM()` function.
#' @param ident A character string specifying the cluster identifier used to group regions, corresponding to a column name in the `SpaMTP@meta.data` slot.
#' @param DE.list A list containing differential expression results from the `FindAllMarkers()` function, with items matching the order of the `analyte_types` argument.
#' @param regpathway A dataframe output from the `SpaMTP::FindRegionalPathways()` function, containing identified regional pathways.
#' @param selected_pathways A character vector specifying the names or IDs of pathways used to construct the network (e.g., `c("Amino acid metabolism", "WP1902", "Aspartate and asparagine metabolism")`). This argument is not case-sensitive.
#' @param path The directory to save the output. If not provided, the default is the current working directory.
#' @param SM_assay A character string specifying the assay name for spatial metabolomics data in `SpaMTP` to extract intensity values (default: `"SPM"`).
#' @param ST_assay A character string specifying the assay name for spatial transcriptomics data in `SpaMTP` to extract RNA count values (default: `"SPT"`).
#' @param SM_slot The slot name containing the spatial metabolomics assay matrix (default: `"counts"`).
#' @param ST_slot The slot name containing the spatial transcriptomics assay matrix (default: `"counts"`).
#' @param colour_palette The color palette used to plot the spatial image in the output HTML file. Default: `grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)`.
#' @param analyte_types A subset of `c("genes", "metabolites")`. Can be `c("genes")`, `c("metabolites")`, or both.
#' @param image Character string specifying which image stored within the SpaMTP object to use for plotting (default = "slice1").
#' @param verbose A logical value indicating whether to display detailed messages during execution (default: `FALSE`).
#'
#' @return An interactive HTML file visualizing the network structure of the specified pathways.
#' @export
#'
#' @examples
#' #PathwayNetworkPlots(SpaMTP, ident = "Custom_ident", regpathway = regpathway, DE.list = DE.list, selected_pathways = selected_pathways)
PathwayNetworkPlots  = function(SpaMTP,
                                ident,
                                regpathway,
                                DE.list,
                                selected_pathways = NULL,
                                path = getwd(),
                                SM_slot = "counts",
                                ST_slot = "counts",
                                colour_palette = NULL,
                                SM_assay = "SPM",
                                ST_assay = "SPT",
                                analyte_types = c("genes", "metabolites"),
                                image = "slice1",
                                verbose = T) {
  if ("genes" %in% analyte_types) {
    if (is.null(SpaMTP@assays[[ST_assay]])) {
      stop(
        paste0(
          "If you are use genetice data with 'gene' in 'analyte_types' (e.g spatial transcriptolomics), please ensure that '",
          ST_assay,
          "' is inside the input seurat object"
        )
      )
    } else{
      gene_matrix = Matrix::t(SpaMTP[[ST_assay]]@layers[[ST_slot]])
    }
  }
  if ("metabolites" %in% analyte_types) {
    if (is.null(SpaMTP@assays[[SM_assay]])) {
      stop(
        paste0(
          "If you are use metabolic data with 'metabolites' in 'analyte_types' (e.g spatial metabolomics), please ensure that '",
          SM_assay,
          "' is inside the input seurat object"
        )
      )
    } else{
      mass_matrix = Matrix::t(SpaMTP[[SM_assay]]@layers[[SM_slot]])
    }
  }
  if (("metabolites" %in% analyte_types) &
      ("genes" %in% analyte_types)) {
    if (ncol(SpaMTP@assays[[SM_assay]]) != ncol(SpaMTP@assays[[ST_assay]])) {
      stop("Please align the spatial data before processing")
    }
  }
  assignment = SpaMTP@meta.data[[ident]]

  # get enrichment dataframe

  SpatialColors <- grDevices::colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
  colour_palette = colour_palette %||% SpatialColors(100)


  importance = regpathway %>% dplyr::group_by(pathwayName) %>%
    dplyr::mutate(group_importance = sum(abs(NES))) %>%
    filter(!duplicated(pathwayName)) %>% arrange(desc(group_importance))
  sub_enriched = regpathway[which(
    tolower(regpathway$pathwayName) %in% tolower(selected_pathways %||% importance$pathwayName[1:10]) |
      tolower(regpathway$sourceId) %in% tolower(selected_pathways %||% importance$sourceId[1:10])
  ), ]
  type = sub_enriched$type
  check_topology_existence = function(type,
                                      pathway_of_interest,
                                      id_of_interest) {
    temp_db = c()
    if (type == "wiki") {
      wikiids = unlist(lapply(RAMP_wikipathway, function(x) {
        return(x[["id"]])
      }))
      index = which((tolower(names(
        RAMP_wikipathway
      )) == tolower(pathway_of_interest)) |
        (wikiids == tolower(id_of_interest)))
      if (length(index) != 0) {
        temp_db = RAMP_wikipathway
      }
    } else if (type == "reactome") {
      reactomeids = unlist(lapply(RAMP_Reactome, function(x) {
        return(x[["id"]])
      }))
      index = which((tolower(names(
        RAMP_Reactome
      )) == tolower(pathway_of_interest)) |
        (reactomeids == tolower(id_of_interest)))
      if (length(index) != 0) {
        temp_db = RAMP_Reactome
      }
    } else if (type == "kegg") {
      keggids = unlist(lapply(RAMP_kegg, function(x) {
        return(x[["id"]])
      }))
      index = which((tolower(names(
        RAMP_kegg
      )) == tolower(pathway_of_interest)) |
        (keggids == tolower(id_of_interest)))
      if (length(index) != 0) {
        temp_db = RAMP_kegg
      }
    } else if (type == "hmdb") {
      hmdbids = unlist(lapply(RAMP_hmdb, function(x) {
        return(x[["id"]])
      }))
      index = which((tolower(names(
        RAMP_hmdb
      )) == tolower(pathway_of_interest)) |
        (tolower(hmdbids) == tolower(id_of_interest)))
      if (length(index) != 0) {
        temp_db = RAMP_hmdb
      }
    } else{
      all_list = c(RAMP_wikipathway,
                   RAMP_Reactome,
                   RAMP_kegg,
                   RAMP_hmdb)
      all_names = names(all_list)
      add_ids = unlist(lapply(all_list, function(x) {
        return(x[["id"]])
      }))
      index = which((tolower(all_names) == tolower(pathway_of_interest)) |
                      (tolower(add_ids) == tolower(id_of_interest)))
      if (length(index) == 0) {
        return(NULL)
      } else{
        return(all_list[index[1]])
      }
    }
    if (exists("temp_db")) {
      return(temp_db[[index[1]]])
    } else{
      return(NULL)
    }
  }

  topodb = apply(
    sub_enriched[!duplicated(sub_enriched$pathwayName), ],
    MARGIN = 1,
    FUN = function(x) {
      check_topology_existence(x$type, x$pathwayName, x$sourceId)
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

  #Get DE list
  DE = list()
  for (i in 1:length(analyte_types)) {
    verbose_message(
      message_text = paste0(
        "Assuming DE.list[",
        i,
        "] contains ",
        analyte_types[i] ,
        " results .... "
      ),
      verbose = verbose
    )
    DE.list[[analyte_types[i]]] <- DE.list[[i]]
  }

  if (is.null(SpaMTP@tools[["db_3"]])) {
    verbose_message(message_text = "Please use MZAnnotation " , verbose = verbose)
  }

  db_3 <- SpaMTP@tools$db_3
  db_3 = db_3 %>%
    tidyr::separate_rows(Isomers, sep = ";")
  verbose_message(message_text = "Query necessary data and establish pathway database" , verbose = verbose)
  input_id = lapply(db_3$Isomers, function(x) {
    x = unlist(x)
    index_hmdb = which(grepl(x, pattern = "HMDB"))
    x[index_hmdb] = paste0("hmdb:", x[index_hmdb])
    index_chebi = which(grepl(x, pattern = "CHEBI"))
    x[index_chebi] = tolower(x[index_chebi])
    index_lm = which(grepl(x, pattern = "LMPK"))
    x[index_lm] = tolower(x[index_lm])
    return(x)
  })
  db_3 = db_3 %>% dplyr::mutate(inputid = input_id) %>%  dplyr::mutate(chem_source_id = input_id)
  rampid = c()
  verbose_message(message_text = "Query db for addtional matching" , verbose = verbose)
  db_3 = merge(chem_props, db_3, by = "chem_source_id")
  ### Adding DE Results
  db_3 = db_3 %>% mutate(mz_name = paste0("mz-", db_3$observed_mz))

  if (length(DE.list) != length(analyte_types)) {
    stop(
      "Number of DE data.frames provided does not match the number of analyte types specified. Please make sure a DE dataframe is provided for each analyte type"
    )
  }
  verbose_message(message_text = "Constructing DE dataframes.... ", verbose = verbose)
  for (i in 1:length(analyte_types)) {
    DE <- DE.list[[i]]
    if (any(c("avg_log2FC", "logFC") %in% colnames(DE)) &&
        any(c("p_val_adj", "FDR") %in% colnames(DE)) &&
        "cluster" %in% colnames(DE) &&
        "gene" %in% colnames(DE)) {
      if ("logFC" %in% colnames(DE)) {
        colnames(DE)[colnames(DE) == "logFC"] <- "avg_log2FC"
      }
      # Rename FDR to p_val_adj if FDR exists
      if ("FDR" %in% colnames(DE)) {
        colnames(DE)[colnames(DE) == "FDR"] <- "p_val_adj"
      }
      if (analyte_types[i] == "metabolites") {
        DE = DE %>% rename(mz_name = gene)
        db_3 = merge(db_3 , DE, by = "mz_name")
        DE.list[[analyte_types[i]]] <- db_3
      } else {
        DE = DE %>% mutate(commonName = toupper(gene))
        source_gene = merge(DE, source_df[which(grepl(source_df$rampId, pattern = "RAMP_G")), ], by = "commonName")
        DE.list[[analyte_types[i]]] <- source_gene
      }
    } else {
      stop(
        "DE dataframe [",
        i,
        "] provided does not have the correct column names ... column names MUST include 'cluster', 'gene', ('avg_log2FC' or 'logFC') and ('p_val_adj' or 'FDR'). Please adjust column names in all DE data.frames to match ..."
      )
    }
  }

  sub_enriched = sub_enriched[which(tolower(sub_enriched$pathwayName) %in% tolower(names(topodb))), ]
  pathway_names = unique(sub_enriched$pathwayName)

  network = paste0('const networks = [')
  matrix_ids = c()
  ucid = naturalsort::naturalsort(unique(sub_enriched$Cluster_id))

  simplified_content = c()
  for (i in 1:length(ucid)) {
    sub_cluster = sub_enriched[which(sub_enriched$Cluster_id == ucid[i]), ]
    if ("genes" %in% analyte_types) {
      sub_expr_rna = DE.list[["genes"]][which(tolower(DE.list[["genes"]]$cluster) ==  tolower(ucid[i])), ]
    }
    if ("metabolites" %in% analyte_types) {
      sub_expr_met = DE.list[["metabolites"]][which(tolower(DE.list[["metabolites"]]$cluster) ==  tolower(ucid[i])), ]
    }
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
            if (grepl(pathway_analyte[k], pattern = "RAMP_G")) {
              rna_expr = sub_expr_rna[which(sub_expr_rna$rampId == pathway_analyte[k]), ]
              rna_expr = rna_expr[which.min(rna_expr$p_val_adj), ]
            }
            if (grepl(pathway_analyte[k], pattern = "RAMP_C")) {
              met_expr = sub_expr_met[which(sub_expr_met$ramp_id == pathway_analyte[k]), ]
              met_expr = met_expr[which.min(met_expr$p_val_adj), ]
            }
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
                ifelse((
                  length(rna_expr$avg_log2FC) == 0
                ), '"NA"', rna_expr$avg_log2FC),
                ifelse((
                  length(met_expr$avg_log2FC) == 0
                ), '"NA"', met_expr$avg_log2FC)
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
                ifelse((length(
                  rna_expr$gene
                ) == 0), 'NA', rna_expr$gene),
                ifelse((
                  length(met_expr$common_name) == 0
                ), 'NA', met_expr$common_name)
              ),
              '", size:',
              ifelse(
                grepl(pathway_analyte[k], pattern = "RAMP_G"),
                ifelse((
                  length(rna_expr$avg_log2FC) == 0
                ), '"NA"', rna_expr$avg_log2FC),
                ifelse((
                  length(met_expr$avg_log2FC) == 0
                ), '"NA"', met_expr$avg_log2FC)
              ),
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
  tab_div = paste0(tab_div,
                   '<div class="tab" id = "default_tab">',
                   pathway_names[1],
                   "</div>")
  for (o in 2:length(pathway_names)) {
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

  fc_vector = c()
  if ("genes" %in% analyte_types) {
    fc_vector = c(fc_vector, DE.list[["genes"]]$avg_log2FC)
  }
  if ("metabolites" %in% analyte_types) {
    fc_vector = c(fc_vector, DE.list[["metabolites"]]$avg_log2FC)
  }

  scale_legend = as.integer(sqrt(max(abs(fc_vector))))

  # Get coordinates
  coordnate = Seurat::GetTissueCoordinates(SpaMTP, image = image)
  non_na_ind = which((!is.na(coordnate[, "x"])) &
                       (!is.na(coordnate[, "y"])))
  coordnate = cbind(coordnate, assign = as.character(assignment))
  coordnate = na.omit(coordnate)
  max_x =  max(na.omit(as.numeric(coordnate[, "x"])))
  max_y =  max(na.omit(as.numeric(coordnate[, "y"])))


  coordi = paste0("const coordinates = [")

  for (t in 1:nrow(coordnate)) {
    coordi = paste0(
      coordi,
      '[',
      as.numeric(coordnate[t, "x"]) * 160 / max_x ,
      ',',
      as.numeric(coordnate[t, "y"]) * 180 / max_y,
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
        met_temp_sub  = DE.list[["metabolites"]][which(DE.list[["metabolites"]]$ramp_id == matrix_ids[a]), ]
        mz_name = met_temp_sub$mz_name[which.min(met_temp_sub$p_val_adj)]
        met_mat_index = which(tolower(rownames(SpaMTP@assays[[SM_assay]]@features)) ==  mz_name)
        mmat_return = mass_matrix[non_na_ind, met_mat_index]
        met_plot = paste0(met_plot,
                          '"',
                          matrix_ids[a],
                          '":[',
                          paste0(mmat_return, collapse = ","),
                          '],')
      }
    } else{
      if ("genes" %in% analyte_types) {
        rna_temp_sub = DE.list[["genes"]]$commonName[which(DE.list[["genes"]]$rampId == matrix_ids[a])]
        gene_mat_index = which(tolower(rownames(SpaMTP@assays[[ST_assay]]@features)) == tolower(unique(rna_temp_sub)[1]))
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
                         paste0(coordnate[, "assign"], collapse = '","'),
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
      height: 92vh;
      position: flex;
    }

    .legend-rectangle {
      width: 90%;
      height: 20px;
      background: linear-gradient(to right, blue, yellow, red);
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
    .parent {
      display: flex; /* Align children horizontally */
      gap: 10px; /* Add spacing between elements */
      width: 100%;
    }
    .shape-text {
      font-size: 12px;
    }

    .scrollable-select {
      width: 200px;
      height: 50px;
      /* Height determines the visible area */
      overflow-y: auto;
      /* Add vertical scroll */
    }
  </style>
  </head>
  <body>
  <div id="main_all" class = "parent">
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
      <button id="saveButton">Save network as svg Image</button>
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
      <canvas id="rc1" class="rasterCanvas" width="160" height="180"></canvas>
      <canvas id="lca1" class="legendCanvas" width="40" height="180"></canvas>
    </div>
    <p width="200" display="inline-block" id="spatial_window">Display spatial data by clicking on vertices</p>
    <div class="rastercontainer">
      <canvas id="rc2" class="rasterCanvas" width="160" height="180"></canvas>
    </div>
    <p width="200" display="inline-block" id="clu_window">Cluster 1</p>
  </div>
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
paste0(
  'const simplified_elements= ["',
  paste0(simplified_content, collapse = '","'),
  '"]'
),
'
',
paste0("const shrink_ratio_x =", 1),
'
',
paste0("const shrink_ratio_y =", 1),
'
',
paste0(
  'const colour_platte = ["',
  paste0(colour_palette, collapse = '","'),
  '"]'
),
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
      if(temp_clu == string_sel.slice(7)|(string_sel == "Baseline" & temp_clu==string_sel)|temp_clu == string_sel){
        var color = "red";
      }else{
        var color = "gray";
      }
      ctx.fillStyle = color;
      ctx.fillRect(coord[0], coord[1], shrink_ratio_x, shrink_ratio_y); // Draw a small square (10x10)
    });
    document.getElementById("clu_window").textContent = string_sel;

document.getElementById("select1").addEventListener("change", function (event) {
  if(slider.textContent==="Simplified net"){
    toggleState()
  };
  selectedNetwork1 = event.target.selectedIndex;
  switchNetwork(network_ind, selectedNetwork1, upper = 1, slider.textContent==="Simplified net"?false:true);
  console.log(selectedNetwork1);
  var rasterCanvas = document.getElementById("rc2");
  var ctx = rasterCanvas.getContext("2d");
  string_sel = document.getElementById("select1").selectedOptions[0].text
  ctx.clearRect(0, 0, rasterCanvas.width, rasterCanvas.height);
    coordinates.forEach((coord, index) => {
      var temp_clu = cluster_info[index];
      if(temp_clu == string_sel.slice(7)|(string_sel == "Baseline" & temp_clu==string_sel)|temp_clu == string_sel){
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
default_tab.classList.add("active");
document.querySelectorAll(".tab").forEach(tab => {
  tab.addEventListener("click", function () {
    if(slider.textContent==="Simplified net"){
    toggleState()
    };
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
document.getElementById("saveButton").addEventListener("click", function () {
  const svgElement = document.getElementById("network");
  const serializer = new XMLSerializer();
  const svgString = serializer.serializeToString(svgElement);

  // Create a Blob from the SVG string
  const svgBlob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(svgBlob);

  // Create a link element to trigger the download
  const link = document.createElement("a");
  link.href = url;
  link.download = "networks.svg"; // Set the filename
  document.body.appendChild(link);
  link.click(); // Trigger the download
  document.body.removeChild(link); // Clean up
  URL.revokeObjectURL(url); // Free up memory
});
</script>
</body>
</html>')

  returnname = paste0(ident, "_",format(Sys.time(), "%Y_%m_%d_%H_%M_%S_%Z"))
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
