#' Calculates Significant Metabolic Pathways using a Fisher Exact Test
#'
#' @param Analyte A list of analytes containing a combination of three possible elements, namely "mzs", "genes" and/or "metabolites". The list must be named with these titles, corresponding to the relative input datasets. Read below for supported input formats.
#' @param min_path_size The min number of  in a specific pathway (default = 5).
#' @param max_path_size The max number of  in a specific pathway (default = 500).
#' @param alternative The hypothesis of the fisher exact test (default = "greater").
#' @param pathway_all_info Whether to included all genes/ screened in the return (default = FALSE).
#' @param pval_cutoff A numerical value defining the adjusted p value cutoff for keeing significant pathways (default = NULL).
#' @param verbose Boolean indicating whether to show informative messages. If FALSE these messages will be suppressed (default = TRUE).
#' @param ... Additional parameters that can be passed through to `annotateTable()` when running `mz`-based analysis. Please see documentation for `annotateTable()` for more details.
#'
#' ### Details
#' * Supported `metabolites` format: strings which contain the metabolite ID with database name. For example = "hmdb:HMDBX", "chebi:X", "pubchem:X","wikidata:X" ,"kegg:X" ,"CAS:X","lipidbank:X","chemspider:X","	LIPIDMAPS:X" (where X stands for upper case of the cooresponding ID in each database)
#' * Supported `genes` data format: strings which contain the gene name and formatting. For example = "entrez:X", "gene_symbol:X", "uniprot:X", "ensembl:X", "hmdb:HMDBPX"
#' * Supported `mzs` format: any string or numeric vector contains the m/z. NOTE: If `mzs` values are provided then `annotateTable()` will be run using default parameters and combining the Chebi_db, Lipidmaps_db and HMDB_db databases.
#'
#'
#' @return a dataframe with the relevant pathway information
#' @export
#'
#' @import dplyr
#' @import stringr
#'
#' @examples
#' ## Running in 'mzs' mode:
#' # FishersPathwayAnalysis(Analyte = list("mzs" = mz_values), ppm_error = 3)
#'
#' ## Running in 'metabolites' mode
#' # FishersPathwayAnalysis(Analyte = list("metabolites" = metabolite_ids))
#'
#' ## Running 'metabolites' and 'genes' combined
#' # FishersPathwayAnalysis(Analyte = list("metabolites" = metabolite_ids, "genes" = gene_names))
FishersPathwayAnalysis <- function (Analyte,
                                    max_path_size = 500,
                                    min_path_size = 5,
                                    alternative = "greater",
                                    pathway_all_info = FALSE,
                                    pval_cutoff = NULL,
                                    verbose = TRUE,
                                    ...)
{
  if (is.null(names(Analyte)) || ! all(names(Analyte) %in% c("mzs", "genes", "metabolites"))){
    stop("Invalid key argument! Name of list was not one of the required values [c('mzs', 'genes', 'metabolites')].  Please specify the names correctly for example: list('mz' = c('mz-100.12','mz-428.32', 'mz-341.201')) ... ")
  }

  verbose_message(message_text = "Running Fisher Testing ......", verbose = verbose)

  pathwayRampId <- rampId <- c()

  if ("metabolites" %in% names(Analyte)) {
    analytes_met = Analyte[["metabolites"]]
    source_met = source_df[which(grepl(source_df$rampId, pattern = "RAMP_C") == T),]
    analytehaspathway_met = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_C") == T),]
    analyte_met = analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),]
  }
  if ("genes" %in% names(Analyte)) {
    analytes_rna = Analyte[["genes"]]
    source_rna = source_df[which(grepl(source_df$rampId, pattern = "RAMP_G") == T),]
    analytehaspathway_rna = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_G") == T),]
    analyte_rna = analyte[which(grepl(analyte$rampId, pattern = "RAMP_G") == T),]
  }

  if ("mzs" %in% names(Analyte)) {

    warning("A list of mzs has been provided! annotateTable() will now be run using input specified via `...`! If no appropriate inputs have been provided default settings will be used including `db = rbind(Chebi_db, Lipidmaps_db, HMDB_db)`. Please see annotateTable() documentation for more ...")

   analytes_mz = Analyte[["mzs"]]

   input_mz = data.frame(cbind(
     row_id = 1:length(analytes_mz),
     mz = as.numeric(stringr::str_extract(analytes_mz, pattern = "\\d+\\.?\\d*"))
   ))

   rownames(input_mz) <- paste0("mz-", input_mz$mz)

   args <- list(...)

   # Set db to 2 if it's not passed in ..., otherwise use the value provided in ...
   db <- if ("db" %in% names(args)) args$db else rbind(Chebi_db,
                                                         Lipidmaps_db,
                                                         HMDB_db)

   remaining_args <- args[setdiff(names(args), "db")]

   db_3 <- do.call(annotateTable, c(list(mz_df= input_mz, db = db, verbose = verbose), remaining_args))


    db_3list = pbapply::pblapply(1:nrow(db_3), function(i){
      if (any(grepl(db_3[i, ], pattern = ";"))) {
        # Only take the first row
        ids = unlist(stringr::str_split(db_3[i, ]$Isomers, pattern = ";"))
        ids[which(grepl(ids, pattern = "HMDB"))] = paste0("HMDB:",ids[which(grepl(ids, pattern = "HMDB"))])
        ids[which(grepl(ids, pattern = "LM"))] = paste0("LIPIDMAPS:",ids[which(grepl(ids, pattern = "LM"))])
        ids = sub(" ","",ids)
        mzs = rep(db_3[i, ]$observed_mz, times = length(ids))
        adducts = rep(db_3[i, ]$Adduct, times = length(ids))
        return(cbind(ids,mzs,adducts))
      } else{
        ids = db_3[i,]$Isomers
        ids[which(grepl(ids, pattern = "HMDB"))] = paste0("HMDB:",ids[which(grepl(ids, pattern = "HMDB"))])
        ids[which(grepl(ids, pattern = "LM"))] = paste0("LIPIDMAPS:",ids[which(grepl(ids, pattern = "LM"))])
        ids = sub(" ","",ids)
        mzs = rep(db_3[i, ]$observed_mz, times = length(ids))
        adducts = rep(db_3[i, ]$Adduct, times = length(ids))
        return(cbind(ids,mzs,adducts))
      }
    })
    expand_db3_df = do.call(rbind, db_3list)
    expand_db3 = expand_db3_df[,1]
    mz_db3 = expand_db3_df[,2]
    adducts_db3 = expand_db3_df[,3]
    analytes_mz = sub(" ", "", expand_db3)
    source_mz = source_df[which(grepl(source_df$rampId, pattern = "RAMP_C") == T),]
    analytehaspathway_mz = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_C") == T),]
  }

  verbose_message(message_text = "Parsing the information of given analytes class" , verbose = verbose)

  analyte_new = analytehaspathway_new = source_new =  data.frame()
  analytes_new = mz_array = adducts_array = c()
  if("mzs" %in% names(Analyte)){
    analyte_new = rbind(analyte_new,
                        analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),])
    analytehaspathway_new = rbind(analytehaspathway_new,
                                  analytehaspathway_mz)
    source_new = rbind(source_new,
                       source_mz)
    analytes_new = c(analytes_new,
                     analytes_mz)
    mz_array = c(mz_array,mz_db3)
    adducts_array = c(adducts_array, adducts_db3)
  }

  if("metabolites" %in% names(Analyte)){
    analyte_new = rbind(analyte_new,
                        analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),])
    analytehaspathway_new = rbind(analytehaspathway_new,
                                  analytehaspathway_met)
    source_new = rbind(source_new,
                       source_met)
    analytes_new = c(analytes_new,
                     analytes_met)
    mz_array = c(mz_array,rep(NA, times = length(analytes_met)))
    adducts_array = c(adducts_array,rep(NA, times = length(analytes_met)))
  }

  if("genes" %in% names(Analyte)){
    analyte_new = rbind(analyte_new,
                        analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),])

    analytehaspathway_new = rbind(analytehaspathway_new,
                                  analytehaspathway_rna)
    source_new = rbind(source_new,
                       source_rna)
    analytes_new = c(analytes_new,
                     analytes_rna)
    mz_array = c(mz_array,rep(NA, times = length(analytes_rna)))
    adducts_array = c(adducts_array,rep(NA, times = length(analytes_rna)))
  }

  # Merge as data.frame to minimise query time
  temp_mz_analyte = data.frame(cbind(mz_array = mz_array,
                                     sourceId = analytes_new,
                                     adduct = adducts_array)) %>% filter(!duplicated(sourceId))

  analytehaspathway_new = unique(analytehaspathway_new)
  source_new = unique(source_new)
  analytes_new = temp_mz_analyte$sourceId
  mzs_new = temp_mz_analyte$mz_array
  adducts_new = temp_mz_analyte$adduct


  ############  pathway analysis ##############
  verbose_message(message_text = "Begin metabolic pathway analysis ......" , verbose = verbose)
  analytes_rampids_df = merge(source_new %>% mutate(sourceId = tolower(sourceId)),
                              temp_mz_analyte%>% mutate(sourceId = tolower(sourceId)),
                              by = "sourceId")

  analytes_rampids = unique(analytes_rampids_df$rampId)

  # (1) Get candidate pathways
  # Get all analytes and number of analytes within a specific pathway

  source_non_duplicated = analytes_rampids_df[!duplicated(analytes_rampids_df$rampId),]

  # rampid = the subset of the database with our query data
  pathway_rampids = analytehaspathway_new[which(analytehaspathway_new$rampId %in% analytes_rampids),]
  pathway_rampids_count = pathway_rampids %>% dplyr::group_by(pathwayRampId) %>% dplyr::mutate(analytes_in_pathways  = n())

  # analytespathway_new. = the subset of the database with all pathways
  analytehaspathway_full = analytehaspathway_new %>%
    group_by(pathwayRampId) %>% dplyr::mutate(total_in_pathways = n())

  # Filter out too large/small pathways
  analytehaspathway_full =analytehaspathway_full[which(analytehaspathway_full$total_in_pathways>= min_path_size & analytehaspathway_full$total_in_pathways <= max_path_size),]

  # Generate a dataframe contains: the list of  IDs, the list of  names, the number of elements in pathway, the number of elements in our dataset, for each pathway
  if(pathway_all_info == T){
    unipathids = unique(pathway_rampids_count$pathwayRampId)
    sub_src = source_non_duplicated[which(source_non_duplicated$rampId  %in% pathway_rampids_count$rampId),]
    src_rid = sub_src$rampId
    src_cn = sub_src$commonName
    src_sid = sub_src$sourceId
    src_adduct = sub_src$adduct
    src_mz = sub_src$mz_array

    enrichment_df = pbapply::pblapply(1:length(unipathids), function(x){

      pathway_id = unipathids[x]
      pathway_info = pathway[which(pathway$pathwayRampId == pathway_id),]
      # get rampids associated with the pathway
      full_list = analytehaspathway_full[which(analytehaspathway_full$pathwayRampId == pathway_id)[1],4]
      screened_List_full = pathway_rampids_count[which(pathway_rampids_count$pathwayRampId == pathway_id),]

      # Get screened index
      met_ind = which(grepl(screened_List_full$rampId,
                            pattern = "RAMP_C_"))
      source_index_met = which(src_rid %in% screened_List_full$rampId[met_ind])
      source_index_gene = which(src_rid %in% screened_List_full$rampId[which(grepl(screened_List_full$rampId,
                                                                                   pattern = "RAMP_G_"))])
      #met
      ananlytes_name_list_met = paste0(src_cn[source_index_met], collapse = ";")
      ananlytes_id_list_met = paste0(src_sid[source_index_met], collapse = ";")

      #met_adductt
      ananlytes_mz_adduct = paste0(paste0(src_mz[source_index_met],"[",
                                          src_adduct[source_index_met]
                                          ,"]"), collapse = ";")
      #gene
      ananlytes_name_list_gene = paste0(src_cn[source_index_gene], collapse = ";")
      ananlytes_id_list_gene = paste0(src_sid[source_index_gene], collapse = ";")

      analytes_in_pathways = screened_List_full[1,4]
      total_in_pathways = full_list
      return_df = data.frame(pathway_name = pathway_info$pathwayName,
                             pathway_id = pathway_info$sourceId,
                             type = pathway_info$type,
                             pathwayCategory = pathway_info$pathwayCategory ,
                             metabolite_name_list=ananlytes_name_list_met,
                             metabolite_id_list= ananlytes_id_list_met,
                             total_in_pathways = total_in_pathways,
                             gene_name_list = ananlytes_name_list_gene,
                             gene_id_list = ananlytes_id_list_gene,
                             analytes_in_pathways = analytes_in_pathways,
                             adduct_info = ananlytes_mz_adduct)
      return(return_df)
    })
    enrichment_df = do.call(rbind, enrichment_df)

  }else{
    unipathids = unique(pathway_rampids_count$pathwayRampId)
    verbose_message(message_text = "Merging datasets" , verbose = verbose)
    analytehaspathway_sub = analytehaspathway_full[which(analytehaspathway_full$pathwayRampId %in% unipathids),] %>% filter(!duplicated(pathwayRampId))

    enrichment_df = base::merge(pathway_rampids_count[which(!duplicated(pathway_rampids_count$pathwayRampId)),], analytehaspathway_sub,
                                by = "pathwayRampId")

    enrichment_df = base::merge(enrichment_df, pathway, by = "pathwayRampId")

  }



  verbose_message(message_text = "Running test" , verbose = verbose)

  # (2) Conduct pathway enrichment
  total_inlist_analytes = length(unique(analytes_rampids_df$rampId))
  total_in_background = length(unique(analytehaspathway_full$rampId))

  verbose_message(message_text = "Calculating p value......" , verbose = verbose)

  enrichment_df = na.omit(enrichment_df)
  enrichment_df = enrichment_df %>% rowwise() %>% mutate(p_val = stats::fisher.test(matrix(
    c(
      # Detected  in pathway, in analytelist
      as.numeric(analytes_in_pathways),
      # Detected  in pathway, not in analytelist
      max(0,as.numeric(total_in_pathways - analytes_in_pathways)),
      # Detected  not in pathway, in analyte list
      max(0,as.numeric(total_inlist_analytes - total_in_pathways)),
      # Detected  not in pathway
      # Pathway elements not detected
      as.numeric(total_in_background -
                   total_inlist_analytes - total_in_pathways + analytes_in_pathways)
    ),
    2,
    2
  ),
  alternative = alternative)$p.value)
  enrichment_df = cbind(enrichment_df,
                        fdr = p.adjust(enrichment_df$p_val, method = "fdr")) %>% mutate(background_analytes_number = total_in_background)

  enrichment_df = enrichment_df%>% mutate(ratio = analytes_in_pathways/total_in_pathways)

  verbose_message(message_text = "Done!" , verbose = verbose)


  if(pathway_all_info == F){
    return =enrichment_df %>% dplyr::select(-c(pathwayRampId,rampId.y, pathwaySource.y)) %>% dplyr::select(pathwayName,
                                                                                                           sourceId,
                                                                                                           type,
                                                                                                           pathwayCategory,
                                                                                                           p_val,
                                                                                                           fdr,ratio,
                                                                                                           analytes_in_pathways,
                                                                                                           total_in_pathways) %>% arrange(p_val)
    colnames(return)[1:4] = c("pathway_name",
                              "pathway_id",
                              "type",
                              "pathwayCategory")

  }else{
    return = data.frame(enrichment_df) %>% dplyr::select(pathway_name,
                                                         pathway_id,
                                                         type,
                                                         pathwayCategory,
                                                         p_val,
                                                         fdr,ratio,
                                                         analytes_in_pathways,
                                                         total_in_pathways,
                                                         metabolite_name_list,
                                                         metabolite_id_list,
                                                         adduct_info,
                                                         gene_name_list,
                                                         gene_id_list)%>% arrange(p_val)

  }

  if (!is.null(pval_cutoff)){
    return = return %>% dplyr::filter(p_val <= pval_cutoff)
  }


  return(return)
}


#' Regional Pathway Enrichment
#'
#' This the function used to compute the gene/metabolites set enrichment for multi-omics spatial data
#'
#' @param SpaMTP A SpaMTP Seurat object contains spatial metabolomics(SM)/transcriptomics(ST) data or both, if contains SM data, it should be annotated via SpaMTP::AnnotateSM function.
#' @param ident A name character to specific the cluster vector for regions in `SpaMTP@meta.data` slot.
#' @param DE.list A list consisting of differential expression data.frames for each input modality. Within each data.frame column names MUST include 'cluster', 'gene', ('avg_log2FC' or 'logFC') and ('p_val_adj' or 'FDR').
#' @param analyte_types Vector of character strings defining which analyte types to use. Options can be c("genes"), c("metabolites") or both (default = c("genes", "metabolites")).
#' @param SM_assay A Character string defining describing slot name for spatial metabolomics data in SpaMTP to extract intensity values from (default = "SPM").
#' @param ST_assay A Character string defining describing slot name for spatial transcriptomics data in SpaMTP to extract RNA count values from (default = "SPT").
#' @param SM_slot The slot name containing the SM assay matrix data (default = "counts").
#' @param ST_slot The slot name containing the ST assay matrix data (default = "counts").
#' @param min_path_size The min number of metabolites in a specific pathway (default = 5).
#' @param max_path_size The max number of metabolites in a specific pathway (default = 500).
#' @param pval_cutoff_mets A numerical value defining the adjusted p value cutoff for significant differentially expressed metabolites. If `NULL` cutoff = `0.05` (default = 0.05).
#' @param pval_cutoff_genes A numerical value defining the adjusted p value cutoff for significant differentially expressed genes. If `NULL` cutoff = `0.05` (default = 0.05).
#' @param verbose Boolean indicating whether to show informative messages. If FALSE these messages will be suppressed (default = TRUE).
#'
#' @return A SpaMTP object with set enrichment on given analyte types.
#' @export
#'
#' @importFrom rlang %||%
#'
#' @examples
#' # SpaMTP = FindRegionalPathways(SpaMTP, polarity = "positive")
FindRegionalPathways = function(SpaMTP,
                                ident,
                                DE.list,
                                analyte_types = c("genes", "metabolites"),
                                SM_assay = "SPM",
                                ST_assay = "SPT",
                                SM_slot = "counts",
                                ST_slot = "counts",
                                min_path_size = 5,
                                max_path_size = 500,
                                pval_cutoff_mets = 0.05,
                                pval_cutoff_genes = 0.05,
                                verbose = TRUE) {
  ## Checks for ident in SpaMTP Object
  if (!(ident %in% colnames(SpaMTP@meta.data))) {
    stop(
      "Ident: ",
      ident,
      " not found in SpaMTP object's @meta.data slot ... Make sure the ident column is in your @metadata and is a factor!"
    )
  }
  cluster_vector = as.factor(SpaMTP@meta.data[[ident]])
  assignment = cluster_vector
  cluster = levels(cluster_vector)
  ## Checks for data in SM and/or ST assay
  if ("genes" %in% analyte_types) {
    if (is.null(SpaMTP@assays[[ST_assay]]@layers[[ST_slot]])) {
      stop(
        paste0(
          "No data exists in object[[",
          ST_assay,
          "]][",
          ST_slot,
          "] .. If you are using transcriptomic data with 'genes' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'genes' from analyte_tpes"
        )
      )
    } else{
      gene_matrix = Matrix::t(SpaMTP[[ST_assay]]@layers[[ST_slot]])
      if (length(cluster_vector) != nrow(gene_matrix)) {
        stop(
          "Please make sure the input ident is a vector the same length as the number of spots/cells in the gene assay!"
        )
      }
    }
  }
  if ("metabolites" %in% analyte_types) {
    if (is.null(SpaMTP@assays[[SM_assay]]@layers[[SM_slot]])) {
      stop(
        paste0(
          "No data exists in object[[",
          SM_assay,
          "]][",
          SM_slot,
          "] .. If you are using metabolic data with 'metabolites' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'metabolites' from analyte_tpes"
        )
      )
    } else{
      mass_matrix = Matrix::t(SpaMTP[[SM_assay]]@layers[[SM_slot]])
      if (length(cluster_vector) != nrow(mass_matrix)) {
        stop(
          "Please make sure the input ident is a vector the same length as the number of spots/cells in the metabolite assay!"
        )
      }
    }
  }
  # (2) Annotation
  if (is.null(SpaMTP@tools$db_3)) {
    stop(
      "@tools$db_3 is empty! No intermediate annotation data saved in SpaMTP object. Please run AnnotateSM() with save.intermediate = TRUE",
      "or save the database by setting filename = '...' and manually assign the annotation dataframe to @tools$db_3 <- [ ..."
    )
  }
  db_3 <- SpaMTP@tools$db_3
  db_3 = db_3 %>%
    tidyr::separate_rows(Isomers_IDs, IsomerNames, sep = "; ")

  verbose_message(message_text = "Query necessary data and establish pathway database" , verbose = verbose)

  db_3 = db_3 %>% dplyr::mutate(inputid = Isomers_IDs) %>%  dplyr::mutate(chem_source_id = inputid)
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
  # Get pathway db
  verbose_message(message_text = "Constructing pathway database ..." , verbose = verbose)
  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")

  pathway_db = split(chempathway$rampId, chempathway$pathwayName)
  pathway_db = pathway_db[which(!duplicated(tolower(names(pathway_db))))]
  pathway_db = pathway_db[lapply(pathway_db, length) >= min_path_size  &
                            lapply(pathway_db, length) <= max_path_size]

  gc()
  gsea_all_cluster = data.frame()
  all_ranks = list()
  pb3 = txtProgressBar(
    min = 0,
    max = length(cluster),
    initial = 0,
    style = 3
  )
  for (i in cluster) {
    i <- as.character(i)
    ranks <- c()
    if ("metabolites" %in% analyte_types) {
      ## metabolites
      DE_met <- DE.list[["metabolites"]]
      sub_db3 = DE_met[which(as.character(DE_met$cluster) == i), ] %>% dplyr::filter(p_val_adj <= pval_cutoff_mets %||% 0.05) %>% dplyr::filter(!duplicated(ramp_id))
      met_ranks = scale(sub_db3$avg_log2FC, center = 0)
      names(met_ranks) = sub_db3$ramp_id
      ranks <- c(ranks, met_ranks)
    }
    if ("genes" %in% names(DE.list)) {
      ## genes
      DE_rna <- DE.list[["genes"]]
      sub_de_gene = DE_rna[which(as.character(DE_rna$cluster) == i), ] %>% dplyr::filter(p_val_adj <= pval_cutoff_genes %||% 0.05) %>% dplyr::filter(!duplicated(rampId))
      ranks_gene_vector = scale(sub_de_gene$avg_log2FC, center = 0)
      names(ranks_gene_vector) = sub_de_gene$rampId
      # Genes and metabolites
      ranks <- c(ranks, ranks_gene_vector)
    }

    ranks = ranks[which(!duplicated(names(ranks)))]
    all_ranks[[i]] = ranks[is.finite(ranks)]

    gsea_result <- c()
    if (length(all_ranks[[i]]) > 0) {
      suppressWarnings({
        gsea_result = fgsea::fgsea(
          pathways =  pathway_db,
          stats = all_ranks[[i]],
          minSize = min_path_size,
          maxSize = max_path_size
        )  %>%  dplyr::mutate(Cluster_id = i)
      })

    } else {
      gsea_result <- data.table::data.table(
        pathway = character(0),
        pval = numeric(0),
        padj = numeric(0),
        log2err = numeric(0),
        ES = numeric(0),
        NES = numeric(0),
        size = integer(0),
        leadingEdge = list(),
        Cluster_id = i
      )
    }
    gsea_result = na.omit(gsea_result) %>% filter(!duplicated(pathway))
    short_source = source_df[which((source_df$rampId %in% names(all_ranks[[i]])) &
                                     !duplicated(source_df$rampId)), ]

    addtional_entry = do.call(rbind, lapply(1:nrow(gsea_result), function(x) {
      temp = unique(unlist(gsea_result$leadingEdge[x]))
      if ("metabolites" %in% analyte_types) {
        temp_ref =   sub_db3[which(sub_db3$ramp_id %in% temp), ] %>% dplyr::mutate(adduct_info = paste0(observed_mz, "[", Adduct, "]")) %>% dplyr::filter(!duplicated(adduct_info))
      }
      if ("genes" %in% analyte_types) {
        temp_rna = short_source[which((short_source$rampId %in% temp) &
                                        (grepl(short_source$rampId, pattern = "RAMP_G"))), ]
      }
      return(
        data.frame(
          adduct_info = if("metabolites" %in% analyte_types){paste0(temp_ref$adduct_info, collapse = ";")}else{""},
          leadingEdge_metabolites = if("metabolites" %in% analyte_types){paste0(sub(";.*", "", temp_ref$IsomerNames), collapse = ";")}else{""},
          leadingEdge_metabolites_id = if("metabolites" %in% analyte_types){paste0(temp_ref$chem_source_id, collapse = ";")}else{""},
          leadingEdge_genes = if("genes" %in% analyte_types){paste0(temp_rna$commonName, collapse = ";")}else{""},
          met_regulation = if("metabolites" %in% analyte_types){paste0(ifelse(ranks[which((names(ranks) %in% temp) &
                                                       (grepl(names(ranks), pattern = "RAMP_C")))] >= 0, "↑", "↓"), collapse = ";")}else{""},
          rna_regulation = if("genes" %in% names(DE.list)){paste0(ifelse(ranks[which((names(ranks) %in% temp) &
                                                       (grepl(names(ranks), pattern = "RAMP_G")))] >= 0, "↑", "↓"), collapse = ";")}else{""}
        )
      )
    }))
    gsea_result = cbind(gsea_result , addtional_entry)
    gsea_all_cluster = rbind(gsea_all_cluster, gsea_result)
    setTxtProgressBar(pb3, as.numeric(which(cluster == i)))
  }
  close(pb3)

  gsea_all_cluster <- na.omit(gsea_all_cluster)%>%
    dplyr::mutate(group_importance = sum(abs(NES)))
  colnames(gsea_all_cluster)[1] = "pathwayName"
  gsea_all_cluster = merge(gsea_all_cluster, pathway, by = "pathwayName")
  return(gsea_all_cluster)
}




#' Runs multilevel Monte-Carlo variant for performing gene sets co-regulation analysis using the RAMP_DB metabolite/gene database.
#'
#' This function is adapted from the [fgsea::geseca](https://github.com/alserglab/fgsea/blob/master/R/geseca-multilevel.R) package to identify significantly expressed RAMP_DB pathways based on an expression/feature embedding matrix.
#'
#' @param E expression matrix, rows corresponds to RAMP_IDs, columns corresponds to cell barcodes.
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded (default = 1).
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded (default = `nrow(E) - 1`).
#' @param center a logical value indicating whether the gene expression should be centered to have zero mean before the analysis takes place (default = TRUE).
#' @param scale a logical value indicating whether the gene expression should be scaled to have unit variance before the analysis takes place (default = FALSE).
#' @param sampleSize sample size for conditional sampling (default = 101).
#' @param eps This parameter sets the boundary for calculating P-values (default = 1e-50).
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply (default = NULL).
#' @param nPermSimple Number of permutations in the simple geseca implementation for preliminary estimation of P-values (default = 1000).
#'
#' @return A table with GESECA results. Each row corresponds to a tested RAMP_DB pathway.
#' @export
#'
#' @examples
#' # E <- SpaMTP@reductions$pca.rev@feature.loadings
#' # sig_pathways <- RunRAMPgeseca(E, minSize=15, maxSize=500)
RunRAMPgeseca <- function(E,
                          minSize     = 1,
                          maxSize     = nrow(E) - 1,
                          center      = TRUE,
                          scale       = FALSE,
                          sampleSize  = 101,
                          eps         = 1e-50,
                          nproc       = 0,
                          BPPARAM     = NULL,
                          nPermSimple = 1000){

  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")

  pathway_db = split(chempathway$rampId, chempathway$pathwayName)
  pathway_db = pathway_db[which(!duplicated(tolower(names(pathway_db))))]
  pathway_db = pathway_db[lapply(pathway_db, length) >= minSize  &
                            lapply(pathway_db, length) <= maxSize]

  gesecaRes <- fgsea::geseca(pathway_db, E, minSize = minSize, maxSize = maxSize, center = center, scale = scale,sampleSize = sampleSize, eps = eps, nproc = nproc, BPPARAM = BPPARAM, nPermSimple = nPermSimple)

  return(gesecaRes)

}


#' Create a Pathway Assay from Gene or Metabolite Data
#'
#' This function creates a new assay within the provided SpaMTP Seurat object which contains features (either genes or metabolites) labeled by their respective RAMP ID. This assay can be used for running feature set co-regulation analysis (based on [GSCA](https://doi.org/10.1093/bioinformatics/btp502)).
#'
#' @param SpaMTP A SpaMTP Seurat object containing either spatial metabolic or transcriptomic data
#' @param analyte_type Character string specifying the type of analytes to process.Must be either "genes" or "metabolites" (default = "metabolites").
#' @param assay Character string specifying the name of the assay to use as source data (default = "Spatial").
#' @param slot Character string specifying which slot in the assay to use as source data (default = "counts").
#' @param new_assay Character string specifying the name of the new assay to create (default = "pathway").
#' @param verbose Boolean logical value indicating whether to print verbose messages during execution. (default = TRUE).
#'
#' @return A SpaMTP object with a new assay added, containing respective gene/metabolite data formatted based on RAMP_db IDs.
#' @export
#'
#' @importFrom dplyr mutate group_by summarise ungroup
#' @importFrom tidyr separate_rows
#' @importFrom data.table as.data.table
#' @importFrom SeuratObject CreateAssay5Object
#'
#' @examples
#' ## Create a pathway assay from metabolite data
#' #spamtp_obj <- CreatePathwayAssay(spamtp_obj, analyte_type = "metabolites", assay = "SPM", new_assay = "pathway")
#'
#' ## Create a pathway assay from gene data with verbose output
#' #spamtp_obj <- CreatePathwayAssay(spamtp_obj, analyte_type = "genes", assay = "SPT", new_assay = "gene_pathway", verbose = TRUE)
CreatePathwayAssay <- function(SpaMTP, analyte_type = "metabolites", assay = "Spatial", slot = "counts", new_assay = "pathway", verbose = TRUE){

  if(!analyte_type %in% c("genes", "metabolites")){
    stop("Incorrect `analyte_type` provided! must be either 'genes' or 'metabolites'. Please provided the correct analyte matching the selected assay data.")
  }

  if (analyte_type == "genes") {
    if (is.null(SpaMTP@assays[[assay]][slot])) {
      stop(
        paste0(
          "No data exists in object[[",
          assay,
          "]][",
          slot,
          "] .. If you are using transcriptomic data with 'genes' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'genes' from analyte_tpes"
        )
      )
    } else{
      matrix <- as.data.frame(SpaMTP[[assay]][slot])
      matrix$commonName <- toupper(rownames(matrix))
      matrix = merge(matrix,
                     unique(source_df[which(grepl(source_df$rampId, pattern = "RAMP_G")), ][c("rampId" ,"commonName")]),
                     by = "commonName")
      dupe_list <- split(which(matrix$rampId %in% matrix$rampId[duplicated(matrix$rampId)]),
                         matrix$rampId[matrix$rampId %in% matrix$rampId[duplicated(matrix$rampId)]])

      meta.data <- matrix[c("rampId" ,"commonName")] %>%
        group_by(rampId) %>%
        summarise(commonName = paste(commonName, collapse = "; ")) %>%
        ungroup()


      matrix$commonName <- NULL

    }
  }
  if (analyte_type == "metabolites") {
    if (is.null(SpaMTP@assays[[assay]][slot])) {
      stop(
        paste0(
          "No data exists in object[[",
          assay,
          "]][",
          slot,
          "] .. If you are using metabolic data with 'metabolites' in 'analyte_types', please ensure this dataslot exists within your SpaMTP object, else remove 'metabolites' from analyte_tpes"
        )
      )
    } else{
      matrix <- as.data.frame(SpaMTP[[assay]][slot])

      # (2) Annotation
      if (is.null(SpaMTP@tools$db_3)) {
        stop(
          "@tools$db_3 is empty! No intermediate annotation data saved in SpaMTP object. Please run AnnotateSM() with save.intermediate = TRUE",
          "or save the database by setting filename = '...' and manually assign the annotation dataframe to @tools$db_3 <- [ ..."
        )
      }
      db_3 <- SpaMTP@tools$db_3
      db_3 = db_3 %>%
        tidyr::separate_rows(Isomers_IDs, IsomerNames, sep = "; ")

      verbose_message(message_text = "Query necessary data and establish pathway database" , verbose = verbose)

      db_3 = db_3 %>% dplyr::mutate(inputid = Isomers_IDs) %>%  dplyr::mutate(chem_source_id = inputid)
      verbose_message(message_text = "Query db for addtional matching" , verbose = verbose)
      db_3 = merge(chem_props, db_3, by = "chem_source_id")
      ### Adding DE Results
      db_3 = db_3 %>% mutate(mz_name = paste0("mz-", db_3$observed_mz))
      db_3 <- db_3[c("mz_name",  "ramp_id")]
      db_3 <- db_3 %>% distinct()
      matrix$mz_name <- rownames(SpaMTP@assays[[assay]])
      matrix = merge(db_3 , matrix, by = "mz_name")

      meta.data <- matrix[c("ramp_id" ,"mz_name")] %>%
        group_by(ramp_id) %>%
        summarise(mz_name = paste(mz_name, collapse = "; ")) %>%
        ungroup()

      meta.data$rampId <- meta.data$ramp_id
      meta.data <- meta.data[c("rampId","mz_name")]

      rm(db_3)

      dupe_list <- split(which(matrix$ramp_id %in% matrix$ramp_id[duplicated(matrix$ramp_id)]),
                         matrix$ramp_id[matrix$ramp_id %in% matrix$ramp_id[duplicated(matrix$ramp_id)]])

      matrix$mz_name <- NULL
      matrix$rampId <- matrix$ramp_id
      matrix$ramp_id <- NULL

    }
  }


  matrix <- data.table::as.data.table(matrix)

  if(length(dupe_list) > 0){

    verbose_message(message_text = paste0("Some RAMP_IDs have multiple mapped analytes. There are: ",
                                         length(names(dupe_list))) , verbose = verbose)

    # Select rows using 'dupe_list' indices
    merged_data <- matrix[unlist(dupe_list),]

    # Compute mean for numeric columns by 'ramp_id'
    merged_data <- merged_data[, lapply(.SD, mean, na.rm = TRUE), by = rampId]


    matrix <- matrix[!rownames(matrix) %in% unlist(unname(dupe_list)),]

    matrix <- rbind(matrix,  merged_data)

  }

  rownames(matrix) <- matrix$rampId
  matrix$rampId <- NULL


  SpaMTP[[new_assay]] <- SeuratObject::CreateAssay5Object(counts = matrix)

  message("Warning: Restoring feature names to contain '_' ...")

  # Manually restore underscores
  rownames(SpaMTP[[new_assay]]) <- gsub("-", "_", x = rownames(SpaMTP[[new_assay]]))

  SpaMTP[[new_assay]]@meta.data$rampId <- rownames(SpaMTP[[new_assay]])
  SpaMTP[[new_assay]]@meta.data <- merge(SpaMTP[[new_assay]]@meta.data, meta.data, by = "rampId", all = TRUE)


  return(SpaMTP)

}



############################# Annotation Estimation using Pathway Results ########################################



#' Create a SpaMTP Seurat Object containing expression values for all present pathways
#'
#' This function computes pathway-level scores from analyte-level expression data
#' and stores the results as a new assay in a Seurat object. Each pathway score
#' is calculated as the scaled mean expression of the analytes associated with
#' that pathway, adjusted by the square root of the pathway size.
#'
#' @param object A SpaMTP Seurat object containing the expression data.
#' @param assay Character. Name of the assay to extract analyte expression from. If no value is assigned, the DefaultAssay of the SpaMTP Seurat Object will be used (default = `DefaultAssay(object)`).
#' @param slot Character. Which data slot to use (e.g., "scale.data") (default = "scale.data").
#' @param new.assay Character. Name of the new assay where pathway scores will be stored (defaults = "pathway").
#' @param remove.nans Logical. Whether to remove pathways with all NaN values (e.g., no matched analytes) (defaults = TRUE).
#'
#' @return A SpaMTP Seurat object with a new assay containing pathway-level expression scores.
#'         Feature names are adjusted to use underscores instead of dashes.
#'
#' @export
#'
#' @examples
#' #object <- CreatePathwayObject(seurat_obj, assay = "RNA", slot = "scale.data")
CreatePathwayObject <- function(object,
                                assay=SeuratObject::DefaultAssay(object),
                                slot = "scale.data",
                                new.assay = "pathway",
                                remove.nans = TRUE
) {

  chempathway = merge(analytehaspathway, pathway, by = "pathwayRampId")
  pathway_db <- split(chempathway$rampId, chempathway$pathwayRampId)
  pathway_db <- pathway_db[!duplicated(tolower(names(pathway_db)))]

  x <- Seurat::GetAssay(object, assay)
  E <- x[slot]

  pathway_sums <- list()
  for (i in seq_along(pathway_db)) {
    pathway <- pathway_db[[i]]
    pathway <- intersect(unique(pathway), rownames(E))
    score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
    score <- scale(score, center=TRUE, scale=TRUE)
    pathway_sums[[names(pathway_db)[i]]] <- score
  }

  pathway_mtx <- do.call(cbind, pathway_sums)
  colnames(pathway_mtx) <- names(pathway_db)

  pathway_mtx <- t(pathway_mtx)

  if (remove.nans){
    nan_rows <- apply(pathway_mtx, 1, function(row) all(is.nan(row)))
    pathway_mtx <- pathway_mtx[!nan_rows,]
  }

  object[[new.assay]] <- SeuratObject::CreateAssay5Object(counts = pathway_mtx)

  message("Warning: Restoring feature names to contain '_' ...")

  rownames(object[[new.assay]]) <- gsub(pattern = "-", replacement = "_", x = rownames(object[[new.assay]]))

  object[[new.assay]]@meta.data <- chempathway %>%
    filter(pathwayRampId %in% rownames(object[[new.assay]])) %>%
    select(pathwayRampId, pathwayName) %>%
    distinct()

  return(object)
}









############################# PATHWAY HELPER FUNCTIONS ########################################

#' Helper function for building a pathway db based on detected
#'
#' @param input_id Vector of characters defining the detected .
#' @param analytehaspathway A dataframe containing RAMP_pathway ID's.
#' @param chem_props A database containing the chemical properties and metadata of each RAMP_DB analyte.
#' @param pathway A dataframe containing RAMP_DB pathways and their relative metadata
#'
#' @return A analyte database containing corresponding pathways associated with each detected metabolite
#'
#' @examples
#' #HELPER FUNCTION
get_analytes_db <- function(input_id,analytehaspathway,chem_props,pathway) {

  rampid = unique(chem_props$ramp_id[which(chem_props$chem_source_id %in% unique(input_id))])
  #
  pathway_ids = unique(analytehaspathway$pathwayRampId[which(analytehaspathway$rampId %in% rampid)])

  analytes_db = lapply(pathway_ids, function(x) {
    content = analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == x)]
    content = content[which(grepl(content, pattern = "RAMP_C"))]
    return(content)
  })
  analytes_db_name = unlist(lapply(pathway_ids, function(x) {
    name = pathway$pathwayName[which(pathway$pathwayRampId == x)]
    return(name)
  }))
  names(analytes_db) = analytes_db_name
  return(analytes_db)
}


