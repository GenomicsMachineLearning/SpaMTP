#' Calculates Significant Metabolic Pathways using a Fisher Exact Test
#'
#' @param Analyte A list of analytes with 3 elements, namely "mz", "genes" and "", each is comprised of the corresponding labels, for ,
#' Supported  format including, X stands for upper case of the cooresponding ID in each database: "hmdb:HMDBX", "chebi:X", "pubchem:X","wikidata:X" ,"kegg:X" ,"CAS:X","lipidbank:X","chemspider:X","	LIPIDMAPS:X"
#' Supported gene data format including: "entrez:X", "gene_symbol:X", "uniprot:X", "ensembl:X", "hmdb:HMDBPX"
#' Supported mz format: any string or numeric vector contains the m/z
#' @param analyte_type = "" or "gene" or "mz", or a vector contains any combinations of them (default = c("mz", "genes")).
#' @param polarity Character string defining the polarity of the MALDI experiment. Inputs must be either 'positive', 'negative' or 'neutral' (default = NULL).
#' @param ppm_error Integer defining the ppm threshold that matched analytes must be between (default = 10).
#' @param max_path_size The max number of  in a specific pathway (default = 500).
#' @param min_path_size The min number of  in a specific pathway (default = 5).
#' @param alternative The hypothesis of the fisher exact test (default = "greater").
#' @param pathway_all_info Whether to included all genes/ screened in the return (default = FALSE).
#' @param pval_cutoff The cut off of raw p value to retain the pathways (default = 0.05).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return a dataframe with the relevant pathway information
#' @export
#'
#' @import dplyr
#' @import stringr
#'
#' @examples
#' # FishersPathwayAnalysis(Analyte = mzs, analyte_type = "mz", ppm_error = 3)
FishersPathwayAnalysis <- function (Analyte,
                                    analyte_type = c("mz"),
                                    polarity = NULL,
                                    ppm_error = 10,
                                    max_path_size = 500,
                                    min_path_size = 5,
                                    alternative = "greater",
                                    pathway_all_info = FALSE,
                                    pval_cutoff = 0.05,
                                    verbose = TRUE)
{

  if((!"mz" %in% analyte_type) & (!"" %in% analyte_type) & (!"genes" %in% analyte_type)){

    stop(
      "analyte_type was not specified correctly.  Please specify one of the following options: , genes"
    )
  }
  now <- proc.time()

  verbose_message(message_text = "Fisher Testing ......", verbose = verbose)

  pathwayRampId <- rampId <- c()
  # Get the RaMP ids for /genes
  convert_to_rows <- function(row,
                              pattern) {
    identifiers = data.frame()
    for (i in which(grepl(row,
                          pattern = pattern))) {
      temp = unlist(strsplit(unlist(row[i]), pattern))
      identifiers <- rbind(identifiers,
                           temp)
    }
    identifiers = t(identifiers)
    colnames(identifiers) = colnames(row)[which(grepl(row,
                                                      pattern = pattern) ==
                                                  T)]
    return(cbind(row[which(grepl(row,
                                 pattern = pattern) == F)][rep(1, times = nrow(identifiers)), ],
                 identifiers))
  }

  if ( "" %in% analyte_type) {
    analytes_met = Analyte[[""]]
    source_met = source_df[which(grepl(source_df$rampId, pattern = "RAMP_C") == T),]
    analytehaspathway_met = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_C") == T),]
    analyte_met = analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),]
  }
  if ("genes" %in% analyte_type) {
    if (!("genes" %in% names(Analyte))){
      stop("Cannot find gene list not present in Analyte object .... Please provide SpaMTP assay containing gene names, or remove 'gene' from analyte_type input!")
    }
    analytes_rna = Analyte[["genes"]]
    source_rna = source_df[which(grepl(source_df$rampId, pattern = "RAMP_G") == T),]
    analytehaspathway_rna = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_G") == T),]
    analyte_rna = analyte[which(grepl(analyte$rampId, pattern = "RAMP_G") == T),]
  }

  if ("mz" %in% analyte_type) {
    analytes_mz = Analyte[["mz"]]

    # Since this file was tested in positive ion mode
    db = rbind(Chebi_db,
               Lipidmaps_db,
               HMDB_db)

    gc()
    if (polarity == "positive") {
      test_add_pos <- adduct_file$adduct_name[which(adduct_file$charge > 0)]

      # 1) Filter DB by adduct.
      db_1 <- db_adduct_filter(db, test_add_pos, polarity = "pos", verbose = verbose)
    } else if (polarity == "negative") {
      test_add_neg <- adduct_file$adduct_name[which(adduct_file$charge < 0)]

      # 1) Filter DB by adduct.
      db_1 <- db_adduct_filter(db, test_add_neg, polarity = "neg", verbose = verbose)
    } else if (polarity == "neutral") {

      # 1) Filter DB by adduct.
      db_1 <- db %>% mutate("M" = `M-H ` + 1.007276)
    }else{
      stop("Please enter correct polarity from: 'positive', 'negative', 'neutral'")
    }

    # 2) only select natural elements
    db_2 <- formula_filter(db_1)

    # 3) search db against mz df return results
    verbose_message(message_text = "search db against mz df return results", verbose = verbose)

    input_mz = data.frame(cbind(
      row_id = 1:length(analytes_mz),
      mz = as.numeric(str_extract(analytes_mz, pattern = "\\d+\\.?\\d*"))
    ))
    ppm_error <- ppm_error
    db_3 <- proc_db(data.frame(input_mz), db_2, ppm_error)
    # Expand the isomer entries

    verbose_message(message_text = "Expanding database to extract all potential ", verbose = verbose)

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
  if("mz" %in% analyte_type){
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

  if("" %in% analyte_type){
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

  if("genes" %in% analyte_type){
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
  # Pathway enrichment

  ############  pathway analysis ##############
  verbose_message(message_text = "Begin metabolic pathway analysis ......" , verbose = verbose)
  analytes_rampids_df = merge(source_new %>% mutate(sourceId = tolower(sourceId)),
                              temp_mz_analyte%>% mutate(sourceId = tolower(sourceId)),
                              by = "sourceId")

  analytes_rampids = unique(analytes_rampids_df$rampId)
  # for(k in 1:length(unique_analytes)){
  #   pattern = unique_analytes[k]
  #   analytes_rampids = c(analytes_rampids,
  #                        unique(source$rampId[which(grepl(source$sourceId, pattern = pattern,
  #                                                         ignore.case = T))]))
  # }
  # analytes_rampids = unique(na.omit(analytes_rampids))
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
      #ananlytes_id_df = analytehaspathway_new[which(analytehaspathway_new$pathwayRampId == unique(pathway_rampids$pathwayRampId)[x]),]
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
    #enrichment_df = enrichment_df %>% dplyr::filter(!duplicated(pathway_name))
  }else{
    unipathids = unique(pathway_rampids_count$pathwayRampId)
    verbose_message(message_text = "Merging datasets" , verbose = verbose)
    analytehaspathway_sub = analytehaspathway_full[which(analytehaspathway_full$pathwayRampId %in% unipathids),] %>% filter(!duplicated(pathwayRampId))

    enrichment_df = base::merge(pathway_rampids_count[which(!duplicated(pathway_rampids_count$pathwayRampId)),], analytehaspathway_sub,
                                by = "pathwayRampId")
    #colnames(enrichment_df)[which(colnames(enrichment_df) == "count")] = "total_in_pathways"
    enrichment_df = base::merge(enrichment_df, pathway, by = "pathwayRampId")
    #enrichment_df = enrichment_df %>% dplyr::filter(!duplicated(pathwayName))
  }

  # abbb = analytehaspathway_new[which(analytehaspathway_new$rampId %in% analytes_rampids)[1:100],] %>% rowwise() %>%
  # dplyr::mutate(ananlytes_id_list = list(analytehaspathway_new$rampId[which(analytehaspathway_new$pathwayRampId == pathwayRampId)])) %>%
  # dplyr::count(pathwayRampId,
  #              ananlytes_id_list,
  #              sort = T,
  #              name = "analytes_in_pathways")

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

  verbose_message(message_text = "P value obtained" , verbose = verbose)

  # (5) Append pathway information to the original df
  gc()
  # (6) Append  information to the original df
  # Paste back the original Ids
  # (7) Reduce the dataframe with respected to the User input pathway size

  verbose_message(message_text = "Done" , verbose = verbose)

  gc()
  #return(enrichment_df_with_both_info %>% select(-c(
  #  pathwayRampId,
  #  ananlytes_id_list,
  #  screened_analytes
  #)))
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
  return(return)
}


#' This the function used to compute the gene/metabolites set enrichment for multi-omics spatial data
#'
#' @param SpaMTP A SpaMTP Seurat object contains spatial metabolomics(SM)/transcriptomics(ST) data or both, if contains SM data, it should be annotated via SpaMTP::AnnotateSM function.
#' @param ident A name character to specific the cluster vector for regions in `SpaMTP@meta.data` slot.
#' @param DE.list A list consist of differetial expression output from FindAllMarkers() function, with items in same order as analyte_types.
#' @param analyte_types Vector of character strings defining which analyte types to use. Options can be c("genes"), c("metabolites") or both (default = c("genes", "metabolites")).
#' @param adduct Vector of character strings defining adducts to use for analysis (e.g. c("M+K","M+H ")). For all possible adducts please visit [here](https://github.com/GenomicsMachineLearning/SpaMTP/blob/main/R/MZAnnotation.R#L305). If NULL will take the full list of SpaMTP::adduct_file$adduct_name (default = NULL).
#' @param SM_assay A Character string defining descrbing slot name for spatial metabolomics data in SpaMTP to extract intensity values from (default = "SPM").
#' @param ST_assay A Character string defining descrbing slot name for spatial transcriptomics data in SpaMTP to extract RNA count values from (default = "SPT").
#' @param SM_slot The slot name containing the SM assay matrix data (default = "counts").
#' @param ST_slot The slot name containing the ST assay matrix data (default = "counts").
#' @param max_path_size The max number of metabolites in a specific pathway (default = 500).
#' @param min_path_size The min number of metabolites in a specific pathway (default = 5).
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion `[ion mass,m/z]`/`[Full width at half height]` (default = 30000).
#' @param pval_cutoff_pathway A numerical value between 0 and 1 describe the cutoff adjusted p value for the permutation test used to compute output pathways
#' @param pval_cutoff_mets A numerical value between 0 and 1 describe the cutoff adjusted p value for the differential expression analysis for metabolites
#' @param pval_cutoff_genes A numerical value between 0 and 1 describe the cutoff adjusted p value for the differential expression analysis for RNAs
#' @param verbose A boolean value indicates whether verbose is shown
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
                                tof_resolution = 30000,
                                min_path_size = 5,
                                max_path_size = 500,
                                pval_cutoff_pathway = NULL,
                                pval_cutoff_mets = NULL,
                                pval_cutoff_genes = NULL,
                                verbose = T) {
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
      temp_ref =   sub_db3[which(sub_db3$ramp_id %in% temp), ] %>% dplyr::mutate(adduct_info = paste0(observed_mz, "[", Adduct, "]")) %>% dplyr::filter(!duplicated(adduct_info))
      temp_rna = short_source[which((short_source$rampId %in% temp) &
                                      (grepl(short_source$rampId, pattern = "RAMP_G"))), ]
      return(
        data.frame(
          adduct_info = paste0(temp_ref$adduct_info, collapse = ";"),
          leadingEdge_metabolites = paste0(sub(";.*", "", temp_ref$IsomerNames), collapse = ";"),
          leadingEdge_metabolites_id = paste0(temp_ref$chem_source_id, collapse = ";"),
          leadingEdge_genes = paste0(temp_rna$commonName, collapse = ";"),
          met_regulation = paste0(ifelse(ranks[which((names(ranks) %in% temp) &
                                                       (grepl(names(ranks), pattern = "RAMP_C")))] >= 0, "↑", "↓"), collapse = ";"),
          rna_regulation = paste0(ifelse(ranks[which((names(ranks) %in% temp) &
                                                       (grepl(names(ranks), pattern = "RAMP_G")))] >= 0, "↑", "↓"), collapse = ";")
        )
      )
    }))
    gsea_result = cbind(gsea_result , addtional_entry)
    gsea_all_cluster = rbind(gsea_all_cluster, gsea_result)
    setTxtProgressBar(pb3, as.numeric(which(cluster == i)))
  }
  close(pb3)
  gsea_all_cluster <- na.omit(gsea_all_cluster)%>% #dplyr::mutate(
    dplyr::mutate(group_importance = sum(abs(NES)))
  colnames(gsea_all_cluster)[1] = "pathwayName"
  gsea_all_cluster = merge(gsea_all_cluster, pathway, by = "pathwayName")
  return(gsea_all_cluster)
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


