
#' Calculates Significant Metabolic Pathways using a Fisher Exact Test
#'
#' @param Analyte A list of analytes with 3 elements, namely "mz", "genes" and "metabolites", each is comprised of the corresponding labels, for metabolites,
#' Supported metabolites format including, X stands for upper case of the cooresponding ID in each database: "hmdb:HMDBX", "chebi:X", "pubchem:X","wikidata:X" ,"kegg:X" ,"CAS:X","lipidbank:X","chemspider:X","	LIPIDMAPS:X"
#' Supported gene data format including: "entrez:X", "gene_symbol:X", "uniprot:X", "ensembl:X", "hmdb:HMDBPX"
#' Supported mz format: any string or numeric vector contains the m/z
#' @param analyte_type = "metabolites" or "gene" or "mz", or a vector contains any combinations of them (default = c("mz", "genes")).
#' @param polarity Character string defining the polarity of the MALDI experiment. Inputs must be either 'positive', 'negative' or 'neutral' (default = NULL).
#' @param ppm_error Integer defining the ppm threshold that matched analytes must be between (default = 10).
#' @param max_path_size The max number of metabolites in a specific pathway (default = 500).
#' @param min_path_size The min number of metabolites in a specific pathway (default = 5).
#' @param alternative The hypothesis of the fisher exact test (default = "greater").
#' @param pathway_all_info Whether to included all genes/metabolites screened in the return (default = FALSE).
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
  data(pathway, package = "SpaMTP")
  if((!"mz" %in% analyte_type) & (!"metabolites" %in% analyte_type) & (!"genes" %in% analyte_type)){
    stop(
      "analyte_type was not specified correctly.  Please specify one of the following options: metabolites, genes"
    )
  }
  now <- proc.time()

  verbose_message(message_text = "Fisher Testing ......", verbose = verbose)

  pathwayRampId <- rampId <- c()
  # Get the RaMP ids for metabolites/genes
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

  if ( "metabolites" %in% analyte_type) {
    analytes_met = Analyte[["metabolites"]]
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
      # Using Chris' pipeline for annotation
      # 1) Filter DB by adduct.
      db_1 <- db_adduct_filter(db, test_add_pos, polarity = "pos", verbose = verbose)
    } else if (polarity == "negative") {
      test_add_neg <- adduct_file$adduct_name[which(adduct_file$charge < 0)]
      # Using Chris' pipeline for annotation
      # 1) Filter DB by adduct.
      db_1 <- db_adduct_filter(db, test_add_neg, polarity = "neg", verbose = verbose)
    } else if (polarity == "neutral") {
      # Using Chris' pipeline for annotation
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

    verbose_message(message_text = "Expanding database to extract all potential metabolites", verbose = verbose)

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

  if("metabolites" %in% analyte_type){
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

  ############ Metabolites pathway analysis ##############
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

  # Generate a dataframe contains: the list of metabolites IDs, the list of metabolites names, the number of elements in pathway, the number of elements in our dataset, for each pathway
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
      # Detected metabolites in pathway, in analytelist
      as.numeric(analytes_in_pathways),
      # Detected metabolites in pathway, not in analytelist
      max(0,as.numeric(total_in_pathways - analytes_in_pathways)),
      # Detected metabolites not in pathway, in analyte list
      max(0,as.numeric(total_inlist_analytes - total_in_pathways)),
      # Detected metabolites not in pathway
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
  # (6) Append metabolites information to the original df
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



##########################################################################################


#' PCA driven pathway analysis of a SpaMTP object
#'
#' @param SpaMTP SpaMTP Seurat class object that contains spatial metabolic information.
#' @param path Character string defining the output path for the visualization (default = getwd()).
#' @param ppm_error is the parts-per-million error tolerance of matching m/z value with potential metabolites
#' @param ion_mode is only needed when ppm_error is not specified, used to access the ppm_error based on polarity. Inputs must be either 'positive' or 'negative'(default = NULL).
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion `[ion mass,m/z]`/`[Full width at half height]` (default = 30000).
#' @param num_retained_component is an integer value to indicated preferred number of PCs to retain (default = NULL).
#' @param variance_explained_threshold Numeric value defining the explained variance threshold (default = 0.9).
#' @param resampling_factor is a numerical value >0, indicate how you want to resample the size of original matrix (default = 1).
#' @param p_val_threshold is the p value threshold for pathways to be significant (default = 0.05).
#' @param byrow is a boolean to indicates whether each column of the matrix is built byrow or bycol (default = FALSE).
#' @param assay Character string defining the SpaMTP assay to extract intensity values from (default = "SPM").
#' @param slot Character string defining the assay slot contatin ght intesity values (default = "counts").
#' @param flip_plot Boolean defining whether to rotate the plot 90 degrees (default = FALSE).
#' @param show_variance_plot Boolean indicating weather to display the variance plot output by this analysis (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#'
#' @return Dataframe containing the pathway enrichment results for each PC
#' @export
#'
#' @import dplyr
#'
#' @examples
#' # PCAPathwayAnalysis(SpaMTP.obj, ion_mode = "positive",resampling_factor = 2,p_val_threshold = 0.05, verbose = FALSE)
PCAPathwayAnalysis <- function(SpaMTP,
                                                path = getwd(),
                                                ppm_error = NULL,
                                                ion_mode = NULL,
                                                tof_resolution = 30000,
                                                num_retained_component = NULL,
                                                variance_explained_threshold = 0.9,
                                                resampling_factor = 1,
                                                p_val_threshold = 0.05,
                                                byrow = FALSE,
                                                assay = "SPM",
                                                slot = "counts",
                                                flip_plot = FALSE,
                                                show_variance_plot= FALSE,
                                                verbose = TRUE) {


  # PCA analysis
  verbose_message(message_text = "Scaling original matrix", verbose = verbose)

  mass_matrix = Matrix::t(SpaMTP[[assay]]@layers[[slot]])

  coords <- GetTissueCoordinates(SpaMTP)[c("x", "y")]

  if (flip_plot){
    colnames(coords) <- c("y", "x")
  }


  mass_matrix_with_coord = cbind(coords,
                                 as.matrix(mass_matrix))
  width = max(coords[c("y")]) #- min(coords[c("y")])
  height = max(coords[c("x")]) # - min(coords[c("x")])

  if (!is.null(resampling_factor) & resampling_factor!= 1) {
    verbose_message(message_text = "Running matrix resampling...." , verbose = verbose)
    pb = txtProgressBar(
      min = 0,
      max = ncol(mass_matrix),
      initial = 0,
      style = 3
    )
    if (!is.numeric(resampling_factor)) {
      stop("Please enter correct resampling_factor")
    }
    new.width = as.integer(width / resampling_factor)
    new.height = as.integer(height / resampling_factor)

    resampled_mat = matrix(nrow =  new.height * new.width)
    for (i in 1:ncol(mass_matrix)) {
      temp_mz_matrix = matrix(mass_matrix[, i],
                              ncol = width,
                              nrow = height,
                              byrow = byrow)
      resampled_temp = ResizeMat(temp_mz_matrix, c(new.width,
                                                   new.height))
      resampled_mat = cbind(resampled_mat, as.vector(resampled_temp))
      setTxtProgressBar(pb, i)
    }
    close(pb)
    resampled_mat = resampled_mat[,-1]
    colnames(resampled_mat) = colnames(mass_matrix)

    verbose_message(message_text = "Resampling finished!" , verbose = verbose)

    gc()
  }else{
    resampled_mat = mass_matrix
    new.width = as.integer(width)
    new.height = as.integer(height)
  }

  verbose_message(message_text = "Running the principal component analysis (can take some time)" , verbose = verbose)

  # Runing PCA

  resampled_mat_standardised = as.matrix(Matrix::t(
    Matrix::t(resampled_mat) - Matrix::colSums(resampled_mat) / nrow(resampled_mat)
  ))

  verbose_message(message_text = "Computing the covariance" , verbose = verbose)
  cov_mat <- cov(resampled_mat_standardised)

  verbose_message(message_text = "Computing the eigenvalue/eigenvectors", verbose = verbose)
  eigen_result <- eigen(cov_mat)
  gc()
  # Extract eigenvectors and eigenvalues
  eigenvectors <- eigen_result$vectors
  eigenvalues <- eigen_result$values

  verbose_message(message_text = "Computing PCA", verbose = verbose)

  pc = pbapply::pblapply(1:ncol(resampled_mat_standardised), function(i) {
    temp = resampled_mat_standardised[, 1] * eigenvectors[1, i]
    for (j in 2:ncol(resampled_mat_standardised)) {
      temp = temp + resampled_mat_standardised[, j] * eigenvectors[j, i]
    }
    return(temp)
  })
  pc = do.call(cbind, pc)
  colnames(pc) = paste0("PC", 1:ncol(eigenvectors))
  # make pca object
  colnames(eigenvectors) = paste0("PC", 1:ncol(eigenvectors))
  rownames(eigenvectors) = colnames(resampled_mat)
  pca = list(
    sdev = sqrt(eigenvalues),
    rotation = eigenvectors,
    center = Matrix::colSums(resampled_mat) / nrow(resampled_mat),
    scale = FALSE,
    x = pc
  )
  pca = list_to_pprcomp(pca)

  verbose_message(message_text = "PCA finished!", verbose = verbose)

  rm(mass_matrix)
  gc()
  eigenvalues = pca$sdev ^ 2
  # Step 5: Compute Principal Components
  # Choose number of principal components, k
  # if not input, use scree test to help find retained components

  if (is.null(num_retained_component)) {
    if (!is.null(variance_explained_threshold)) {
      tryCatch({
        cumulative_variance = cumsum(eigenvalues) / sum(eigenvalues)
        threshold = variance_explained_threshold  # Example threshold

        if (show_variance_plot){

          par(mfrow = c(1, 1))
          par(mar = c(2, 2, 1, 1))
          # Plot cumulative proportion of variance explained
          plot(
            cumulative_variance,
            type = 'b',
            main = "Cumulative Variance Explained",
            xlab = "Number of Principal Components",
            ylab = "Cumulative Proportion of Variance Explained"
          )

          # Add a horizontal line at the desired threshold

          abline(h = threshold,
                 col = "red",
                 lty = 2)
        }

        # Find the number of principal components to retain based on the threshold
        retained =  which(cumulative_variance >= threshold)[1] - 1
      },
      error = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      },
      warning = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      })

    } else{
      # if threshold not inputted, use Kaiser's criterion
      verbose_message(message_text = "Both variance_explained_threshold and num_retained_component not inputted, use Kaiser's criterion for determination", verbose = verbose)


      plot(
        eigenvalues,
        type = 'b',
        main = "Scree Plot",
        xlab = "Principal Component",
        ylab = "Eigenvalue"
      )

      # Add a horizontal line at 1 (Kaiser's criterion)
      abline(h = 1,
             col = "red",
             lty = 2)

      # Add a vertical line at the elbow point
      elbow_point <- which(diff(eigenvalues) < 0)[1]
      abline(v = elbow_point,
             col = "blue",
             lty = 2)
      retained = length(which(eigenvalues >= 1))
    }
  } else{
    retained = as.integer(num_retained_component)
    if (is.na(retained)) {
      stop("Please enter correct number of principle components to retain")
    }
  }


  tryCatch({
    input_mz = data.frame(cbind(
      row_id = 1:ncol(resampled_mat),
      mz = as.numeric(stringr::str_extract(row.names(SpaMTP[[assay]]@features), pattern = "\\d+\\.?\\d*"))
    ))
  },
  error = function(cond) {
    stop(
      "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
    )
  },
  warning = function(cond) {
    stop(
      "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
    )
  })


  #### Load the Cleaned and summarized DB ####
  #load("data/Chebi_db.rda")
  #load("data/HMDB_db.rda")

  # Set the db that you want to search against
  db = rbind(HMDB_db, Chebi_db)
  # set which adducts you want to search for
  #load("data/adduct_file.rda")

  if (is.null(ion_mode)) {
    stop("Please enter correct polarity:'positive' or 'negative'")
  } else if (ion_mode == "positive") {
    test_add = sub(" ", "", adduct_file$adduct_name[which(adduct_file$charge >= 0)])
  } else if (ion_mode == "negative") {
    test_add = sub(" ", "", adduct_file$adduct_name[which(adduct_file$charge <= 0)])
  }
  # Using Chris' pipeline for annotation
  # 1) Filter DB by adduct.
  db_1 = db_adduct_filter(db, test_add, polarity = ifelse(ion_mode == "positive",
                                                          "pos", "neg"), verbose = verbose)

  # 2) only select natural elements
  db_2 = formula_filter(db_1)

  # 3) search db against mz df return results
  # Need to specify ppm error
  # If ppm_error not specified, use function to estimate
  # Set error tolerance
  ppm_error = 1e6 / tof_resolution / sqrt(2 * log(2))
  db_3 = proc_db(input_mz, db_2, ppm_error) %>% mutate(entry = stringr::str_split(Isomers,
                                                                                  pattern = "; "))

  verbose_message(message_text = "Query necessary data and establish pathway database" , verbose = verbose)

  input_id = lapply(db_3$entry, function(x) {
    x = unlist(x)
    index_hmdb = which(grepl(x, pattern = "HMDB"))
    x[index_hmdb] = paste0("hmdb:", x[index_hmdb])
    index_chebi = which(grepl(x, pattern = "CHEBI"))
    x[index_chebi] = tolower(x[index_chebi])
    return(x)
  })

  #load("/data/chem_props.rda")

  db_3 = db_3 %>% mutate(inputid = input_id)
  rampid = c()
  chem_source_id = unique(chem_props$chem_source_id)

  verbose_message(message_text = "Query db for addtional matching" , verbose = verbose)

  pb2 = txtProgressBar(
    min = 1,
    max = nrow(db_3),
    initial = 0,
    style = 3
  )
  for (i in 1:nrow(db_3)) {
    rampid[i] = (chem_props$ramp_id[which(chem_source_id %in% db_3$inputid[i][[1]])])[1]
    setTxtProgressBar(pb2, i)
  }
  close(pb2)
  db_3 = cbind(db_3, rampid)

  verbose_message(message_text = "Query finished!" , verbose = verbose)

  ####################################################################################################

  # get rank pathway database

  verbose_message(message_text = "Getting reference pathways...." , verbose = verbose)

  #load("data/analytehaspathway.rda")
  #load("data/pathway.rda")
  #load("data/source.rda")


  pathway_db = get_analytes_db(unlist(input_id), analytehaspathway,
                               chem_props, pathway)

  pathway_db = pathway_db[which(!duplicated(names(pathway_db)))]
  # get names for the ranks
  name_rank = lapply(input_mz$mz, function(x) {
    return(unique(na.omit(db_3$rampid[which(db_3$observed_mz == x)])))
  })

  #Set progress bar

  verbose_message(message_text = "Runing set enrichment analysis", verbose = verbose)
  pca_sea_list = list()
  pb_new = txtProgressBar(
    min = 0,
    max = retained,
    initial = 0,
    style = 3
  )
  for (i in 1:retained) {
    # get the absolute value and sign of the loading
    loading = data.frame(cbind(
      abs_loading = abs(pca[["rotation"]][, i]),
      sign_loading = sign(pca[["rotation"]][, i])
    )) %>% arrange(desc(abs_loading))
    # run set enrichment analysis

    ranks = unlist(lapply(1:length(pca[["rotation"]][, i]), function(x) {
      pc_new = rep(pca[["rotation"]][x, i], times = length(name_rank[[x]])) +
        sample(-5:5, length(name_rank[[x]]), replace = T) * 1e-7
      names(pc_new) = name_rank[[x]]
      return(pc_new)
    }))
    ranks = ranks[which(!duplicated(names(ranks)))]
    suppressWarnings({
      gsea_result = fgsea::fgsea(
        pathways =  pathway_db,
        stats = ranks,
        minSize = 5,
        maxSize = 500
      ) %>% filter(pval <= p_val_threshold) %>% mutate(principle_component = paste0("PC", i)) %>% mutate(leadingEdge_metabolites = lapply(leadingEdge, function(x) {
        temp = unique(unlist(x))
        metabolites_name = c()
        for (z in 1:length(temp)) {
          index <- find_index(name_rank , temp[z])
          direction = sign(sum(pca[["rotation"]][index, i]))
          metabolites_name = c(metabolites_name,
                               paste0(
                                 tolower(chem_props$common_name[which(chem_props$ramp_id == temp[z])])[1],
                                 ifelse(direction >= 0, "↑", "↓")
                               ))
        }

        return(metabolites_name)
      }))
    })
    setTxtProgressBar(pb_new, i)
    # Make sure sign of loading is positive to make it positively correlate with the PC
    pca_sea_list = rlist::list.append(pca_sea_list,
                                      gsea_result)
  }
  close(pb_new)
  names(pca_sea_list) = paste0("PC", 1:retained)

  ######################### Creating html
  image_matrix = matrix(
    Matrix::rowSums(resampled_mat),
    ncol = new.width,
    nrow = new.height,
    byrow = T
  )

  # # plot(mass_matrix_with_coord[,1],
  # #      mass_matrix_with_coord[,2])
  # image( matrix(
  #  (resampled_mat[,1]),
  #   ncol = new.width,
  #   nrow = new.height,
  #   byrow = T
  # ))
  # Assign different colours to different layers

  #image_matrix =  Matrix::rowSums(resampled_mat)

  quantiles <-
    quantile(as.numeric(as.vector(image_matrix)), probs = seq(0, 1, 0.2))
  quantiles <-
    unique(quantiles + seq(0, 1e-10, length.out = length(quantiles)))
  col_index <-
    cut(as.vector(image_matrix),
        breaks = quantiles,
        labels = FALSE)


  colors <- c("red", "blue", "green", "yellow", "orange")

  # (1) Tissue image
  quantiles = quantile(as.numeric(as.vector(image_matrix)),
                       probs = seq(0, 1, 0.2))
  quantiles <-
    unique(quantiles + seq(0, 1e-10, length.out = length(quantiles)))

  col_index <-
    cut(t(image_matrix), breaks = quantiles, labels = FALSE)

  matrix_data_melted = data.frame(cbind(x = reshape2::melt(t(image_matrix))[,1],
                                        y = reshape2::melt(t(image_matrix))[,2],
                                        as.vector(image_matrix))) %>%
    mutate(color = col_index)

  matrix_data_melted = matrix_data_melted %>% rowwise() %>% mutate(col_nam  = colors[color])

  #index = which(matrix_data_melted$image_matrix!= 0)

  tissue_image <- list(
    y = matrix_data_melted$y,
    x = matrix_data_melted$x,
    value = matrix_data_melted$image_matrix,
    colour = matrix_data_melted$col_nam,
    width = new.width,
    height = new.height
  )


  #a = reshape2::melt(c(1: nrow(pca$x)))
  # Convert to JSON

  # Optionally, write to a file



  #PCA_result
  # Convert list to JSON
  # Extract relevant components
  num = 5
  prcomp_data <- list(
    rotation = t(pca$rotation[, 1:num]),
    center = pca$center,
    scale = pca$scale,
    sdev = pca$sdev,
    x = t(pca$x[, 1:num]),
    mz =  row.names(SpaMTP@assays[[assay]]@features),
    coordinate = t(cbind(
      matrix_data_melted$x,
      matrix_data_melted$y
    )),
    name =  sprintf(
      "(%d,%d)",
      matrix_data_melted$x,
      matrix_data_melted$y
    )
  )

  ################### GSEA results
  var_gsea_table = c()
  for (i in 1:num) {
    temp = pca_sea_list[[i]]
    if (nrow(temp) == 0) {
      empty_row = t(c("No Signifcant Pathway Found",
                      rep(NA, times = ncol(temp) - 1)))
      temp = rbind(temp, empty_row, use.names = FALSE)
    }
    metabolites = c()
    if (nrow(temp) == 0) {
      var_gsea_table = c(var_gsea_table,
                         paste0('var PC', i, '=[', paste(rep('[],', times = 4),
                                                         collapse = ""), "]"))
      next
    }
    for (j in 1:length(temp$leadingEdge_metabolites)) {
      metabolites = c(metabolites,
                      paste0(unique(temp$leadingEdge_metabolites[[j]]), collapse = ", "))
    }
    temp_table = cbind(temp[, 1:7], metabolites = metabolites)
    var_gsea_table = c(
      var_gsea_table,
      paste0(
        'var PC',
        i,
        '=[',
        jsonlite::toJSON(temp_table$pathway),
        ',',
        jsonlite::toJSON(temp_table$pval),
        ',',
        jsonlite::toJSON(temp_table$NES),
        ',',
        jsonlite::toJSON(temp_table$metabolites),
        ']',
        collapse = ""
      )
    )
  }


  json_data <- jsonlite::toJSON(prcomp_data)
  json_matrix_data <- jsonlite::toJSON(tissue_image)

  html_temp = paste0(
    '<!DOCTYPE html>
  <html lang="en">
  <head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <script src="gmm_class.js"></script>
  <script src="kmeans.js"></script>
  <title>Prcomp Visualization</title>
  <style>
    .container {
      display: flex;
      justify-content: start; /* Optional: Adjust as needed */
    }
    .box {
      border: 1px solid black; /* Optional: Add borders for visibility */
      position: relative;
    }

    .input-container {
            display: inline-block;
        }

        /* Optional styling for label */
        label {
            margin-right: 5px; /* Adjust spacing between label and input */
        }
        .red-text {
    color: red;
    position: absolute;
    top: 0;
    left: 0;
}
  </style>
  </head>
  <body>
  <h1>Prcomp Data</h1>
  <div class="input-container"></div>
  <label>x-axis selection
    <select id="xSelect">
        <option value="" disabled selected>Select your PC</option>
        <option value="PC1" selected="selected">PC1</option>
        <option value="PC2">PC2</option>
        <option value="PC3">PC3</option>
        <option value="PC4">PC4</option>
        <option value="PC5">PC5</option>
    </select>
    </label>
    <label>y-axis selection
      <select id="ySelect">
          <option value="" disabled selected>Select your PC</option>
          <option value="PC1">PC1</option>
          <option value="PC2"selected="selected">PC2</option>
          <option value="PC3">PC3</option>
          <option value="PC4">PC4</option>
          <option value="PC5">PC5</option>
      </select>
      </label>
      <label for="integerInput">Enter an Integer for Cluster Numbers:</label>
      <input type="number" id="numcluster" name="integerInput" step="1" min = "2">
      <label>Clustering method
        <select id="clustering">
            <option value="" disabled selected>Select clustering method</option>
            <option value="GMM">GMM</option>
            <option value="Kmean">Kmean</option>
        </select>
        <button type="button" id= "button1">Regenerate Clusters</button>
        </label>
      </div>
    <div id = "Warning" class = "red-text"> </div>
    <div class="container">
      <div>
        <div id="pcaPlot" class="box" style="width: 600px; height: 400px;">
        </div>
        <div class="box" id="plot" style="position: relative;width: 600px; height: 400px;z-index: 2;">
          <canvas id="rasterCanvas" style="position: absolute; top: 0; left: 0;width:',600/new.width*98,'%; height: ',400/new.height*98,'%;z-index: 0;"></canvas>
          <canvas id="pointCanvas" style="position: absolute; top: 0; left: 0;width: ',600/new.width*98,'%; height: ',400/new.height*98,'%;z-index: 1;"></canvas>
        </div>
      </div>
      <div>
        <div id="table1" class="box" style="width: 600px; height: 400px;"></div>
        <div class="box" id="table2" style="width: 600px; height: 400px;"></div>
      </div>
    </div>

  <script>


  // Embed the JSON data directly \n',
paste0(var_gsea_table, collapse = "\n")
,
'
const prcompData = [',
json_data ,
'];
const matrixData = [',
json_matrix_data,
'];
// Example of using the data
// Example of using the data
//PCA plot

 var data_table1 = [{
  type: "table",
  header: {
    values: [["<b>Pathway</b>"], ["<b>p_value</b>"],
				 ["<b>NES</b>"], ["<b>metabolites</b>"]],
    align: "center",
    line: {width: 1, color: "black"},
    fill: {color: "grey"},
    font: {family: "Arial", size: 12, color: "white"}
  },
  cells: {
    values: PC1,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 11, color: ["black"]}
  }
}]

var layout_table1 = {
  title: "PC1 Pathway Information"
}

Plotly.newPlot("table1", data_table1, layout_table1);

// Table 2
var data_table2  = [{
  type: "table",
  header: {
    values: [["<b>Pathway</b>"], ["<b>p_value</b>"],
				 ["<b>NES</b>"], ["<b>metabolites</b>"]],
    align: "center",
    line: {width: 1, color: "black"},
    fill: {color: "grey"},
    font: {family: "Arial", size: 12, color: "white"}
  },
  cells: {
    values: PC2,
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 11, color: ["black"]}
  }
}]
var layout_table2 = {
  title: "PC2 Pathway Information"
}
Plotly.newPlot("table2", data_table2 ,layout_table2);


// pca PLOT
// Function to initialize parameters of GMM


function createArray(length) {
    // Calculate the value that each slot should hold
    const slotValue = 1 / length;

    // Create an array and populate it with the calculated value
    const array = Array(length).fill(slotValue);

    // Adjust the last element to ensure the sum equals 1
    array[length - 1] += 1 - array.reduce((sum, value) => sum + value, 0);

    return array;
}
function generateColors(n) {
    const colors = [];
    const increment = 360 / n;

    for (let i = 0; i < n; i++) {
        const hue = Math.round(increment * i);
        const saturation = 70; // You can adjust this value if needed
        const lightness = 70; // You can adjust this value if needed
        const color = `hsl(${hue}, ${saturation}%, ${lightness}%)`;
        colors.push(color);
    }

    return colors;
}
function transpose(array) {
    if (array.length === 0 || !Array.isArray(array[0])) return [];

    // Create an array to hold the result with the same number of sub-arrays as there are elements in each sub-array
    let result = [];
    for (let i = 0; i < array[0].length; i++) {
        result[i] = [];
        for (let j = 0; j < array.length; j++) {
            result[i].push(array[j][i]);
        }
    }

    return result;
}


// Assuming your data structure is stored in a variable called data_structure
var xCoordinates = matrixData[0].x;
var yCoordinates = matrixData[0].y;
var colors = matrixData[0].colour;
var value = matrixData[0].value
const canvas = document.getElementById("rasterCanvas");
const ctx = canvas.getContext("2d");
const container_tissue = document.getElementById("plot");
canvas.width = container_tissue.offsetWidth;
canvas.height = container_tissue.offsetHeight;
        // Calculate scaling factors
const scaleX = canvas.width/matrixData[0].height/1.1; // Assuming original width is 600
const scaleY = canvas.height/matrixData[0].width/1.1;
// draw raster image
  // Set canvas size


for (var i = 0; i < xCoordinates.length; i++) {
  //
            // Convert coordinates to canvas coordinates
            var canvasX = xCoordinates[i]  ;
            var canvasY = yCoordinates[i]  ;
            // Draw point
            ctx.fillStyle = colors[i];
            ctx.fillRect(canvasX, canvasY, 1, 1);
        }

function byRowToArrayByColumn(arr, numRows, numCols) {
    const result = [];
    for (let col = 0; col < numCols; col++) {
        for (let row = 0; row < numRows; row++) {
            result.push(arr[row * numCols + col]);
        }
    }
    return result;
}
function transposeArray(array) {
    return array[0].map((_, colIndex) => array.map(row => row[colIndex]));
}
var pc = [PC1,PC2,PC3,PC4,PC5]
console.log(pc)
function updatepcaPlot() {
  var num_to_cluster = Number(document.getElementById("numcluster").value)
  console.log(num_to_cluster)
  var candidate = ["PC1","PC2","PC3","PC4","PC5"]
  var xValue = document.getElementById("xSelect").value;
  var yValue = document.getElementById("ySelect").value;
  var new_pc1_index = candidate.findIndex(num => num === xValue)
  var new_pc2_index = candidate.findIndex(num => num === yValue)
  var method = document.getElementById("clustering").value;
  var pca_input = prcompData[0].x[new_pc1_index].map((element, index) => [element, prcompData[0].x[new_pc2_index][index]])
  //var pca_input = prcompData[0].x[0].map((element, index) => [element, prcompData[0].x[1][index]])


  //console.log(kmean_result)
  //console.log(kmean_result.centroids)

  let xMax = pca_input.reduce((max, subarray) => Math.max(max, subarray[0]), -Infinity);
  let yMax = pca_input.reduce((max, subarray) => Math.max(max, subarray[1]), -Infinity);
  let xMin = pca_input.reduce((min, subarray) => Math.min(min, subarray[0]), Infinity);
  let yMin = pca_input.reduce((min, subarray) => Math.min(min, subarray[1]), Infinity);
  let dx = xMax-xMin;
	let dy = yMax-yMin;
	let covariances = Array(num_to_cluster).fill(0)
		.map(_ => [[dx*dx*.01, 0], [0, dy*dy*.01]]);
  let colors = generateColors(num_to_cluster)
  if(method === "GMM"){
    let gmm = new GMM({
	  weights: Array(num_to_cluster).fill(1/num_to_cluster),
	  means: Array(num_to_cluster).fill(0).map(_ => [xMin + Math.random()*dx, yMin + Math.random()*dy]),
	  covariances: covariances
});
  //console.log(pca_input)
  pca_input.forEach(p => gmm.addPoint(p));
  gmm.runEM(10);
  let assignments = []

  for (let i = 0; i < pca_input.length; i++) {
    let probNorm = gmm.predictNormalize(pca_input[i]);
    let index = probNorm.indexOf(Math.max(...probNorm));
    assignments.push(colors[index])
    // Do something with probNorm if needed
 }
  //console.log(assignments)
  //prcompData[0].x[0].map((element, index) => [element, prcompData[0].x[1][index]])
  var update_pca = {
    "marker.color": [assignments]
  };
  //const canvas = document.querySelector("cycles_canvas");
  //const draw = new Draw(canvas, xMin, xMax, yMin, yMax);
  Plotly.restyle("pcaPlot", update_pca, 0);
  //for(let i=0; i<gmm.clusters; i++) {
	//		draw.ellipse(gmm.means[i], gmm.covariances[i], clusterColors[i]);
	//	}

  // update colours
  //console.log(assignments)
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  //let new_colour = assignments
  for (var i = 0; i < yCoordinates.length; i++) {
  //
            // Convert coordinates to canvas coordinates
            var canvasX = xCoordinates[i]  ;
            var canvasY = yCoordinates[i]  ;
            // Draw point
            ctx.fillStyle = assignments[i];
            ctx.fillRect(canvasX, canvasY, 1, 1);
        }

  }else if(method === "Kmean"){
    let assignments = []
    let kmean_result = kmeans(pca_input, num_to_cluster)
    console.log(kmean_result)
    for (let i = 0; i < pca_input.length; i++) {
    assignments.push(colors[kmean_result[1][i]])
    // Do something with probNorm if needed
 }
    // Do something with probNorm if needed
    //console.log(assignments)

 var update_pca = {
    "marker.color": [assignments]
  };
  //const canvas = document.querySelector("cycles_canvas");
  //const draw = new Draw(canvas, xMin, xMax, yMin, yMax);
  Plotly.restyle("pcaPlot", update_pca, 0);

  ctx.clearRect(0, 0, canvas.width, canvas.height);

  for (var i = 0; i < yCoordinates.length; i++) {
  //
            // Convert coordinates to canvas coordinates
            var canvasX = xCoordinates[i]  ;
            var canvasY = yCoordinates[i]  ;
            // Draw point
            ctx.fillStyle = assignments[i];
            ctx.fillRect(canvasX, canvasY, 1, 1);
        }
  }
}
document.getElementById("clustering").addEventListener("change", updatepcaPlot);
var button = document.getElementById("button1")

button.onclick = function(){
  updatepcaPlot();
}
  // Simulated PCA data (you would use your actual PCA data here)
  var pcaData = {
    names: prcompData[0].name,
    pc1: prcompData[0].x[0],
    pc2: prcompData[0].x[1]
  };

  var trace1pc = {
    x: pcaData.pc1,
    y: pcaData.pc2,
    mode: "markers",
    type: "scatter",
    marker: { size: 0.3 },
    text: pcaData.names,
    hoverinfo: "text+x+y"
  };

  var data_pc = [trace1pc];

  var layout_pc = {
    title: "PCA Plot",
    xaxis: { title: "PC1" },
    yaxis: { title: "PC2" },
    width: 600,
    height: 400
  };

  Plotly.newPlot("pcaPlot", data_pc, layout_pc);

  function updatePlot() {
      var candidate = ["PC1","PC2","PC3","PC4","PC5"]
      var xValue = document.getElementById("xSelect").value;
      var yValue = document.getElementById("ySelect").value;

      var new_pc1_index = candidate.findIndex(num => num === xValue)
      var new_pc2_index = candidate.findIndex(num => num === yValue)
      var pcaData = {
        names: prcompData[0].name,
        pc1: prcompData[0].x[new_pc1_index],
        pc2: prcompData[0].x[new_pc2_index]
      };

      var trace1pc = {
        x: pcaData.pc1,
        y: pcaData.pc2,
        mode: "markers",
        type: "scatter",
        marker: { size: 0.3 },
        text: pcaData.names,
        hoverinfo: "text+x+y"
      };
      var layout_pc = {
      title: "PCA Plot",
      xaxis: { title: candidate[new_pc1_index] },
      yaxis: { title: candidate[new_pc2_index] },
      width: 600,
      height: 400
     };
      var data_pc = [trace1pc];

      Plotly.newPlot("pcaPlot", data_pc, layout_pc);
    // Update table with new content
    var data_table1_new = [{
  type: "table",
  header: {
    values: [["<b>Pathway</b>"], ["<b>p_value</b>"],
				 ["<b>NES</b>"], ["<b>metabolites</b>"]],
    align: "center",
    line: {width: 1, color: "black"},
    fill: {color: "grey"},
    font: {family: "Arial", size: 12, color: "white"}
  },
  cells: {
    values: pc[new_pc1_index],
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 11, color: ["black"]}
  }
}]

var layout_table1_new = {
  title: "PC"+(new_pc1_index+1)+" Pathway Information"
}

Plotly.react("table1", data_table1_new, layout_table1_new);
Plotly.relayout("table1", layout_table1_new);
// Table 2
var data_table2_new  = [{
  type: "table",
  header: {
    values: [["<b>Pathway</b>"], ["<b>p_value</b>"],
				 ["<b>NES</b>"], ["<b>metabolites</b>"]],
    align: "center",
    line: {width: 1, color: "black"},
    fill: {color: "grey"},
    font: {family: "Arial", size: 12, color: "white"}
  },
  cells: {
    values: pc[new_pc2_index],
    align: "center",
    line: {color: "black", width: 1},
    font: {family: "Arial", size: 11, color: ["black"]}
  }
}]
var layout_table2_new = {
  title: "PC"+(new_pc2_index+1)+" Pathway Information"
}
Plotly.react("table2", data_table2_new ,layout_table2_new);
Plotly.relayout("table2", layout_table2_new);
}



    // Event listeners to update the plot when selections change
    document.getElementById("xSelect").addEventListener("change", updatePlot);
    document.getElementById("ySelect").addEventListener("change", updatePlot);
    document.getElementById("clustering").addEventListener("change", updatepcaPlot);

// Example of using the data
console.log(prcompData[0].rotation[0]);






  // Plot each pixel

  // let xCoords = Array.from({length: 100}, (_, i) => 2*i);
  // let yCoords = Array.from({length: 100}, (_, i) => i);
  // let colors2 = new Array(10000).fill().map(() => `rgb(${Math.random()*255}, ${Math.random()*255}, ${Math.random()*255})`);

  // // Set canvas size
  // canvas.width = 100;
  // canvas.height = 100;
  // for (let x = 0; x < xCoords.length; x++) {
  //     for (let y = 0; y < yCoords.length; y++) {
  //         // Index for the colors array
  //         let index = x * yCoords.length + y;
  //         ctx.fillStyle = colors2[index];
  //         ctx.fillRect(xCoords[x], yCoords[y], 1, 1);
  //     }
  // }
// Create a trace for the heatmap
var trace = {
  y: xCoordinates,
  x: yCoordinates,
  z: value,
  mode: "markers",
  type: "heatmap",
  color: colors,
  //marker: {
  //  color:colors,
  //  size: 0.2
  //},
  //colorscale: "Viridis"
};

let lastX = null;
let lastY = null;
const canvas2 = document.getElementById("pointCanvas");
const ctx2 = canvas2.getContext("2d");
canvas2.width = container_tissue.offsetWidth;
canvas2.height = container_tissue.offsetHeight;
// Create the plot

//Plotly.newPlot("plot", [trace,trace_heatmap_scatter], layout,{staticPlot: true});
document.getElementById("pcaPlot").on("plotly_hover", function(data){
    // Get index of clicked point
    var pointIndex = data.points[0].pointIndex;
    // console.log(pointIndex)
    let x_new_scatter = xCoordinates[pointIndex]
    let y_new_scatter = yCoordinates[pointIndex]

        ctx2.clearRect(0, 0, canvas2.width, canvas2.height);

        // Draw new point
        ctx2.fillStyle = "#42f5c2";  // Color of the new point
        ctx2.beginPath();
        ctx2.arc(x_new_scatter, y_new_scatter , 3, 0, Math.PI * 2, true);  // Small circle for the point
        console.log(x_new_scatter)
        ctx2.fill();
        // Update last point coordinates
        lastX = x_new_scatter;
        lastY = y_new_scatter;
});

document.addEventListener("DOMContentLoaded", function() {
    var divElement = document.getElementById("Warning");
    var newText = document.createElement("span"); // Creates a new <span> element
    newText.className = "Warning"; // Assigns the class for styling
    newText.textContent = "Note that ↑ and ↓ of PC pathway information are the sign of the loading corresponds to the metabolites m/z, since the sign is abitarily assigned during PCA depends on the setting, the arrow is only comparable within each PC and do not have biological implications, they are just used to describe how the metabolite contribute the PC"; // Adds text

    divElement.style.position = "relative"; // Ensures that the div can hold absolute positioned elements
    divElement.appendChild(newText); // Adds the newly created <span> to the div
});


</script>
</body>
</html>')
  setwd(path)
  gmm = readLines(system.file("js_data_files", "gmm_class.js", package = "SpaMTP"))
  writeLines(gmm,"gmm_class.js")
  kmean = readLines(system.file("js_data_files", "kmeans.js", package = "SpaMTP"))
  writeLines(kmean,"kmeans.js")

  writeLines(html_temp,"pca_analysis.html")
  return(list(pathway_enrichment_pc = pca_sea_list,
              new.width = as.integer(width/as.numeric(resampling_factor)),
              new.height = as.integer(height/as.numeric(resampling_factor))))
}









#' Helper function for building a pathway db based on detected metabolites
#'
#' @param input_id Vector of characters defining the detected metabolites.
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

######################################## HELPER FUNCTIONS ########################################

#' Creates a pprcomp object based on an input list
#'
#' @param lst List containing PCA results
#'
#' @return A pprcomp object contating results from PCA analysis
#'
#' @examples
#' #HELPER FUNCTION
list_to_pprcomp <- function(lst) {
  # Create an empty object with class pprcomp
  obj <- structure(list(), class = "prcomp")
  # Assign components from the list to the object
  obj$sdev <- lst$sdev
  obj$rotation <- lst$rotation
  obj$center <- lst$center
  obj$scale <- lst$scale
  obj$x <- lst$x
  # Add other components as needed

  # Return the constructed pprcomp object
  return(obj)
}


#' Finds the index values of the m/z values with their respective GSEA result
#'
#' @param lst List containing relative mz analytes and pathways
#' @param value Value returned based on the GSEA results
#'
#' @return returns a vector of indices that match the relative GSEA results to the m/z list
#'
#' @examples
#' #HELPER FUNCTION
find_index <- function(lst, value) {
  indices <- which(sapply(lst, function(x) value %in% x))
  if (length(indices) == 0) {
    return(NULL)  # If value not found, return NULL
  } else {
    return(indices)
  }
}


##########################################################################################################################################################
### Bellow helper functions are sourced from https://stackoverflow.com/questions/11123152/function-for-resizing-matrices-in-r by Vyga.


#' Helper function to rescale a sampled matrix
#'
#' @param x Vector defining new matrix coordinates
#' @param newrange Vector defining range of old coordinates
#'
#' @return Vector containing rescaled coordinates
#'
#' @examples
#' #HELPER FUNCTION
rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

#' Helper function to resize a matrix back to its original layout after sampling
#'
#' @param mat matrix defining the sampled object
#' @param ndim number of dimentions to resize the sampled matrix to (default = dim(mat)).
#'
#' @return Returns a resized sampled matrix to match the dimentions of the original
#'
#' @examples
#' #HELPER FUNCTION
ResizeMat <- function(mat, ndim=dim(mat)){
  #if(!require(fields)) stop("`fields` required.")

  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)

  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)

  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))

  # interpolation
  ans[ncord] <- fields::interp.surface(obj, loc)

  ans
}




#' Helper function that generated PCA analysis results for a SpaMTP Seurat Object
#'
#' @param SpaMTP SpaMTP Seurat class object that contains spatial metabolic information.
#' @param num_retained_component is an integer value to indicated preferred number of PCs to retain
#' @param variance_explained_threshold Numeric value defining the explained variance threshold.
#' @param resampling_factor is a numerical value > 0, indicate how you want to resample the size of original matrix.
#' @param p_val_threshold is the p value threshold for pathways to be significant.
#' @param byrow is a boolean to indicates whether each column of the matrix is built byrow or bycol.
#' @param assay Character string defining the SpaMTP assay to extract intensity values from.
#' @param slot Character string defining the assay slot containing the intensity values.
#' @param flip_plot Boolean defining whether to rotate the plot 90 degrees.
#' @param show_variance_plot Boolean indicating weather to display the variance plot output by this analysis.
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed.
#'
#'
#' @return PCA object contating embeddings and projections
#'
#' @import dplyr
#'
#' @examples
#' # HELPER FUNCTION
getPCA <- function(SpaMTP,
                   num_retained_component,
                   variance_explained_threshold,
                   resampling_factor,
                   byrow,
                   assay,
                   slot,
                   flip_plot,
                   show_variance_plot,
                   verbose) {
  # PCA analysis
  verbose_message(message_text = "Scaling original matrix", verbose = verbose)

  mass_matrix = Matrix::t(SpaMTP[[assay]]@layers[[slot]])

  coords <- GetTissueCoordinates(SpaMTP)[c("x", "y")]

  if (flip_plot){
    colnames(coords) <- c("y", "x")
  }


  mass_matrix_with_coord = cbind(coords,
                                 as.matrix(mass_matrix))
  width = max(coords[c("y")]) #- min(coords[c("y")])
  height = max(coords[c("x")]) # - min(coords[c("x")])

  if (!is.null(resampling_factor) & resampling_factor!= 1) {
    verbose_message(message_text = "Running matrix resampling...." , verbose = verbose)
    pb = txtProgressBar(
      min = 0,
      max = ncol(mass_matrix),
      initial = 0,
      style = 3
    )
    if (!is.numeric(resampling_factor)) {
      stop("Please enter correct resampling_factor")
    }
    new.width = as.integer(width / resampling_factor)
    new.height = as.integer(height / resampling_factor)

    resampled_mat = matrix(nrow =  new.height * new.width)
    for (i in 1:ncol(mass_matrix)) {
      temp_mz_matrix = matrix(mass_matrix[, i],
                              ncol = width,
                              nrow = height,
                              byrow = byrow)
      resampled_temp = ResizeMat(temp_mz_matrix, c(new.width,
                                                   new.height))
      resampled_mat = cbind(resampled_mat, as.vector(resampled_temp))
      setTxtProgressBar(pb, i)
    }
    close(pb)
    resampled_mat = resampled_mat[,-1]
    colnames(resampled_mat) = colnames(mass_matrix)

    verbose_message(message_text = "Resampling finished!" , verbose = verbose)

    gc()
  }else{
    resampled_mat = mass_matrix
    new.width = as.integer(width)
    new.height = as.integer(height)
  }

  verbose_message(message_text = "Running the principal component analysis (can take some time)" , verbose = verbose)

  # Runing PCA

  resampled_mat_standardised = as.matrix(Matrix::t(
    Matrix::t(resampled_mat) - Matrix::colSums(resampled_mat) / nrow(resampled_mat)
  ))

  verbose_message(message_text = "Computing the covariance" , verbose = verbose)
  cov_mat <- cov(resampled_mat_standardised)

  verbose_message(message_text = "Computing the eigenvalue/eigenvectors", verbose = verbose)
  eigen_result <- eigen(cov_mat)
  gc()
  # Extract eigenvectors and eigenvalues
  eigenvectors <- eigen_result$vectors
  eigenvalues <- eigen_result$values

  verbose_message(message_text = "Computing PCA", verbose = verbose)

  pc = pbapply::pblapply(1:ncol(resampled_mat_standardised), function(i) {
    temp = resampled_mat_standardised[, 1] * eigenvectors[1, i]
    for (j in 2:ncol(resampled_mat_standardised)) {
      temp = temp + resampled_mat_standardised[, j] * eigenvectors[j, i]
    }
    return(temp)
  })
  pc = do.call(cbind, pc)
  colnames(pc) = paste0("PC", 1:ncol(eigenvectors))
  # make pca object
  colnames(eigenvectors) = paste0("PC", 1:ncol(eigenvectors))
  rownames(eigenvectors) = colnames(resampled_mat)
  pca = list(
    sdev = sqrt(eigenvalues),
    rotation = eigenvectors,
    center = Matrix::colSums(resampled_mat) / nrow(resampled_mat),
    scale = FALSE,
    x = pc
  )
  pca = list_to_pprcomp(pca)

  verbose_message(message_text = "PCA finished!", verbose = verbose)

  rm(mass_matrix)
  gc()
  eigenvalues = pca$sdev ^ 2
  # Step 5: Compute Principal Components
  # Choose number of principal components, k
  # if not input, use scree test to help find retained components

  if (is.null(num_retained_component)) {
    if (!is.null(variance_explained_threshold)) {
      tryCatch({
        cumulative_variance = cumsum(eigenvalues) / sum(eigenvalues)
        threshold = variance_explained_threshold  # Example threshold

        if (show_variance_plot){

          par(mfrow = c(1, 1))
          par(mar = c(2, 2, 1, 1))
          # Plot cumulative proportion of variance explained
          plot(
            cumulative_variance,
            type = 'b',
            main = "Cumulative Variance Explained",
            xlab = "Number of Principal Components",
            ylab = "Cumulative Proportion of Variance Explained"
          )

          # Add a horizontal line at the desired threshold

          abline(h = threshold,
                 col = "red",
                 lty = 2)
        }

        # Find the number of principal components to retain based on the threshold
        retained =  which(cumulative_variance >= threshold)[1] - 1
      },
      error = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      },
      warning = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      })

    } else{
      # if threshold not inputted, use Kaiser's criterion
      verbose_message(message_text = "Both variance_explained_threshold and num_retained_component not inputted, use Kaiser's criterion for determination", verbose = verbose)


      plot(
        eigenvalues,
        type = 'b',
        main = "Scree Plot",
        xlab = "Principal Component",
        ylab = "Eigenvalue"
      )

      # Add a horizontal line at 1 (Kaiser's criterion)
      abline(h = 1,
             col = "red",
             lty = 2)

      # Add a vertical line at the elbow point
      elbow_point <- which(diff(eigenvalues) < 0)[1]
      abline(v = elbow_point,
             col = "blue",
             lty = 2)
      retained = length(which(eigenvalues >= 1))
    }
  } else{
    retained = as.integer(num_retained_component)
    if (is.na(retained)) {
      stop("Please enter correct number of principle components to retain")
    }
  }
  return(pca)
}



#' Generates PCA analysis results for a SpaMTP Seurat Object
#'
#' @param SpaMTP SpaMTP Seurat class object that contains spatial metabolic information.
#' @param num_retained_component is an integer value to indicated preferred number of PCs to retain
#' @param variance_explained_threshold Numeric value defining the explained variance threshold (default = 0.9).
#' @param resampling_factor is a numerical value > 0, indicate how you want to resample the size of original matrix (default = 1).
#' @param p_val_threshold is the p value threshold for pathways to be significant (default = 0.05).
#' @param byrow is a boolean to indicates whether each column of the matrix is built byrow or bycol (default = FALSE).
#' @param assay Character string defining the SpaMTP assay to extract intensity values from (default = "SPM").
#' @param slot Character string defining the assay slot containing the intensity values (default = "counts").
#' @param flip_plot Boolean defining whether to rotate the plot 90 degrees (default = FALSE).
#' @param show_variance_plot Boolean indicating weather to display the variance plot output by this analysis (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#'
#' @return SpaMTP object with pca results stored in the
#' @export
#'
#' @examples
#' # HELPER FUNCTION
RunMetabolicPCA <- function(SpaMTP,
                   num_retained_component = NULL,
                   variance_explained_threshold = 0.9,
                   resampling_factor = 1,
                   byrow = FALSE,
                   assay = "SPM",
                   slot = "counts",
                   flip_plot = FALSE,
                   show_variance_plot= FALSE,
                   verbose = TRUE)
{
  verbose_message(message_text = "Running PCA Analysis ... ", verbose = verbose)

  pca <- getPCA(SpaMTP = SpaMTP,
                num_retained_component = num_retained_component,
                variance_explained_threshold = variance_explained_threshold,
                resampling_factor = resampling_factor,
                byrow = byrow,
                assay = assay,
                slot = slot,
                flip_plot = flip_plot,
                show_variance_plot= show_variance_plot,
                verbose = verbose)

  SpaMTP_pca <- pca

  rownames(SpaMTP_pca$rotation) <- rownames(SpaMTP[[assay]]@features)
  rownames(SpaMTP_pca$x) <- rownames(SpaMTP@meta.data)

  SpaMTP_pcas <- SeuratObject::CreateDimReducObject(embeddings = SpaMTP_pca$x, loadings = SpaMTP_pca$rotation, assay = assay, key = "pca_")

  SpaMTP[["pca"]] <- SpaMTP_pcas

  return(SpaMTP)
}






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



#' Plotting of differentially expressed metabolic pathways per grouping identity
#'
#' @param SpaMTP A seurat object contains spatial metabolomics/transcriptomics features or both.
#' @param ident Character string defining the metadata column used for grouping. NOTE: Metadata column must be a factor.
#' @param polarity The polarity of the MALDI experiment. Inputs must be either NULL, 'positive' or 'negative'. If NULL, pathway analysis will run in neutral mode (default = NULL).
#' @param assay Character string defining the SpaMTP assay to extract intensity values from (default = "SPM").
#' @param slot Character string specifying which slot to pull the m/z inensity values from (default = "counts").
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion `[ion mass,m/z]`/`[Full width at half height]` (default = 30000).
#' @param ppm_error is the parts-per-million error tolerance of matching m/z value with potential metabolites (default = NULL).
#' @param pval_cutoff Numeric value defining the p-value cutoff that defines significant pathways (default = 0.01).
#' @param plot.sig Numeric value defining the p-value cutoff that highlights the significant pathways on the plot (default = 0.05).
#' @param verbose Boolean
#'
#' @return A list consists of a ggplot2 object and a dataframe containing the set enrichment results
#' @export
#'
#' @examples
#' # PathwaysPerRegion(SpaMTP, ident = "clusters", polarity = "positive")
PathwaysPerRegion = function(SpaMTP,
                      ident,
                      polarity = NULL,
                      assay = "Spatial",
                      slot = "counts",
                      tof_resolution =30000,
                      ppm_error = NULL,
                      pval_cutoff = 0.01,
                      plot.sig = 0.05,
                      verbose = TRUE
){
  #(1) Get the rank entry for each cluster
  assignment = SpaMTP@meta.data[[ident]]
  cluster = levels(assignment)
  mass_matrix = Matrix::t(SpaMTP[[assay]]@layers[[slot]])

  # (2) Annotation
  tryCatch({
    input_mz = data.frame(cbind(
      row_id = 1:ncol(mass_matrix),
      mz = as.numeric(stringr::str_extract(row.names(SpaMTP[[assay]]@features), pattern = "\\d+\\.?\\d*"))
    ))
  },
  error = function(cond) {
    stop(
      "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
    )
  },
  warning = function(cond) {
    stop(
      "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
    )
  })
  # Set the db that you want to search against
  db = rbind(HMDB_db, Chebi_db)
  # set which adducts you want to search for
  #load("data/adduct_file.rda")

  if (is.null(polarity)) {
    stop("Please enter correct ion_mode:ion_mode = 'positive' or ion_mode ='negative'")
  } else if (polarity == "positive") {
    test_add = sub(" ", "", adduct_file$adduct_name[which(adduct_file$charge >= 0)])
  } else if (polarity == "negative") {
    test_add = sub(" ", "", adduct_file$adduct_name[which(adduct_file$charge <= 0)])
  }
  # Using Chris' pipeline for annotation
  # 1) Filter DB by adduct.
  db_1 = db_adduct_filter(db, test_add, polarity = ifelse(polarity== "positive",
                                                          "pos", "neg"), verbose = verbose)

  # 2) only select natural elements
  db_2 = formula_filter(db_1)

  # 3) search db against mz df return results
  # Need to specify ppm error
  # If ppm_error not specified, use function to estimate
  # Set error tolerance
  ppm_error = ppm_error %||% (1e6 / tof_resolution / sqrt(2 * log(2)))
  db_3 = proc_db(input_mz, db_2, ppm_error) %>% mutate(entry = stringr::str_split(Isomers,
                                                                                  pattern = "; "))

  verbose_message(message_text = "Query necessary data and establish pathway database" , verbose = verbose)

  input_id = lapply(db_3$entry, function(x) {
    x = unlist(x)
    index_hmdb = which(grepl(x, pattern = "HMDB"))
    x[index_hmdb] = paste0("hmdb:", x[index_hmdb])
    index_chebi = which(grepl(x, pattern = "CHEBI"))
    x[index_chebi] = tolower(x[index_chebi])
    return(x)
  })

  #load("/data/chem_props.rda")

  db_3 = db_3 %>% mutate(inputid = input_id)
  rampid = c()
  chem_source_id = unique(chem_props$chem_source_id)

  verbose_message(message_text = "Query db for addtional matching" , verbose = verbose)

  pb2 = txtProgressBar(
    min = 1,
    max = nrow(db_3),
    initial = 0,
    style = 3
  )
  for (i in 1:nrow(db_3)) {
    rampid[i] = (chem_props$ramp_id[which(chem_source_id %in% db_3$inputid[i][[1]])])[1]
    setTxtProgressBar(pb2, i)
  }
  close(pb2)
  db_3 = cbind(db_3, rampid)

  verbose_message(message_text = "Query finished!" , verbose = verbose)
  # get names for the ranks
  name_rank = lapply(input_mz$mz, function(x) {
    return(unique(na.omit(db_3[which(db_3$observed_mz == x),])))
  })


  name_rank_ids = lapply(input_mz$mz, function(x) {
    return(unique(na.omit(db_3$rampid[which(db_3$observed_mz == x)])))
  })
  # Get pathway db
  verbose_message(message_text = "Constructing pathway database ..." , verbose = verbose)
  pathway_db = get_analytes_db(unlist(input_id), analytehaspathway,
                               chem_props, pathway)

  pathway_db = pathway_db[which(!duplicated(names(pathway_db)))]
  verbose_message(message_text = "Run set enrichment for each cluster" , verbose = verbose)
  pb3 = txtProgressBar(
    min = 0,
    max = length(cluster),
    initial = 0,
    style = 3
  )
  gsea_all_cluster = data.frame()
  for(i in 1:length(cluster)){
    # Get coordinates for the elements in the cluster
    clu_wise = which(assignment  == cluster[i])
    sub_mass_mat = Matrix::colSums(mass_matrix[clu_wise,])/length(clu_wise)
    other_region = Matrix::colSums(mass_matrix[-clu_wise,])/(nrow(mass_matrix) - length(clu_wise))
    log_fc = scale(log(sub_mass_mat/other_region+0.1), center = 0)

    ranks = unlist(lapply(1:length(log_fc), function(x) {
      pc_new = rep(log_fc[x], times = length(name_rank[[x]]$rampid)) +
        sample(-5:5, length(name_rank[[x]]$rampid), replace = T) * 1e-7
      names(pc_new) = name_rank[[x]]$rampid
      return(pc_new)
    }))

    ranks = ranks[which(!duplicated(names(ranks)))]
    suppressWarnings({
      gsea_result = fgsea::fgsea(
        pathways =  pathway_db,
        stats = ranks,
        minSize = 5,
        maxSize = 500
      )  %>% mutate(Cluster_id =  paste0(cluster[i])) %>% mutate(leadingEdge_metabolites = lapply(leadingEdge, function(x) {
        temp = unique(unlist(x))
        metabolites_name = c()
        for (z in 1:length(temp)) {
          index <- find_index(name_rank_ids, temp[z])
          direction = log_fc
          metabolites_name = c(metabolites_name,
                               paste0(
                                 tolower(chem_props$common_name[which(chem_props$ramp_id == temp[z])])[1],
                                 ifelse(direction >= 0, "↑", "↓")
                               ))
        }
        return(metabolites_name)
      }))
    })
    gsea_all_cluster = rbind(gsea_all_cluster,gsea_result)
    setTxtProgressBar(pb3, i)
  }
  close(pb3)
  gsea_p_val <- (pval_cutoff %||% 0.01)
  plot.sig <- (plot.sig %||% 0.05)
  gsea_all_cluster_sig =   gsea_all_cluster %>% filter(pval<=gsea_p_val) %>% mutate(Significance = ifelse(pval <= plot.sig,
                                                                                                                       paste0("> ", plot.sig),
                                                                                                                       "Non-Significant"))
  colnames(gsea_all_cluster_sig )[1] = "pathwayName"
  gsea_all_cluster_sig = merge(gsea_all_cluster_sig,
                               pathway, by = "pathwayName")

  ########################################################
  #Plot#
  #dendrogram based on jaccard distance
  gsea_all_cluster_sig = gsea_all_cluster_sig %>% mutate(pathnameid = paste0(pathwayName,
                                                                             "(",sourceId,")"))
  pathwaynames = unique(gsea_all_cluster_sig$pathnameid)
  n = length(unique(gsea_all_cluster_sig$pathnameid))
  jaccard_matrix = matrix(nrow = n,
                          ncol = n)
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
    scale_x_continuous(expand = c(0, 0.1),
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
  #### ggplot
  gg_dot = ggplot(data = gsea_all_cluster_sig, aes(x = factor(Cluster_id,
                                                              levels = sort(unique(Cluster_id))
  ),
  y = factor(pathnameid,levels = ggdendro::dendro_data(hc)$labels$label))
  ) +
    geom_point(aes(colour = as.numeric(NES), size = as.numeric(size))) +
    scale_colour_gradient2(
      name = "Normalised enrichment score",
      low = "blue",   # Original low color
      high = "red",   # Original high color
      limits = c(min(gsea_all_cluster_sig$NES), max(gsea_all_cluster_sig$NES))
    ) +
    scale_size_continuous(name = "Number of altered metabolites \n in the pathway") +
    ggnewscale::new_scale_colour() +
    geom_point(shape = 1, aes(colour = Significance, size = as.numeric(size)+0.2)) +
    scale_color_manual(
      values = c(
        setNames("black", paste0("> ", plot.sig)),
        "Non-Significant" = NA
      )
    ) +
    labs(title = "",
         y = "Pathways", x = "Idents") +
    theme(
      title = element_text(size = 8, face = 'bold'),
      axis.text.x = element_text(
        size = 12,
        angle = 90,
        vjust = 0.5,
        hjust = 1,
      ),
      #20
      axis.title = element_text(size = 9),
      #40
      legend.key.size = unit(1, 'cm'),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    ) +
    ggnewscale::new_scale_colour() + theme(panel.background = element_rect(fill = "white"),
                               panel.grid = element_line(color = "grey")) +   theme(
                                 legend.position = "left",
                                 axis.text.y = element_blank())

  combined_plot = gg_dot| plt_dendr


  return(list(ggplot_item = combined_plot,
              pathway_df = gsea_all_cluster_sig))
}













