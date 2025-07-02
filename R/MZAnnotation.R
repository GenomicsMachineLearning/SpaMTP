#### SpaMTP m/z Annotation Functions #####################################################################################################################################################################################

#' Subset a SpaMTP Seurat Spatial Metabolomic object by a list of m/z's
#'
#' @param data A Seurat Spatial Metabolomic Object for subsetting.
#' @param features A list of character strings defining the features/mz values to subset against.
#' @param assay A character string identifying the Seurat assay which contains the count data being subset.
#'
#' @returns A subset Seurat object containing only m/z values that were specified
#' @export
#'
#' @examples
#' # SubsetMZFeatures(SeuratObj, c("mz-160","mz-170","mz-180"))
SubsetMZFeatures <- function(data, features, assay = "Spatial"){
  feature.metadata <- data[[assay]]@meta.data
  sub.data <- subset(data, features = features)

  sub.data[[assay]]@meta.data <- feature.metadata[feature.metadata[["mz_names"]] %in% features,]
  return(sub.data)
}



#' Filters the annotation list to only include the first n number of annotations per m/z
#'
#' @param annotation_column Vector of the meta.data column containing the m/z annotations.
#' @param n Numeric value defining the number of annotations to keep (default = 3).
#'
#' @return Vector containing the first n number of annotations
#'
#' @examples
#' # labels_to_show(`SeuratObject[["Spatial"]]@meta.data$annotations`, n = 3)
labels_to_show <- function(annotation_column, n = 3) {
  new_column <- sapply(strsplit(annotation_column, "; "), function(x) {
    # Filter out NA values and select the first three entries
    non_na_values <- x[!is.na(x)]
    paste(non_na_values[1:min(n, length(non_na_values))], collapse = "; ")
  })
  return(new_column)
}




#' Annotates m/z values stored in a SpaMTP Object
#'
#' This function assigns each valid m/z peak with one/multiple metabolite names based on the mass difference between the observed value and the theoretical value documented in the reference database.
#' SpaMTP contains 4 cleaned reference databases to choose from these include HMDB, Lipid Maps, ChEBI and GNPS. These databases can also be combined for increased coverage.
#'
#' @param data Seurat Spatial Metabolomic Object containing m/z values for annotation.
#' @param db Reference metabolite dataset in the form of a Data.Frame.
#' @param assay Character string defining the Seurat assay which contains the mz counts being annotated (default = "Spatial").
#' @param raw.mz.column Character string defining the Seurat assay slot which contains the raw mz values, this is without the 'mz-' and are a vector of integers. This is setup by default when running the cardinal_to_seurat() function (default = "raw_mz").
#' @param ppm_error Numeric value indicating the size of the ppm error allowed when matching molecular weights between Seurat object and reference dataset. If only want exact matches set ppm = 0 (default = NULL).
#' @param adducts List of adducts to use for searching the database (e.g. "M+NH4","M+Na","M+CH3OH+H","M+K" etc.). For all possible adducts please visit [here](https://github.com/GenomicsMachineLearning/SpaMTP/blob/main/R/MZAnnotation.R#L305). If NULL will take the full list of adducts (default = NULL).
#' @param polarity Character string defining the polarity of adducts to use, either "positive", "negative" or "neutral" (default = "positive").
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion `[ion mass,m/z]`/`[Full width at half height]`. This value is used to estimate ppm_error when set to NULL (default = 30000).
#' @param filepath Character string of the directory to store the _annotated_mz_peaks.csv. If set to NULL no dataframe will be saved (default = NULL).
#' @param return.only.annotated Boolean value indicating if the annotated Seurat Object should only include m/z values that were successfully annotated (default = TRUE).
#' @param save.intermediate Boolean indicating whether to save an intermediate file in the `@tools` slot of the SpaMTP object required for later analysis functions such as `FindRegionalPathways()` (default = TRUE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A Seurat Object with m/z values annotated. These annotations are stored in the relative assay's meta.data (e.g. SeuratObj`[["Spatial"]][[]]`)
#' @export
#'
#' @examples
#' # HMDB_db <- load("data/HMDB_1_names.rds")
#' # Annotated_SeuratObj <- AnnotateSM(SeuratObj, HMDB_db)
AnnotateSM <- function(data, db, assay = "Spatial", raw.mz.column = "raw_mz", ppm_error = NULL, adducts = NULL, polarity = "positive", tof_resolution = 30000, filepath = NULL, return.only.annotated = TRUE, save.intermediate = TRUE, verbose = TRUE){

  if (is.null(data@assays[[assay]])) {
    stop(paste0("No assay '",assay,"'exists in SpaMTP object! Please check assay name input ..."))
  }

  ## Extracting m/z values from SpaMTP@assay@meta.data
  mz_df <- data[[assay]][[raw.mz.column]]
  mz_df$mz <- mz_df[[raw.mz.column]]
  mz_df$row_id <- seq(1, length(mz_df[[raw.mz.column]]))
  mz_df <- mz_df[c("row_id", "mz")]

  db_3 <- annotateTable(mz_df= mz_df, db = db, ppm_error = ppm_error, adducts = adducts, polarity = polarity,tof_resolution = tof_resolution,verbose = verbose)


  if (save.intermediate){
    data@tools$db_3 <- db_3
  }

  if (!(is.null(filepath))){
    path <- paste0(filepath,"_annotated_mz_peaks.csv")
    message(paste0("Writing annotated m/z file to: '", path, "'"))
    data.table::fwrite(db_3, path)
  }

  verbose_message(message_text = "Adding annotations to Seurat Object .... ", verbose = verbose)

  result_df <- db_3 %>%
    dplyr::group_by(observed_mz) %>%
    dplyr::summarise(
      all_IsomerNames = paste(IsomerNames, collapse = "; "),
      all_Isomers = paste(Isomers, collapse = "; "),
      all_Isomers_IDs = paste(Isomers_IDs, collapse = "; "),
      all_Adducts = paste(unique(Adduct), collapse = "; "),
      all_Formulas = paste(unique(Formula), collapse = "; "),
      all_Errors = paste(round(Error,4), collapse = "; ")
    )
  rownames(result_df) <- paste0("mz-",result_df$observed_mz)

  result_df$mz_names <- rownames(result_df)

  data[[assay]][["mz_names"]] <- rownames(data[[assay]][[]])


  result_df <- result_df %>% dplyr::mutate(present = TRUE)


  # Perform left join and replace NAs with "No Annotation"
  feature_metadata <- data[[assay]][[]] %>%
    dplyr::left_join(result_df, by = "mz_names") %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ifelse(is.na(.), "No Annotation", .)))

  feature_metadata <- dplyr::select(feature_metadata, -present)

  data[[assay]][[]] <- feature_metadata
  if (return.only.annotated == TRUE){

    verbose_message(message_text = "Returning Seurat object that include ONLY SUCCESSFULLY ANNOTATED m/z features", verbose = verbose)

    if(length(Assays(data)) != 1){
      features = c()
      for(non_met_assay in Assays(data)[which(Assays(data)!=assay)]){
        features =  c(features, rownames(data@assays[[non_met_assay]]@features))
      }
      data <- suppressWarnings({SubsetMZFeatures(data, assay = assay, features = c(result_df$mz_names, features))})
    }else{
      data <- suppressWarnings({SubsetMZFeatures(data, assay = assay, features = result_df$mz_names)})
    }
  }

  return(data)
}




#' Annotates m/z values sotred in a data.frame based on reference metabolite dataset
#'
#' Helper function for `AnnotatesSM()` and `FishersPathwayAnalysis()`.
#'
#' @param mz_df dataframe containing m/z values for annotation.
#' @param db Reference metabolite dataset in the form of a Data.Frame.
#' @param ppm_error Numeric value indicating the size of the ppm error allowed when matching molecular weights between Seurat object and reference dataset. If only want exact matches set ppm = 0 (default = NULL).
#' @param adducts List of adducts to use for searching the database (e.g. "M+NH4","M+Na","M+CH3OH+H","M+K" etc.). For all possible adducts please visit [here](https://github.com/GenomicsMachineLearning/SpaMTP/blob/main/R/MZAnnotation.R#L305). If NULL will take the full list of adducts (default = NULL).
#' @param polarity Character string defining the polarity of adducts to use, either "positive", "negative" or "neutral" (default = "positive").
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion `[ion mass,m/z]`/`[Full width at half height]`. This value is used to estimate ppm_error when set to NULL (default = 30000).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns Generates an intermediate annotated m/z dataframe
#'
#' @importFrom rlang %||%
#'
#' @examples
#'
#' ### HelperFunction
annotateTable <- function(mz_df, db, ppm_error = NULL, adducts = NULL, polarity = "positive", tof_resolution = 30000, verbose = TRUE){


  # Uses:
  # db: db that you want to search against

  # test_add_pos: which adducts you want to search for
  # Note; "M+NH4","M+Na","M+CH3OH+H","M+K" etc. Look at the formula filter func to get the rest of the possible adducts.

  # ppm_error: the ppm error/threshold for searching

  # Three main steps relates to the three main functions
  # Steps 1) & 2) are aimed at condensing the databases by applying 1) a filter to only consider the adducts that the user specifies. 2) Filtering the molecular formulas to contain only elements that the user specifies. # Step 3) This last function then does the database matching and searching.
  # 1) Filter DB by adduct.
  verbose_message(message_text = paste0("Filtering provided database by ", paste0(adducts, collapse = ", "), " adduct/s"), verbose = verbose)

  if (polarity == "positive") {
    test_add_pos = adduct_file$adduct_name[which(adduct_file$charge > 0)]
    test_add_pos <- gsub(" ", "", test_add_pos)
    test_add_pos = test_add_pos[which(test_add_pos %in% (adducts %||% test_add_pos))]
    # Using Chris' pipeline for annotation
    # 1) Filter DB by adduct.
    db_1 <- db_adduct_filter(db,
                            test_add_pos,
                            polarity = "pos",
                            verbose = verbose)
  } else if (polarity == "negative") {
    test_add_neg = adduct_file$adduct_name[which(adduct_file$charge < 0)]
    test_add_neg <- gsub(" ", "", test_add_neg)
    test_add_neg = test_add_neg[which(test_add_neg %in% (adducts %||% test_add_neg))]
    # Using Chris' pipeline for annotation
    # 1) Filter DB by adduct.
    db_1 <- db_adduct_filter(db,
                            test_add_neg,
                            polarity = "neg",
                            verbose = verbose)
  } else if (polarity == "neutral") {
    # Using Chris' pipeline for annotation
    # 1) Filter DB by adduct.
    db_1 <- db %>% mutate("M" = `M-H ` + 1.007276)
  } else{
    stop("Please enter correct polarity from: 'positive', 'negative', 'neutral'")
  }


  # 2) only select natural elements
  db_2 <- formula_filter(db_1)

  # 3) search db against mz df return results
  verbose_message(message_text = "Searching database against input m/z's to return annotaiton results", verbose = verbose)

  if (is.null(ppm_error) && is.null(tof_resolution)){
    stop("ppm_error and tof_resolution cannot both = NULL! Please set one of these variables to calculate the ppm threshold for annotations ... ")
  }

  ppm_error = ppm_error %||% (1e6 / tof_resolution / sqrt(2 * log(2)))

  db_3 <- proc_db(mz_df, db_2, ppm_error)

  ## Add in database labels
  db_3 = db_3 %>% dplyr::mutate(entry = stringr::str_split(Isomers, pattern = "; "))
  input_id = lapply(db_3$entry, function(x) {
    x = unlist(x)
    index_hmdb = which(grepl(x, pattern = "HMDB"))
    x[index_hmdb] = paste0("hmdb:", x[index_hmdb])
    index_chebi = which(grepl(x, pattern = "CHEBI"))
    x[index_chebi] = tolower(x[index_chebi])
    index_lipidm = which(grepl(x, pattern = "^LM"))
    x[index_lipidm] = paste0("LIPIDMAPS:", x[index_lipidm])
    return(x)
  })
  db_3$entry <- NULL

  db_3$Isomers_IDs <- input_id

  db_3 <- db_3 %>%
    dplyr::mutate(Isomers_IDs = sapply(Isomers_IDs, function(x) stringr::str_c(x, collapse = "; ")))

  return(db_3)
}



#' Refines and reduces m/z annotations
#'
#' Used to subset dataset to only include annotations that have n number of entries.
#' For example some peaks can have multiple annotations. Peaks which have above n number of annotations assigned will be removed from Seurat Object.
#'
#' @param obj Seurat object needing annotation refinement. This object must have annotations present in 'obj`[[assay]]@meta.data`'
#' @param assay Character string defining the Seurat object assay where the annotation data is stored (default = "Spatial").
#' @param n Integer defining the number of entries an annotation can have assigned. Any higher counts will be removed (default = 1).
#'
#' @return Refined Seurat object that only contains annotated mz values that have n number of annotations assigned (per mz value)
#' @export
#'
#' @examples
#' # HMDB_db <- load("data/HMDB_1_names.rds")
#' # AnnotatedSeuratObj <- AnnotateSeuratMALDI(SeuratObj, HMDB_db)
#'
#' # getRefinedAnnotations(AnnotatedSeuratObj, n = 2)
getRefinedAnnotations <- function(obj, assay = "Spatial",n = 1){
  n <- n-1
  subset.metadata <- obj[[assay]]@meta.data %>% dplyr::filter(stringr::str_count(all_IsomerNames, ";") %in% c(0:n))
  subset.obj <- SubsetMZFeatures(data = obj, features = subset.metadata$mz_names,assay = assay)
  return(subset.obj)

}

########################################################################################################################################################################################################################


#### SpaMTP Find m/z Annotation Functions #####################################################################################################################################################################################

#' Adds in backslashes required to take into account special using grepl such as brackets and +
#'
#' @param input_string Character string requiring additional backslashes
#'
#' @return Character string containing double backslashes around special features
add_backslashes_to_specialfeatures <- function(input_string) {
  # Use gsub to replace ( with \\( and ) with \\)
  result_string <- gsub("\\(", "\\\\(", input_string)
  result_string <- gsub("\\)", "\\\\)", result_string)
  result_string <- gsub("\\[", "\\\\[", result_string)
  result_string <- gsub("\\]", "\\\\]", result_string)
  result_string <- gsub("\\+", "\\\\+", result_string)
  result_string <- gsub(" ", "[ ]", result_string)

  return(result_string)
}


#' Find Annotation
#'
#' Searches through annotated m/z values to return all which contain the metabolite search term provided
#'
#' @param data Seurat Spatial Metabolomic Object containing annotated m/z values.
#' @param metabolite Character string of metabolite search term.
#' @param assay Character string defining the Seurat assay that contains the annotated metadata corresponding to the m/z values (default = "Spatial").
#' @param search.exact Boolean value defining if to only return m/z values which contain the exact match to the metabolite search term (default = FALSE).
#' @param column.name Character string defining the column name where the annotations are stored in the slot meta.data (default = "all_IsomerNames").
#'
#' @return A Data.Frame containing the peak metadata corresponding to the metabolite search term provided
#' @export
#'
#' @examples
#' # SearchAnnotations(SeuratObj, "Glucose", search.exact = TRUE)
SearchAnnotations <- function (data, metabolite, assay = "Spatial",search.exact = FALSE, column.name = "all_IsomerNames"){

  ## Takes into account '( )' in the string name
  search_term <- add_backslashes_to_specialfeatures(metabolite)

  indexs <- which(grepl(search_term, data[[assay]]@meta.data[column.name][[1]], ignore.case = TRUE))

  if (search.exact){
    df <- data[[assay]]@meta.data[indexs,]
    indexs <- c()
    for (row in rownames(df)){
      annotation_list <- df[row,][[column.name]]
      split_list <- unlist(strsplit(annotation_list, "; "))

      if (metabolite %in% split_list){
        indexs <- c(indexs, row)
      }

    }
  }

  return(data[[assay]]@meta.data[indexs,])
}



#' Finds if any metabolite is duplicated across multiple m/z values.
#'
#' @param data Seurat Spatial Metabolomic Object containing annotated m/z values.
#' @param assay Character string defining the Seurat assay that contains the annotated metadata corresponding to the m/z values (default = "Spatial").
#'
#' @return Vector of character strings describing metabolites that are assigned to multiple m/z values
#' @export
#'
#' @examples
#' # FindDuplicateAnnotations(SeuratObj)
FindDuplicateAnnotations <- function (data, assay = "Spatial"){
  all_annotations <- data[[assay]]@meta.data$all_IsomerNames
  all_terms <- unlist(strsplit(all_annotations, ";"))
  all_terms <- trimws(all_terms)
  terms_counts <- table(all_terms)
  return(names(terms_counts[terms_counts > 1]))
}

#' Gets values from a single metadata column for a respective m/z value.
#'
#' @param obj SpaMTP Spatial Metabolomic Seurat Object containing annotated m/z values.
#' @param mz Character string specifying the m/z value to return the metadata for.
#' @param assay Character string defining the Seurat assay that contains the annotated metadata corresponding to the m/z values (default = "Spatial").
#' @param metadata.column Character string corresponding to the `@meta.data` column to extract the data from (default = "all_IsomerNames").
#' @param separate Boolean indicating whether to separate the metadata string and return a vector. Note, if `TRUE` the metadata column should contain values separate by "; " (default = TRUE).
#'
#' @return Vector of character strings containing the metadata for a specific mz value.
#' @export
#'
#' @examples
#  ##### Example for getting metabolite annotation for a m/z value
#' # GetMZMetadata(SpaMTP, mz = "mz-100", metadata.column = "all_IsomerNames")
#'
#' ##### Example for getting metabolite annotation IDs for a m/z value
#' # GetMZMetadata(SpaMTP, mz = "mz-100", metadata.column = "all_Isomers")
GetMZMetadata <- function(obj, mz, assay = "Spatial", metadata.column = "all_IsomerNames", separate = TRUE){
  df <- obj[[assay]]@meta.data
  mz_row <- df[df$mz_names == mz,metadata.column]
  if (separate){
    return(strsplit(mz_row ,split = "; ")[[1]])
  } else {
    return(mz_row)
  }
}

########################################################################################################################################################################################################################

#!! ALL CODE BELOW WAS WRITEN BY Christopher Fitzgerald github https://github.com/ChemCharles !!#


#' Checks if the complete adduct is in the data base, else returns a truncated adduct
#'
#' @param adduct Character string defining the adduct to be checked.
#' @param db DataFrame of the current reference database.
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return Character string of the complete or truncated adduct
#'
#' @examples
#' # HMDB_db <- load("data/HMDB_1_names.rds")
#' # check_and_truncate_adduct_vector(c("M+H"), HMDB_db)
check_and_truncate_adduct_vector <- function(adduct, db, verbose = TRUE) {
  element_exists <- adduct %in% colnames(db)
  missing_elements <- adduct[!element_exists]
  if (length(missing_elements) > 0) {
    for (missing_element in missing_elements) {
      verbose_message(message_text =  paste0("Adduct ",
                      missing_element,
                      " is not in the DB, it has been removed from the search."), verbose = verbose)

    }
    truncated_adduct <- adduct[element_exists]
    return(truncated_adduct)
  } else {
    return(adduct)
  }
}




#' Filters the provided metabolomic database by polarity and adducts
#'
#' @param adduct Character string defining the adduct to be checked.
#' @param db DataFrame of the current reference database.
#' @param polarity Character string defining the polarity of the adducts (default = "neg").
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A filtered reference metabolomic database DataFrame
#'
#' @examples
#' # HMDB_db <- load("data/HMDB_1_names.rds")
#' # db_adduct_filter(HMDB_db,c("M+H"), polarity = "pos")
db_adduct_filter <- function(db, adduct, polarity = "neg", verbose = TRUE) {
  # only include adducts from either neg or pos polarity
  if (polarity == "neg") {
    neg_adducts_1 <-
      c(
        "M-H2O-H",
        "M-H",
        "M+Na-2H",
        "M+Cl",
        "M+K-2H",
        "M+FA-H",
        "M+Hac-H",
        "M+Br",
        "M+TFA-H",
        "2M-H",
        "2M+FA-H",
        "2M+Hac-H",
        "3M-H"
      )
    pol <- neg_adducts_1
  }
  else if (polarity == "pos") {
    pos_adducts_1 <-
      c(
        "M+H",
        "M+NH4",
        "M+Na",
        "M+CH3OH+H",
        "M+K",
        "M+ACN+H",
        "M+2Na-H",
        "M+IsoProp+H",
        "M+ACN+Na",
        "M+2K+H",
        "M+DMSO+H",
        "M+2ACN+H",
        "M+IsoProp+Na+H",
        "2M+H",
        "2M+NH4",
        "2M+Na",
        "2M+K",
        "2M+ACN+H",
        "2M+ACN+Na"
      )
    pol <- pos_adducts_1
  } else{
    stop("Invalid polarity. Choose 'neg' or 'pos'.")
  }

  # get rid of spaces in the adduct names
  # in col names
  db <- db %>%
    dplyr::rename_all( ~ gsub(" ", "", .))

  # Filter the db by polarity
  db <- db %>%
    dplyr::select(formula, exactmass, isomers, isomers_inchikey, isomers_names, pol)

  # in adduct entry
  adduct <- gsub(" ", "", adduct)

  adduct <- check_and_truncate_adduct_vector(adduct, db, verbose = verbose)

  db_filtered <- db %>%
    dplyr::select(formula,
           exactmass,
           isomers,
           isomers_inchikey,
           isomers_names,
           tidyr::any_of(adduct)) %>%
    as.data.frame()
  return(db_filtered)
}




#' Checks if a formula contains only the allowed elements
#'
#' @param formula Character string defining t
#' @param allowed_elements Vector of character strings defining allowed elements
#'
#' @examples
#' ### Helper function ###
is_formula_valid <- function(formula,allowed_elements) {
  elements <-
    stringr::str_extract_all(formula, "[A-Z][a-z]*")[[1]] # defines an element as a Uppercase immediately followed by none or more lowercase letters.

  all(elements %in% allowed_elements) # Then checks if they are in the vector of allowed elements.
}




#' Filters reference Database to only select natural elements
#'
#' @param df DataFrame of the reference database.
#' @param elements Vector of character strings of elements to include (default = c("H", "C", "N", "O", "S", "Cl", "Br", "F", "Na", "P", "I", "Si")).
#'
#' @return A refined DataFrame which only includes annotations containing the specified elements
#'
#' @examples
#' ## Get filtered DB by adduct
#' # db_1 <- db_adduct_filter(db, test_add_pos, polarity = "pos")
#'
#' ## Refine to only include natural elements
#' # db_2 <- formula_filter(db_1)
formula_filter <- function(df, elements = NULL) {
  if (is.null(elements)) {
    elements <- c("H", "C", "N", "O", "S", "Cl", "Br",
                  "F", "Na", "P", "I","Si")
  }

  # Elements to allow
  allowed_elements <- elements

  # Filter rows based on allowed elements
  filtered_df <- df %>%
    dplyr::filter(sapply(formula, allowed_elements = allowed_elements,is_formula_valid))

  return(filtered_df)

}




#' Calculates the mz range of the observed_df
#'
#' @param input_df DataFrame of the observed dataframe being annotated
#'
#' @return A list containing the lower and upper mz range for the provided sample
#'
#' @examples
#' # mz_df <- SeuratObject`[["Spatial"]][["mz"]]`
#' # mz_df$row_id <- seq(1, length(mz_df$mz))
#'
#' # mass_range <- calculate_bounds(mz_df)
#' # lower_bound <- mass_range$lower_bound
#' # upper_bound <- mass_range$upper_bound
calculate_bounds <- function(input_df) {

  lower_bound <- min(input_df$mz, na.rm = TRUE)
  upper_bound <- max(input_df$mz, na.rm = TRUE)

  bounds <-
    list(lower_bound = lower_bound, upper_bound = upper_bound)

  return(bounds)
}





#' Calculates the ppm error as a valve
#'
#' @param observed_mz Numeric value defining the observed mz value.
#' @param reference_mz Numeric value defining the reference mz value.
#' @param ppm Numeric value defining the maximum acceptable ppm_error/threshold for searching.
#'
#' @return Numeric value defining the ppm_error between the observed and reference mz value
#'
#' @examples
#' ### Helper Function ###
ppm_error <- function(observed_mz, reference_mz, ppm) {
  abs_diff_ppm <-
    abs(observed_mz - reference_mz) / abs(reference_mz) * 1e6
  if (abs_diff_ppm <= ppm) {
    return(abs_diff_ppm)
  }
  else{
    return("Out")
  }
}




#' Calculates the ppm range and check if mz values are within the range
#'    -  Returns TRUE if match is found and false if no match.
#'
#' @param observed_mz Numeric value defining the observed mz value.
#' @param reference_mz Numeric value defining the reference mz value.
#' @param ppm Numeric value defining the maximum acceptable ppm_error/threshold for searching.
#'
#' @return Boolean value indicating if a match is found (TRUE) or not (FALSE)
#'
#' @examples
#' ### Helper Function ###
ppm_range_match <- function(observed_mz, reference_mz, ppm) {
  abs_diff_ppm <-
    abs(observed_mz - reference_mz) / abs(reference_mz) * 1e6
  abs_diff_ppm <= ppm
}



#' Searches observed mz values against the data base list and returns matching annotations
#'
#' @param observed_df DataFrame containing the observed mz values
#' @param reference_df DataFrame contating the reference mz values and relative annotations.
#' @param ppm_threshold Numeric value defining the maximum acceptable ppm_error/threshold allowed between observed and reference mz values
#'
#' @return A DataFrame containing matched mz values between the observed and reference dataframes
#'
#' @examples
#' # HMDB_db <- load("data/HMDB_1_names.rds")
#' # mz_df <- SeuratObject[["Spatial"]][["mz"]]
#' # mz_df$row_id <- seq(1, length(mz_df$mz))
#'
#' ## 1) Filter DB by adduct.
#' # db_1 <- db_adduct_filter(HMDB_db, c("M+H"), polarity = "pos")
#'
#' ## 2) only select natural elements
#' # db_2 <- formula_filter(db_1)
#'
#' ## 3) search db against mz df return results
#' # db_3 <- proc_db(mz_df, db_2, ppm_threshold = 5)
proc_db <- function(observed_df,
                    reference_df,
                    ppm_threshold = 10) {
  # create an empty list to store matched results.
  result_list <- list()

  # Check if observed_df has only one row
  if (nrow(observed_df) == 1) {
    # Create a dummy entry with all values set to 0
    dummy_row <-
      data.frame(row_id = 0, mz = 0)  # Modify this line based on your column names

    # Combine the dummy entry with the original observed_df
    observed_df <- rbind(observed_df, dummy_row)
  } else if (nrow(observed_df) == 0) {
    # Handle the case when observed_df is empty
    verbose_message(message_text =  "Warning: No entries detected in input mz list.\n Please check format of input list,\n it must contain row_id and mz as the headers", verbose = verbose)

    return(NULL)
  }

  # extract out bounds
  lower_bound <- calculate_bounds(observed_df)$lower_bound
  upper_bound <- calculate_bounds(observed_df)$upper_bound

  #  For loop to go through each column of the reference_df that is provided.
  #  probably a good idea to filter reference_df to only adducts that you want before putting it into this function.
  # the -c(1:4) essentially makes it loop over the numeric portions of the DBs.

  for (col_name in names(reference_df)[-c(1:5)]) {
    result_col <- list()
    for (i in seq_len(nrow(observed_df))) {
      # Check if the reference mz is within the range
      within_range <-
        reference_df[[col_name]] >= lower_bound &
        reference_df[[col_name]] <= upper_bound

      # Condition 1: Only proceed if there are values within the range
      if (!any(within_range)) {
        next
      }

      # check for matches
      matches <- ppm_range_match(
        observed_mz = observed_df$mz[i],
        reference_mz = reference_df[[col_name]],
        ppm = ppm_threshold
      )

      # Condition 2: skip to the next iteration if no match
      if (!any(matches)) {
        next
      }

      # Condition 3: index the matches in the reference df
      # which keeps things that are TRUE.
      # This allows referencing to be quicker below.
      matching_indices <- which(matches)
      # print(matching_indices)

      # Calculate ppm error for each match.
      error <- mapply(
        ppm_error,
        observed_mz = observed_df$mz[i],
        reference_mz = reference_df[[col_name]][matching_indices],
        ppm = ppm_threshold
      )

      # Extract out relevant info for each match.
      filtered_matches <- data.frame(
        ID = observed_df$row_id[i],
        Match = matches[matching_indices],
        observed_mz = observed_df$mz[i],
        Reference_mz = reference_df[[col_name]][matching_indices],
        Error = error,
        Adduct = col_name,
        # DB_ID = reference_df$compound_id[matching_indices],
        Formula = reference_df$formula[matching_indices],
        Exactmass = reference_df$exactmass[matching_indices],
        Isomers = reference_df$isomers[matching_indices],
        InchiKeys = reference_df$isomers_inchikey[matching_indices],
        IsomerNames = reference_df$isomers_names[matching_indices]
      )
      # store each results df for each mz value in observed_df
      result_col[[i]] <- filtered_matches
    }
    # store results df in a list for each mz value for each adduct
    result_list[[col_name]] <- result_col
  }

  # Combine the individual matches per adduct df from the list into a dataframe
  combined_results <-
    do.call(rbind, lapply(result_list, function(result_col)
      do.call(rbind, result_col)))
}



########################################################################################################################################################################################################################

#' Assign custom annotations to m/z values
#'
#' Adds custom metabolite annotations to respective m/z values (ideal for specific matrices such as FMP10).
#'
#' @param data SpaMTP Seurat object containing m/z intensity values.
#' @param annotations data.frame containing two columns named 'annotation' and 'mass'. These columns should contain the custom metabolite annotation and the relative m/z mass respectively.
#' @param assay Character string defining the Seurat object assay to store the respective annotations in the feature meta.data dataframe (default = "Spatial").
#' @param return.only.annotated Boolean defining whether to return a SpaMTP Seurat object containing only successfully annotated m/z values (default = FALSE).
#' @param mass.threshold Numeric value defining the acceptable threshold (plus-minus) between the custom annotations and the actual m/z values contained within the SpaMTP object (default = 0.05).
#' @param annotation.column Character string defining the feature meta.data column name that will contain the assigned annotations (default = "all_IsomerNames").
#'
#' @return SpaMTP Seurat object containing the custom annotations stored in the feature metadata dataframe.
#' @export
#'
#' @examples
#' # annotated_data <- AddCustomMZAnnotations(SpaMTP.obj, annotation.df)
AddCustomMZAnnotations <- function(data, annotations, assay = "Spatial", return.only.annotated = FALSE, mass.threshold = 0.05, annotation.column = "all_IsomerNames"){

  if (!all(c("annotation", "mass") %in% colnames(annotations))) {
    stop("Error: The annotation columns provided does not match the required format. Must bet 'annotation' and 'mass'")
  }

  true_mzs <- c()
  for (mass in annotations$mass){
    true_mz <- FindNearestMZ(data, mass, assay = assay)
    true_mzs <- c(true_mzs,true_mz)
  }

  annotations$mz_names <- true_mzs
  annotations$true_mzs <- gsub("mz-", "", annotations$mz_names)
  annotations$ppm_diff <- abs(annotations$mass - as.numeric(annotations$true_mzs))

  if (!is.null(mass.threshold)){
    annotations <- annotations %>% dplyr::filter(ppm_diff < mass.threshold)
  }

  annotations <- annotations %>% dplyr::rename(!!annotation.column := annotation)
  annotations <- annotations %>%
    dplyr::group_by(mz_names) %>%
    dplyr::summarise(
      across(
        everything(),  # Apply to all columns
        ~ paste(unique(.), collapse = "; ")  # For each column, collapse unique values into a single string
      )
    )


  data[[assay]]@meta.data <- data[[assay]]@meta.data %>%
    dplyr::left_join(annotations, by = "mz_names") %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, "No Annotation")))


  if (return.only.annotated == TRUE){

    annotated_mzs <- data[[assay]]@meta.data$mz_names[data[[assay]]@meta.data[[annotation.column]] != "No Annotation"]
    #verbose_message(message_text = "Returning Seurat object that include ONLY SUCCESSFULLY ANNOTATED m/z features", verbose = verbose)

    if(length(Assays(object = data)) != 1){
      features = c()
      for(non_met_assay in Assays(object = data)[which(Assays(object = data)!=assay)]){
        features =  c(features, rownames(data@assays[[non_met_assay]]@features))
      }
      data <- suppressWarnings({SubsetMZFeatures(data, assay = assay, features = c(annotated_mzs, features))})
    }else{
      data <- suppressWarnings({SubsetMZFeatures(data, assay = assay, features = annotated_mzs)})
    }
  }

  return(data)
}


#' Annotates FMP10 matrix data
#'
#' Adds metabolite annotations to respective m/z values generated using an FMP10 matrix. This is done based on a curated FMP10 matrix database.
#'
#' @param obj SpaMTP Seurat object containing m/z intensity values for annotation. Data should be generated with a FMP10 matrix.
#' @param only.fmp.adduct Boolean indicating if only metabolites with FMP10+ adducts (`+FMP10`, +`2FMP10`, etc.) should be assigned to m/z values (default = FALSE).
#' @param add.custom.annotation data.frame containing addition FMP10 metabolite annotations that are not in the current FMP10 database. Note: this data.frame must contain these column c("mass", "annotation", "Adduct", "Formula", "Isomers", "Isomers_IDs"). If set to NULL, only reference FMP10 database will be used (default = NULL).
#' @param assay Character string defining the Seurat object assay to store the respective annotations in the feature meta.data dataframe (default = "Spatial").
#' @param return.only.annotated Boolean defining whether to return a SpaMTP Seurat object containing only successfully annotated m/z values (default = FALSE).
#' @param mass.threshold Numeric value defining the acceptable threshold (plus-minus) between the custom annotations and the actual m/z values contained within the SpaMTP object (default = 0.05).
#' @param annotation.column Character string defining the feature meta.data column name that will contain the assigned annotations (default = "all_IsomerNames").
#'
#' @return SpaMTP Seurat object containing the relative metabolite annotations stored in the feature metadata dataframe.
#' @export
#'
#' @examples
#' # AddFMP10Annotations(spamtp, only.fmp.adduct = FALSE)
AddFMP10Annotations <- function(obj,  only.fmp.adduct = FALSE,
                                add.custom.annotation = NULL,
                                assay = "Spatial",
                                return.only.annotated = FALSE,
                                mass.threshold = 0.05,
                                annotation.column = "all_IsomerNames"){

  if(only.fmp.adduct){
    filtered_fmp10 <- filtered_fmp10[grepl(x = filtered_fmp10$Adduct, pattern = "FMP10"),]
  }

  if(!is.null(add.custom.annotation)){
    required_columns <- c("mass", "annotation", "Adduct", "Formula", "Isomers", "Isomers_IDs")

    if (all(required_columns %in% colnames(add.custom.annotation))) {
      add.custom.annotation$IsomerNames <- add.custom.annotation$annotation
      add.custom.annotation <- add.custom.annotation[c(required_columns, "IsomerNames")]
      filtered_fmp10 <- dplyr::bind_rows(filtered_fmp10, add.custom.annotation)
    } else {
      stop("data.frame provided in `add.custom.annotation` does not contain the correct column names. Column names must be: c('mass', 'annotation', 'Adduct', 'Formula', 'Isomers', 'Isomers_IDs')")
    }

  }



  obj <- AddCustomMZAnnotations(obj, filtered_fmp10,
                                assay = assay,
                                return.only.annotated = return.only.annotated,
                                annotation.column = annotation.column,
                                mass.threshold = mass.threshold)


  # Set up @tools$db_3
  mz_list <- lapply(obj[[assay]]@meta.data$mass, function(x) unlist(strsplit(x, "; ")))
  error_list <-  lapply(obj[[assay]]@meta.data$error, function(x) unlist(strsplit(x, "; ")))

  names(mz_list) <- obj[[assay]]@meta.data$raw_mz

  db_3 <- filtered_fmp10
  db_3$observed_mz <- ""

  # Iterate through the list to find matches and update the new column
  for (mz_name in names(mz_list)) {
    for (value in mz_list[[mz_name]]) {
      # Find matching values in the observed_mz column
      matches <- which(db_3$mass == value)
      # If there are matches, append the mz name to the new column
      if (length(matches) > 0) {
        db_3$observed_mz[matches] <- paste(db_3$observed_mz[matches], mz_name, sep = "; ")
        db_3$observed_mz[matches] <- paste(db_3$observed_mz[matches], mz_name, sep = "; ")

      }
    }
  }

  # Clean up the new_column by removing the leading "; "
  db_3$observed_mz <- sub("^; ", "", db_3$observed_mz)

  db_3 = db_3 %>%
    tidyr::separate_rows(observed_mz, sep = "; ")

  db_3 <- db_3[db_3$observed_mz != "",]

  #Add to obj
  obj@tools$db_3 <- db_3

  obj[[assay]]@meta.data$observed_mz <- obj[[assay]]@meta.data$mass
  obj[[assay]]@meta.data$all_Isomers <- obj[[assay]]@meta.data$Isomers
  obj[[assay]]@meta.data$all_Isomers_IDs <- obj[[assay]]@meta.data$Isomers_IDs
  obj[[assay]]@meta.data$all_Adducts <- obj[[assay]]@meta.data$Adduct
  obj[[assay]]@meta.data$all_Formulas <- obj[[assay]]@meta.data$Formula
  obj[[assay]]@meta.data$all_Errors <- obj[[assay]]@meta.data$ppm_diff

  obj[[assay]]@meta.data[c("mass", "Isomers", "Isomers_IDs", "Adduct","Formula", "ppm_diff", "IsomerNames", "observed_mz")] <- NULL

  return(obj)
}




#' Annotates vector of m/z values
#'
#' This function assigns each valid m/z peak with one/multiple metabolite names based on the mass difference between the observed value and the theoretical value documented in the reference database.
#' This function is to be used when dealing with large datasets as a preprocessing step. Users can annotate m/z values first and then subset their data accordinly before loading it into a SpaMTP Seurat Object.
#'
#' @param mzs Vector containing m/z values for annotation.
#' @param db Reference metabolite dataset in the form of a Data.Frame. SpaMTP provides 4 pre-cleaned and optimised databases (`HMDB_db`, `Lipidmaps_db`, `Chebi_db`, `GNPS_db`).
#' @param ppm_error Numeric value indicating the size of the ppm error allowed when matching molecular weights between Seurat object and reference dataset. If only want exact matches set ppm = 0 (default = NULL).
#' @param adducts List of adducts to use for searching the database (e.g. "M+NH4","M+Na","M+CH3OH+H","M+K" etc.). For all possible adducts please visit [here](https://github.com/GenomicsMachineLearning/SpaMTP/blob/main/R/MZAnnotation.R#L305). If NULL will take the full list of adducts (default = NULL).
#' @param polarity Character string defining the polarity of adducts to use, either "positive", "negative" or "neutral" (default = "positive").
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion `[ion mass,m/z]`/`[Full width at half height]`. This value is used to estimate ppm_error when set to NULL (default = 30000).
#' @param verbose Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).
#'
#' @returns A data.frame containing all successfully annotated m/z values, with their corresponding annotation.
#' @export
#'
#' @examples
#' #cardinal <- readImzML("./Test_Data/Spotted/test_data1")
#' #mzs <- data.frame(Cardinal::featureData(cardinal))$mz
#' #results <- AnnotateBigData(mzs, db = HMDB_db, ppm_error = 3, adducts = c("M-H", "M+Cl"), polarity = "negative")
#' #cardinal_subset <- Cardinal::subset(cardinal, mz %in% results$observed_mz)
#' #SpaMTP_data <- CardinalToSeurat(cardinal_subset)
AnnotateBigData <- function(mzs, db,  ppm_error = NULL, adducts = NULL,polarity = "positive", tof_resolution = 30000,verbose = TRUE){
  mz_df <- data.frame(mz = mzs)
  mz_df$row_id <- seq(1, length(mz_df[["mz"]]))
  mz_df <- mz_df[c("row_id", "mz")]

  annotations <- annotateTable(mz_df, db = db, ppm_error = ppm_error, adducts = adducts, polarity = polarity,tof_resolution = tof_resolution,verbose = verbose)

  result_df <- annotations %>%
    dplyr::group_by(observed_mz) %>%
    dplyr::summarise(
      all_IsomerNames = paste(IsomerNames, collapse = "; "),
      all_Isomers = paste(Isomers, collapse = "; "),
      all_Isomers_IDs = paste(Isomers_IDs, collapse = "; "),
      all_Adducts = paste(unique(Adduct), collapse = "; "),
      all_Formulas = paste(unique(Formula), collapse = "; "),
      all_Errors = paste(round(Error,4), collapse = "; ")
    )

  rownames(result_df) <- paste0("mz-",result_df$observed_mz)
  result_df$mz_names <- rownames(result_df)
  result_df <- result_df %>% dplyr::mutate(present = TRUE)
  result_df <- data.frame(result_df)
  rownames(result_df) <- 1:length(result_df$observed_mz)
  return(result_df)
}
