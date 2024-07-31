
#### SpaMTP Lipid Simplifier Functions #################################################################################################################################################################################

#' Uses common lipid nomenclature to simplify lipid annotations
#'
#' @param data data.frame containing a column with lipid annotations
#' @param annotation.column Character string matching the column name containing the lipid annotations in the provided df (default = "annotations").
#' @param database Character string defining the database matching the lipid annotations. Possible entries include c('Shorthand2020','Goslin','FattyAcids','LipidMaps','SwissLipids','HMDB') (default = "HMDB").
#' @param lipid_info Character string defining what level of information to return for each lipid. Options are either "all" or "simple". "Simple" returns a smaller, simplified list of annotations for each lipid. Please visit https://bioconductor.org/packages/release/bioc/vignettes/rgoslin/inst/doc/introduction.html to see all possible outputs (default = "simple")).
#'
#' @return Data.frame containing additional columns with simplified lipid names
#' @export
#'
#' @examples
#' # RefineLipids(DEPs_df)
RefineLipids <- function(data, annotation.column = "annotations", database = "HMDB", lipid_info = "simple"){


  if (!(lipid_info %in% c("simple", "all"))){
    warning("Invalid input for 'lipid_info'!... Must be either 'all' or 'simple'. Will use default == 'simple'")
    lipid_info <- "simple"
  }

  rows_list <- list()


  total_rows = length(rownames(data))
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                              max = total_rows, # Maximum value of the progress bar
                              style = 3,    # Progress bar style (also available style = 1 and style = 2)
                              width = 50,   # Progress bar width. Defaults to getOption("width")
                              char = "=")

  for (idx in 1:total_rows){
    row_name <- rownames(data)[idx]
    row <- data[row_name,][[annotation.column]]
    annotations <- strsplit(row, "; ")[[1]]

    suppressMessages({
      lipid.df <- rgoslin::parseLipidNames(annotations, grammar = database)
    })

    if (lipid_info == "all"){
      add_infomation <- c(colnames(lipid.df), "Species.Name.Simple")
    } else {
      add_infomation <- c("Lipid.Maps.Category", "Lipid.Maps.Main.Class", "Species.Name", "Species.Name.Simple")
    }

    pass <- unique(lipid.df$Grammar)

    if (length(pass) == 1 && pass == "NOT_PARSEABLE"){
      empty_df <- data.frame(matrix(NA, ncol = length(add_infomation), nrow = 1))
      colnames(empty_df) <- add_infomation
      rows_list[[row_name]] <- empty_df
    } else {
      lipid.df$Species.Name.Simple <- ifelse(
        is.na(lipid.df$Lipid.Maps.Main.Class) | is.na(lipid.df$Total.C) | is.na(lipid.df$Total.DB),
        NA,
        paste0(lipid.df$Lipid.Maps.Main.Class, "(", lipid.df$Total.C, ":", lipid.df$Total.DB, ")"))

      lipid.df  <- lipid.df[add_infomation]
      new.row.data <- data.frame(lapply(lipid.df, function(col) paste(unique(col), collapse = "; ")))
      rows_list[[row_name]] <- new.row.data
    }

    utils::setTxtProgressBar(pb, idx)
  }

  close(pb)


  final_df <- dplyr::bind_rows(rows_list)
  data[colnames(final_df)] <- final_df

  # For Returning columns with NA's
  #data <- data %>%
  #  dplyr::mutate(across(all_of(add_infomation), ~ ., .names = "All.{.col}"))

  data <- data %>%
    dplyr::mutate_at(vars(add_infomation), ~ stringr::str_replace_all(., "\\bNA\\b", ""))

  for (col in add_infomation) {

    data[[col]] <- gsub("(^|\\W)\\;\\s+|\\s+\\;($|\\W)", "\\1\\2", data[[col]])

    data[[col]] <- gsub("\\;\\s+$", "", data[[col]])

    if (col == "Species.Name.Simple"){
      data[[col]] <- gsub("[; ]", "", data[[col]])
      data[[col]] <- gsub("\\)", "); ", data[[col]])
      data[[col]] <-  sub(";\\s*$", "", data[[col]])
    }

  }

  data <-data %>%
    dplyr::mutate_at(vars(add_infomation), ~ stringr::str_trim(.))


  return(data)
}
