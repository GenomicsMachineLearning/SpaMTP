% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MZAnnotation.R
\name{proc_db}
\alias{proc_db}
\title{Searches observed mz values against the data base list and returns matching annotations}
\usage{
proc_db(observed_df, reference_df, ppm_threshold = 10)
}
\arguments{
\item{observed_df}{DataFrame containing the observed mz values}

\item{reference_df}{DataFrame contating the reference mz values and relative annotations.}

\item{ppm_threshold}{Numeric value defining the maximum acceptable ppm_error/threshold allowed between observed and reference mz values}
}
\value{
A DataFrame containing matched mz values between the observed and reference dataframes
}
\description{
Searches observed mz values against the data base list and returns matching annotations
}
\examples{
# HMDB_db <- load("data/HMDB_1_names.rds")
# mz_df <- SeuratObject[["Spatial"]][["mz"]]
# mz_df$row_id <- seq(1, length(mz_df$mz))

## 1) Filter DB by adduct.
# db_1 <- db_adduct_filter(HMDB_db, c("M+H"), polarity = "pos")

## 2) only select natural elements
# db_2 <- formula_filter(db_1)

## 3) search db against mz df return results
# db_3 <- proc_db(mz_df, db_2, ppm_threshold = 5)
}
