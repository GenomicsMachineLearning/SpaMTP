% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MZAnnotation.R
\name{AnnotateBigData}
\alias{AnnotateBigData}
\title{Annotates vector of m/z values}
\usage{
AnnotateBigData(
  mzs,
  db,
  ppm_error = NULL,
  adducts = NULL,
  polarity = "positive",
  tof_resolution = 30000,
  verbose = TRUE
)
}
\arguments{
\item{mzs}{Vector containing m/z values for annotation.}

\item{db}{Reference metabolite dataset in the form of a Data.Frame. SpaMTP provides 4 pre-cleaned and optimised databases (\code{HMDB_db}, \code{Lipidmaps_db}, \code{Chebi_db}, \code{GNPS_db}).}

\item{ppm_error}{Numeric value indicating the size of the ppm error allowed when matching molecular weights between Seurat object and reference dataset. If only want exact matches set ppm = 0 (default = NULL).}

\item{adducts}{List of adducts to use for searching the database (e.g. "M+NH4","M+Na","M+CH3OH+H","M+K" etc.). For all possible adducts please visit \href{https://github.com/GenomicsMachineLearning/SpaMTP/blob/main/R/MZAnnotation.R#L305}{here}. If NULL will take the full list of adducts (default = NULL).}

\item{polarity}{Character string defining the polarity of adducts to use, either "positive", "negative" or "neutral" (default = "positive").}

\item{tof_resolution}{is the tof resolution of the instrument used for MALDI run, calculated by ion \verb{[ion mass,m/z]}/\verb{[Full width at half height]}. This value is used to estimate ppm_error when set to NULL (default = 30000).}

\item{verbose}{Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).}
}
\value{
A data.frame containing all successfully annotated m/z values, with their corresponding annotation.
}
\description{
This function assigns each valid m/z peak with one/multiple metabolite names based on the mass difference between the observed value and the theoretical value documented in the reference database.
This function is to be used when dealing with large datasets as a preprocessing step. Users can annotate m/z values first and then subset their data accordinly before loading it into a SpaMTP Seurat Object.
}
\examples{
#cardinal <- readImzML("./Test_Data/Spotted/test_data1")
#mzs <- data.frame(Cardinal::featureData(cardinal))$mz
#results <- AnnotateBigData(mzs, db = HMDB_db, ppm_error = 3, adducts = c("M-H", "M+Cl"), polarity = "negative")
#cardinal_subset <- Cardinal::subset(cardinal, mz \%in\% results$observed_mz)
#SpaMTP_data <- CardinalToSeurat(cardinal_subset)
}
