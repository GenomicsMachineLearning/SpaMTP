% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MZAnnotation.R
\name{annotateTable}
\alias{annotateTable}
\title{Annotates m/z values sotred in a data.frame based on reference metabolite dataset}
\usage{
annotateTable(
  mz_df,
  db,
  ppm_error = NULL,
  adducts = NULL,
  polarity = "positive",
  tof_resolution = 30000,
  verbose = TRUE
)
}
\arguments{
\item{mz_df}{dataframe containing m/z values for annotation.}

\item{db}{Reference metabolite dataset in the form of a Data.Frame.}

\item{ppm_error}{Numeric value indicating the size of the ppm error allowed when matching molecular weights between Seurat object and reference dataset. If only want exact matches set ppm = 0 (default = NULL).}

\item{adducts}{List of adducts to use for searching the database (e.g. "M+NH4","M+Na","M+CH3OH+H","M+K" etc.). For all possible adducts please visit \href{https://github.com/GenomicsMachineLearning/SpaMTP/blob/main/R/MZAnnotation.R#L305}{here}. If NULL will take the full list of adducts (default = NULL).}

\item{polarity}{Character string defining the polarity of adducts to use, either "positive", "negative" or "neutral" (default = "positive").}

\item{tof_resolution}{is the tof resolution of the instrument used for MALDI run, calculated by ion \verb{[ion mass,m/z]}/\verb{[Full width at half height]}. This value is used to estimate ppm_error when set to NULL (default = 30000).}

\item{verbose}{Boolean indicating whether to show the message. If TRUE the message will be show, else the message will be suppressed (default = TRUE).}
}
\value{
Generates an intermediate annotated m/z dataframe
}
\description{
Helper function for \code{AnnotatesSM()} and \code{FishersPathwayAnalysis()}.
}
\examples{

### HelperFunction
}
