#' Common adduct constants
#'
#' A small table of legacy adduct definitions retained for compatibility with
#' earlier SpaMTP workflows. New annotation code should use [AdductRules()],
#' which includes charge, stoichiometry, ion-mode, and chemical-validity fields.
#'
#' @format A data frame with 47 rows and 6 variables:
#' \describe{
#'   \item{adduct_name}{Adduct notation.}
#'   \item{ion.mass}{Legacy ion-mass expression.}
#'   \item{charge}{Ion charge.}
#'   \item{mult}{Analyte stoichiometric multiplier.}
#'   \item{add_mass}{Exact mass shift.}
#'   \item{pol}{Ion polarity.}
#' }
#' @return A data frame of legacy adduct definitions.
#' @seealso [AdductRules()], [MALDIMatrixRules()]
#' @keywords datasets
"adduct_file"

#' Pathway-network reaction styles
#'
#' A small lookup table mapping RaMP reaction codes to visual properties used
#' by the interactive pathway-network viewer.
#'
#' @format A data frame with 11 rows and 5 variables:
#' \describe{
#'   \item{reaction_type}{Integer reaction code.}
#'   \item{reaction_name}{Human-readable reaction type.}
#'   \item{linetype}{Line style.}
#'   \item{arrowhead}{Arrowhead style.}
#'   \item{colour}{Edge colour.}
#' }
#' @return A data frame of pathway-network reaction styles.
#' @keywords datasets
"reaction_type"
