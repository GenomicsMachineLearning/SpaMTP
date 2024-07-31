#### THIS FILE ONLY CONTAINS DOCUMENTATION FOR THE ATTACHED CLEANED METABOLITE DATASETS ####


#' @title HMDB_db: A cleaned version of the reference metabolomics dataset from the Human Metabolome Database (HMDB)
#'
#' @description A dataset containing metabolite information from the HMDB,
#' including chemical formulas, exact masses, isomers, InChIKeys, and names.
#'
#' @format ## A data frame with 26190 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source <https://hmdb.ca/>
#'
#' @references Wishart DS, et al. (2018). HMDB 4.0: the human metabolome database for 2018. Nucleic Acids Res. 46(D1):D608-D617.
"HMDB_db"



#' @title Lipidmaps_db: A cleaned version of the lipid database from LIPID MAPS
#'
#' @description This object contains a collection of lipids from the LIPID MAPS Structure Database.
#'
#' @format ## A data frame with 9493 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source <https://www.lipidmaps.org/>
"Lipidmaps_db"




#' @title Chebi_db: Cleaned ChEBI `(Chemical entities of biological interest)` reference dataset
#'
#' @description A dataset containing information on chemical entities, including their formulas, exact masses, isomers, InChIKeys, and names.
#'
#' @format ## A data frame with 46297 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source <https://www.ebi.ac.uk/chebi/>
"Chebi_db"


#' @title GNPS_db: A cleaned database of metabolites from GNPS
#'
#' @description This object contains a collection of metabolites from the Global Natural Products Social Molecular Networking (GNPS) platform.
#'
#' @format ## A data frame with 489 rows and 37 variables:
#' \describe{
#'   \item{formula}{Chemical formula of the entity (character)}
#'   \item{exactmass}{Exact mass of the entity (numeric)}
#'   \item{isomers}{Presence of isomers (character)}
#'   \item{isomers_inchikey}{InChIKeys for isomers (character)}
#'   \item{isomers_names}{Names of isomers (character)}
#'   \item{2M-H}{Mass of the entity as a 2M-H ion (double)}
#'   \item{2M+ACN+H}{Mass of the entity as a 2M+ACN+H ion (double)}
#'   \item{2M+ACN+Na}{Mass of the entity as a 2M+ACN+Na ion (double)}
#'   \item{2M+FA-H}{Mass of the entity as a 2M+FA-H ion (double)}
#'   \item{2M+H}{Mass of the entity as a 2M+H ion (double)}
#'   \item{2M+Hac-H}{Mass of the entity as a 2M+Hac-H ion (double)}
#'   \item{2M+K}{Mass of the entity as a 2M+K ion (double)}
#'   \item{2M+NH4}{Mass of the entity as a 2M+NH4 ion (double)}
#'   \item{2M+Na}{Mass of the entity as a 2M+Na ion (double)}
#'   \item{3M-H}{Mass of the entity as a 3M-H ion (double)}
#'   \item{M+2ACN+H}{Mass of the entity as a M+2ACN+H ion (double)}
#'   \item{M+2K+H}{Mass of the entity as a M+2K+H ion (double)}
#'   \item{M+2Na-H}{Mass of the entity as a M+2Na-H ion (double)}
#'   \item{M+ACN+H}{Mass of the entity as a M+ACN+H ion (double)}
#'   \item{M+ACN+Na}{Mass of the entity as a M+ACN+Na ion (double)}
#'   \item{M+Br}{Mass of the entity as a M+Br ion (double)}
#'   \item{M+CH3OH+H}{Mass of the entity as a M+CH3OH+H ion (double)}
#'   \item{M+Cl}{Mass of the entity as a M+Cl ion (double)}
#'   \item{M+DMSO+H}{Mass of the entity as a M+DMSO+H ion (double)}
#'   \item{M+FA-H}{Mass of the entity as a M+FA-H ion (double)}
#'   \item{M+H}{Mass of the entity as a M+H ion (double)}
#'   \item{M+Hac-H}{Mass of the entity as a M+Hac-H ion (double)}
#'   \item{M+IsoProp+H}{Mass of the entity as a M+IsoProp+H ion (double)}
#'   \item{M+IsoProp+Na+H}{Mass of the entity as a M+IsoProp+Na+H ion (double)}
#'   \item{M+K}{Mass of the entity as a M+K ion (double)}
#'   \item{M+K-2H}{Mass of the entity as a M+K-2H ion (double)}
#'   \item{M+NH4}{Mass of the entity as a M+NH4 ion (double)}
#'   \item{M+Na}{Mass of the entity as a M+Na ion (double)}
#'   \item{M+Na-2H}{Mass of the entity as a M+Na-2H ion (double)}
#'   \item{M+TFA-H}{Mass of the entity as a M+TFA-H ion (double)}
#'   \item{M-H}{Mass of the entity as a M-H ion (double)}
#'   \item{M-H2O-H}{Mass of the entity as a M-H2O-H ion (double)}
#'   ...
#' }
#'
#' @source \url{https://gnps.ucsd.edu/: https://gnps.ucsd.edu/}
"GNPS_db"



#' @title adduct_file: A dataframe containing possible adducts used for pathway analysis
#'
#' @description This object contains a collection of adducts with their relative ion.mass, change and polarity
#'
#' @format ## A data frame with 47 rows and 6 variables:
#' \describe{
#'   \item{adduct_name}{Chemical formula of the adduct (character)}
#'   \item{ion.mass}{Mass of the molecular ion formula (character)}
#'   \item{charge}{Relative charge of the adduct (integer)}
#'   \item{mult}{Multiplication factor based on the original molecule mass (double)}
#'   \item{add_mass}{The true mass of the ion being added or reduced (double)}
#'   \item{pol}{Polarity of the ion (reduced = negative, added = positive) (character)}
#' }
#'
"adduct_file"


#' @title analyte: A dataframe containing ID's of possible RAMP analytes
#'
#' @description This object contains a collection of Ramp_DB analytes, along with an annotation of their analyte type
#'
#' @format ## A data frame with 250,961 rows and 2 variables:
#' \describe{
#'   \item{rampId}{Ramp_DB analyte ID/Name (character)}
#'   \item{type}{Analyte type (character)}
#' }
#'
"analyte"


#' @title analytehaspathway: A dataframe containing RAMP_pathway ID's
#'
#' @description This object contains a collection of RAMP_DB ID's, Pathway ID's and relative pathway database source
#'
#' @format ## A data frame with 819,103 rows and 3 variables:
#' \describe{
#'   \item{rampId}{Ramp_DB analyte ID/Name (character)}
#'   \item{pathwayRampId}{Ramp_DB pathway ID (character)}
#'   \item{pathwaySource}{Relative source database for the respective pathway (character)}
#' }
#'
"analytehaspathway"


#' @title chem_props: A database containing the chemical properties and metadata of each RAMP_DB analyte
#'
#' @description This object contains a collection of RAMP_DB analytes with their corresponding metadata including chemical structure key (smiles), isotop mass, common name and molecular fomular
#'
#' @format ## A data frame with 283,003 rows and 11 variables:
#' \describe{
#'   \item{rampId}{Ramp_DB analyte ID/Name (character)}
#'   \item{chem_data_source}{Relative source database for the respective analyte (character)}
#'   \item{chem_source_id}{Relative source database ID for analyste (character)}
#'   \item{iso_smiles}{Smile structure for relative analyte (character)}
#'   \item{inchi_key_prefix}{InChlKey prefix of Internation Chemical Identifier (InChl) (character)}
#'   \item{inchi_key}{Full InChlKey for corresponding analyte (character)}
#'   \item{inchi}{Full InChl identifier structure for analyte (character)}
#'   \item{mw}{Molecular weight for analyte (double)}
#'   \item{monoisotop_mass}{Relative monoisotopic mass for analyate (double)}
#'   \item{common_name}{Analytes common name (character)}
#'   \item{mol_formula}{Analytes simplified molecular fomula (character)}
#' }
#'
"chem_props"


#' @title pathway: A dataframe containing RAMP_DB pathways and their relative metadata
#'
#' @description This object contains a collection of RAMP_DB pathways, their source ID, pathway catagory and common name
#'
#' @format ## A data frame with 53,952 rows and 5 variables:
#' \describe{
#'   \item{pathwayRampId}{Relative source database ID for analyste (character)}
#'   \item{rampId}{RAMP_DB pathway ID (character)}
#'   \item{IDtype}{Relative source database for the respective pathway (character)}
#'   \item{geneOrCompound}{Catagory grouping of respecitive pathway (character)}
#'   \item{pathwayName}{Common name of pathway (character)}
#' }
#'
"pathway"



#' @title source_df: A dataframe containing source information about RAMP_ID analyte used for analysis
#'
#' @description This object contains a collection of RAMP_DB analytes, their source ID, common name and pathway count
#'
#' @format ## A data frame with 789,665 rows and 8 variables:
#' \describe{
#'   \item{sourceId}{Ramp_DB analyte ID(character)}
#'   \item{rampId}{Analyte source ID (character)}
#'   \item{type}{Analyte type (character)}
#'   \item{commonName}{Common name of analyte (character)}
#'   \item{priorityHMDBStatus}{Priority level of the analyte (character)}
#'   \item{dataSource}{Relative source databases where analyte is represented (character)}
#'   \item{pathwayCount}{Number of pathways analyte is present in (integer)}
#' }
#'
"source_df"
