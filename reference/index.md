# Package index

## Loading Spatial Metabolic Data into a SpaMTP Seurat Object

Functions that allow the user to load in SM data in different formats

- [`LoadSM()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/LoadSM.md)
  : Loads spatial metabolic data into a SpaMTP Seurat Object
- [`ReadSM_mtx()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ReadSM_mtx.md)
  : Read Spatial Metabolomics matrix file (.csv format)
- [`Load_METASPACE()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/Load_METASPACE.md)
  : Load METASPACE Data and Create Seurat Object

## Converting Between Data Objects

Functions that allow the user to convert between SpaMTP Seurat Objects
and Cardinal Objects

- [`CardinalToSeurat()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/CardinalToSeurat.md)
  : Converts a Cardinal Object into a SpaMTP Seurat Object
- [`ConvertSeuratToCardinal()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ConvertSeuratToCardinal.md)
  : Converts a SpaMTP Seurat object to a Cardinal object, including
  annotations and metadata

## Binning Spatial Metabolomic Data

Functions that bin m/z values into a lower resolution/wider peak.

- [`BinSpaMTP()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/BinSpaMTP.md)
  : Bin SpaMTP Object

## Annotating m/z Masses

Functions required for performing and handling m/z annotation using a
reference metabolic database

- [`AnnotateSM()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AnnotateSM.md)
  : Annotates m/z values stored in a SpaMTP Object
- [`AnnotateMZ()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AnnotateMZ.md)
  : Annotate one or more observed m/z values
- [`AdductRules()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AdductRules.md)
  : Return the validated SpaMTP adduct rule table
- [`MALDIMatrixProfiles()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MALDIMatrixProfiles.md)
  : List MALDI matrix and on-tissue derivatization profiles
- [`MALDIMatrixRules()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MALDIMatrixRules.md)
  : Select annotation rules from the MALDI matrix/reagent profile
- [`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/BuildMZAnnotationIndex.md)
  : Build a reusable indexed metabolite annotation search space
- [`QueryMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/QueryMZAnnotationIndex.md)
  : Query an indexed metabolite annotation search space
- [`AnnotationInfo()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AnnotationInfo.md)
  : Inspect the metabolite annotation provenance stored in a SpaMTP
  object
- [`print(`*`<spamtp_mz_index>`*`)`](https://genomicsmachinelearning.github.io/SpaMTP/reference/print.spamtp_mz_index.md)
  : Print an indexed metabolite annotation search space
- [`AddCustomMZAnnotations()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AddCustomMZAnnotations.md)
  : Assign custom annotations to m/z values
- [`AddFMP10Annotations()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AddFMP10Annotations.md)
  : Annotates FMP10 matrix data
- [`SearchAnnotations()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/SearchAnnotations.md)
  : Find Annotation
- [`GetMZMetadata()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/GetMZMetadata.md)
  : Gets values from a single metadata column for a respective m/z
  value.
- [`FindDuplicateAnnotations()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindDuplicateAnnotations.md)
  : Finds if any metabolite is duplicated across multiple m/z values.
- [`SubsetMZFeatures()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/SubsetMZFeatures.md)
  : Subset a SpaMTP Seurat Spatial Metabolomic object by a list of m/z's
- [`AnnotateBigData()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AnnotateBigData.md)
  : Annotates vector of m/z values
- [`getRefinedAnnotations()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/getRefinedAnnotations.md)
  : Refines and reduces m/z annotations

## Simplify Annotations

Functions required for simplifying m/z annotation results

- [`CalculateAnnotationStatistics()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/CalculateAnnotationStatistics.md)
  : Calculate annotation statistics for all m/z value suggesting the
  most likely metabolite based on correlated pathway expression.
- [`CalculateSingleAnnotationStatistics()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/CalculateSingleAnnotationStatistics.md)
  : Calculate annotation statistics for a single m/z value suggesting
  the most likely metabolite based on correlated pathway expression.
- [`Pseudo_msms()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/Pseudo_msms.md)
  : Find all library spectra10 above a cosine threshold for each imaging
  peak, restricting to a global precursor-mz range.
- [`CreatePathwayObject()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/CreatePathwayObject.md)
  : Create a SpaMTP Seurat Object containing expression values for all
  present pathways

## Simplifying Lipid Nomenclature

Function used to simplify lipid names into general lipid categories,
classes, and more

- [`RefineLipids()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RefineLipids.md)
  : Uses common lipid nomenclature to simplify lipid annotations

## Analysis of Differentially Expressed Peaks

Functions required for performing pseudo-bulking differential expression
analysis

- [`FindAllDEMs()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindAllDEMs.md)
  : Finds differentially expressed m/z values/metabolites between all
  comparison groups.

## Visualising DEPs Analysis

Functions used to generate a heatmap from DEP results

- [`DEMsHeatmap()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/DEMsHeatmap.md)
  : Heatmap of Differentially Expressed Metabolites

## Metabolic and Transcriptomic Pathway Analysis

Functions used to perform pathway analysis, both PCA and metabolite/gene
set-based (GSEA)

- [`FishersPathwayAnalysis()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FishersPathwayAnalysis.md)
  : Calculates Significant Metabolic Pathways using a Fisher Exact Test
- [`FindRegionalPathways()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindRegionalPathways.md)
  : Regional Pathway Enrichment
- [`RunRAMPgeseca()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RunRAMPgeseca.md)
  : Runs multilevel Monte-Carlo variant for performing gene sets
  co-regulation analysis using the RAMP_DB metabolite/gene database.
- [`CreatePathwayAssay()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/CreatePathwayAssay.md)
  : Create a Pathway Assay from Gene or Metabolite Data

## Metabolic and Transcriptomic Pathway Visualisation

Functions used to visualise pathway analysis results

- [`VisualisePathways()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/VisualisePathways.md)
  : Visualise Significant Pathways
- [`PlotRegionalPathways()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/PlotRegionalPathways.md)
  : Plot significantly enriched pathways per region
- [`PlotPathways()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/PlotPathways.md)
  : Plots the expression profile of a feature set corresponding to
  specified pathways onto a 2D scatter plot based on a dimensionality
  reduction technique.
- [`PlotPathwaysSpatially()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/PlotPathwaysSpatially.md)
  : Plots the spatial expression profile of a feature set corresponding
  to specified pathways

## Pathway Network Plotting

Functions used to generate network plots for specified patheays

- [`PathwayNetworkPlots()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/PathwayNetworkPlots.md)
  : Construct an interactive pathway network visualization

## Dimentionality Reduction Analysis

Functions that are used for calculating PCA embeddings and projections
based on SM and/or ST data

- [`RunMetabolicPCA()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RunMetabolicPCA.md)
  : Generates PCA analysis results for a SpaMTP Seurat Object
- [`RunSpatialGraphPCA()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RunSpatialGraphPCA.md)
  : Perform Dimensionality Reduction using Graph-Regularised PCA on
  Spatial Data
- [`GetKmeanClusters()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/GetKmeanClusters.md)
  : Perform K-means clustering on a specified reduction

## SpaMTP Metabolic Data Visualisation

Functions that can be used to visualize data from SpaMTP object

- [`ImageMZPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ImageMZPlot.md)
  : Plot expression of m/z values spatially
- [`ImageMZAnnotationPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ImageMZAnnotationPlot.md)
  : Plot expression of annotated metabolites spatially
- [`SpatialMZPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/SpatialMZPlot.md)
  : Plot expression of m/z values spatially for a Spatial SpaMTP Seurat
  Objects.
- [`SpatialMZAnnotationPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/SpatialMZAnnotationPlot.md)
  : Plot expression of metabolites in spatially from a Spatial SpaMTP
  Objects.
- [`Plot3DFeature()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/Plot3DFeature.md)
  : Generates a 3D spatial feature plot from a SpaMTP object
- [`MassIntensityPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MassIntensityPlot.md)
  : Plot mass intensity spectra
- [`DensityMap()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/DensityMap.md)
  : Generates interactive 3D spatial density plot for m/z values
- [`CheckAlignment()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/CheckAlignment.md)
  : Check multi-modal coordinate alignment

## Interactive Spatial Binning Visualisation

Interactive plot that displays spatial changes to m/z intensity values
based on changes to bin size.

- [`InteractiveSpatialPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/InteractiveSpatialPlot.md)
  : Interactive Spatial Plot for visualising different m/z bin sizes

## Additional SpaMTP Functions

Functions that can be used to find the closest metabolite and bin the
expression of multiple metabolites into one

- [`FindNearestMZ()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindNearestMZ.md)
  : Finds the nearest m/z peak to a given value in a SpaMTP Object
- [`BinMetabolites()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/BinMetabolites.md)
  : Sums the intensity values of multiple m/z values into one

## Spatial Analysis of Metabolomic Data

Functions used to identify spatially correlated features
(metabolites/genes)

- [`FindCorrelatedFeatures()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindCorrelatedFeatures.md)
  : Find top features and metabolites that are strongly correlated with
  a given feature
- [`FindSpatiallyVariableMetabolites()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/FindSpatiallyVariableMetabolites.md)
  : Find Spatially Variable Metabolites
- [`GetSpatiallyVariableMetabolites()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/GetSpatiallyVariableMetabolites.md)
  : Get top spatially variable metabolites
- [`RowVar()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RowVar.md)
  : Compute the row variances for each m/z value

## Multi-Omic Data Integration

Functions used to Align, Map and Integrate Spatial Metabolomic and
Transcriptomic data

- [`MapSpatialOmics()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MapSpatialOmics.md)
  : Maps Spatial Metabolomic (MALDI) data to corresponding Spatial
  Transcriptomics data and coordinates.
- [`ApplySpatialAlignment()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ApplySpatialAlignment.md)
  : Apply automatic or precomputed spatial alignment
- [`AlignSpatialOmics()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AlignSpatialOmics.md)
  : Interactive app for SM and ST coordinate alignment
- [`MultiOmicIntegration()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MultiOmicIntegration.md)
  : Mult-Omic data integration
- [`CreateMergedModalityAssay()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/CreateMergedModalityAssay.md)
  : Create a singular multiomics assay by merging data from multiple
  assays.

## Adding Image Data

Functions used to align and add an image to a spatial metabolic SpaMTP
object.

- [`AddSMImage()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AddSMImage.md)
  : Manually align an image (e.g. H&E, Immuno) to a SM SpaMTP dataset

## ROI Selection

Functions used subset or select specific spatial regions of interest.

- [`SelectROIs()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/SelectROIs.md)
  : Launch an Interactive ROI Annotation App for Seurat Spatial Data

## Pre-Processing SpaMTP Metabolic Data

Functions for normalising and visualising the pre-processing of SpaMTP
datasets

- [`NormalizeSMData()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/NormalizeSMData.md)
  : Normalizes m/z intensity data stored in a SpaMTP Seurat Object
- [`TMMNormalize()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/TMMNormalize.md)
  : Performs TMM normalization between categories based on a specified
  ident
- [`MZRidgePlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MZRidgePlot.md)
  : Generates a ridge plot of spatial metabolic intensity data
- [`MZVlnPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MZVlnPlot.md)
  : Generates a violin plot of spatial metabolic intensity data
- [`MZBoxPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/MZBoxPlot.md)
  : Generates a Boxplot of spatial metabolic intensity data

## Exporting SpaMTP Data

Function to export SpaMTP data in .mtx, barcodes.csv, features.csv,
metadata.csv, and feature.metadata.csv files

- [`SaveSpaMTPData()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/SaveSpaMTPData.md)
  : Saves SpaMTP Object

## Reference Metabolite Datasets

Metabolite datasets used for annotating m/z masses

- [`HMDB_db`](https://genomicsmachinelearning.github.io/SpaMTP/reference/HMDB_db.md)
  : HMDB_db: A cleaned version of the reference metabolomics dataset
  from the Human Metabolome Database (HMDB)

- [`Lipidmaps_db`](https://genomicsmachinelearning.github.io/SpaMTP/reference/Lipidmaps_db.md)
  : Lipidmaps_db: A cleaned version of the lipid database from LIPID
  MAPS

- [`Chebi_db`](https://genomicsmachinelearning.github.io/SpaMTP/reference/Chebi_db.md)
  :

  Chebi_db: Cleaned ChEBI `(Chemical entities of biological interest)`
  reference dataset

- [`GNPS_db`](https://genomicsmachinelearning.github.io/SpaMTP/reference/GNPS_db.md)
  : GNPS_db: A cleaned database of metabolites from GNPS

- [`filtered_fmp10`](https://genomicsmachinelearning.github.io/SpaMTP/reference/filtered_fmp10.md)
  : filtered_fmp10: data.frame containing FMP10+ metabolite mappings

## Metabolic Pathway Datasets

Various datasets required for pathway analysis

- [`adduct_file`](https://genomicsmachinelearning.github.io/SpaMTP/reference/adduct_file.md)
  : adduct_file: A dataframe containing possible adducts used for
  pathway analysis
- [`analyte`](https://genomicsmachinelearning.github.io/SpaMTP/reference/analyte.md)
  : analyte: A dataframe containing ID's of possible RAMP analytes
- [`analytehaspathway`](https://genomicsmachinelearning.github.io/SpaMTP/reference/analytehaspathway.md)
  : analytehaspathway: A dataframe containing RAMP_pathway ID's
- [`chem_props`](https://genomicsmachinelearning.github.io/SpaMTP/reference/chem_props.md)
  : chem_props: A database containing the chemical properties and
  metadata of each RAMP_DB analyte
- [`pathway`](https://genomicsmachinelearning.github.io/SpaMTP/reference/pathway.md)
  : pathway: A dataframe containing RAMP_DB pathways and their relative
  metadata
- [`source_df`](https://genomicsmachinelearning.github.io/SpaMTP/reference/source_df.md)
  : source_df: A dataframe containing source information about RAMP_ID
  analyte used for analysis
- [`ramp_db_metadata`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ramp_db_metadata.md)
  : Metadata for the bundled pruned RaMP snapshot

## Pathway Network Datasets

Various datasets required for generated pathway network plots

- [`RAMP_hmdb`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RAMP_hmdb.md)
  : RAMP_hmdb: A list containing network plot information about pathways
  from the HMDB database
- [`RAMP_Reactome`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RAMP_Reactome.md)
  : RAMP_Reactome: A list containing network plot information about
  pathways from the Reactome database
- [`RAMP_kegg`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RAMP_kegg.md)
  : RAMP_kegg: A list containing network plot information about pathways
  from the KEGG database
- [`RAMP_wikipathway`](https://genomicsmachinelearning.github.io/SpaMTP/reference/RAMP_wikipathway.md)
  : RAMP_wikipathway: A list containing network plot information about
  pathways from the Wiki database
- [`reaction_type`](https://genomicsmachinelearning.github.io/SpaMTP/reference/reaction_type.md)
  : reaction_type: data.frame containing reaction type mappings

## Cardinal Wrapper Functions

Functions used to alter Cardinal Objects.

- [`add_ssc_annotation()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/add_ssc_annotation.md)
  : Adds Cardinal ssc segmentation annotation to m/z count data object

## Additional Worker Functions

Various helper functions used by SpaMTP

- [`verbose_message()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/verbose_message.md)
  : Helper function for suppressing function progress messages
- [`subset_SPM()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/subset_SPM.md)
  : Subsets SpaMTP Seurat Object containing FOVs
- [`check_cardinal_version()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/check_cardinal_version.md)
  : Check Cardinal Version

## Utility Functions

### Binning Helper Functions

Helper functions for binning m/z intensity spectra for SpaMTP objects

- [`BinnedCardinalToSeurat()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/BinnedCardinalToSeurat.md)
  : Converts a SpaMTP binned Cardinal Object into a SpaMTP Seurat Object
- [`spectral_binning()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/spectral_binning.md)
  : Spectral binning of intensity values stored in a Matrix object,
  converted from matter.

### Plotting Helper Functions

Helper functions used for SpaMTP data visualisation

- [`bin.mz()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/bin.mz.md)
  : Bins multiple m/z values into one.
- [`plusminus()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/plusminus.md)
  : Identifies all mz peaks within a plus-minus range of the target_mz
- [`plot_plus_minus()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/plot_plus_minus.md)
  : Helper Function to generate merged counts within the plus minus
  range provided
- [`check_column_type()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/check_column_type.md)
  : Helper function to determine if a column contains categorical or
  continuous Data
- [`pixelPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/pixelPlot.md)
  : Helper function for converting Seurat Class ggplots from spot to
  pixel layout

### Annotation Helper Functions

Helper functions required for running ‘AnnotateSM()’

- [`annotateTable()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/annotateTable.md)
  : Annotates m/z values stored in a data.frame based on reference
  metabolite dataset

- [`labels_to_show()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/labels_to_show.md)
  : Filters the annotation list to only include the first n number of
  annotations per m/z

- [`add_backslashes_to_specialfeatures()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/add_backslashes_to_specialfeatures.md)
  : Adds in backslashes required to take into account special using
  grepl such as brackets and +

- [`check_and_truncate_adduct_vector()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/check_and_truncate_adduct_vector.md)
  : Checks if the complete adduct is in the data base, else returns a
  truncated adduct

- [`db_adduct_filter()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/db_adduct_filter.md)
  : Filters the provided metabolomic database by polarity and adducts

- [`is_formula_valid()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/is_formula_valid.md)
  : Checks if a formula contains only the allowed elements

- [`formula_filter()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/formula_filter.md)
  : Filters reference Database to only select natural elements

- [`calculate_bounds()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/calculate_bounds.md)
  : Calculates the mz range of the observed_df

- [`ppm_error()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ppm_error.md)
  : Calculates the ppm error as a valve

- [`ppm_range_match()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/ppm_range_match.md)
  :

  Calculates the ppm range and check if mz values are within the range

  - Returns TRUE if match is found and false if no match.

- [`proc_db()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/proc_db.md)
  : Searches observed mz values against the data base list and returns
  matching annotations

### METASPACE Loader helper functions

Helper functions required for running ‘Load_METASPACE()’

- [`metaspace_client()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/metaspace_client.md)
  : Metaspace R Client Constructor
- [`get_metaspace()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/get_metaspace.md)
  : Core METASPACE Data Retrieval
- [`metaspace_to_feature_matrix()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/metaspace_to_feature_matrix.md)
  : Convert METASPACE Images to Feature Matrix

### Differential Abundance Helper Functions

Helper functions for calculating and plotting differentially expressed
metabolites

- [`run_pooling()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/run_pooling.md)
  : Pools SpaMTP Seurat object into random pools for pseudo-bulking.
- [`run_DE()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/run_DE.md)
  : Runs EdgeR analysis for pooled data
- [`save_pheatmap_as_pdf()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/save_pheatmap_as_pdf.md)
  : Saves a DEMsHeatmap as a PDF

### Pathway Analysis Helper Functions

Helper functions used for running pathway analysis

- [`PlotSinglePathway()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/PlotSinglePathway.md)
  : Plots the expression profile of a feature set corresponding to
  specified pathways onto a 2D scatter plot based on a dimensionality
  reduction technique.
- [`PlotSinglePathwaySpatially()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/PlotSinglePathwaySpatially.md)
  : Plot expression profile of a single RAMP pathway spatially
- [`addGesecaScores()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/addGesecaScores.md)
  : Add GESECA Scores to SpaMTP Object
- [`get_analytes_db()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/get_analytes_db.md)
  : Helper function for building a pathway db based on detected
- [`list_to_pprcomp()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/list_to_pprcomp.md)
  : Creates a pprcomp object based on an input list

### Multi-Omic Helper Functions

Helper functions used in multi-omic analysis functions

- [`kneighbors_graph()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/kneighbors_graph.md)
  : Construct a k-Nearest Neighbour Graph from Spatial Coordinates.
- [`get_square_coordinates()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/get_square_coordinates.md)
  : This function computes the coordinates of a square's four corners
  based on a given center point and width.
- [`lowresMapping()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/lowresMapping.md)
  : Maps SM pixels to low resolution ST data
- [`hiresMapping()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/hiresMapping.md)
  : Maps SM pixels to high resolution ST data
- [`translate()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/translate.md)
  : Creates a transformation matrix that translates an object in 2D
- [`mirror()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/mirror.md)
  : Creates a transformation matrix that mirrors an object in 2D along
  either the x axis or y axis around its center of mass
- [`stretch()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/stretch.md)
  : Stretch along angle
- [`rigid.rot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/rigid.rot.md)
  : Creates a transformation matrix for rotation
- [`rigid.transf()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/rigid.transf.md)
  : Creates a transformation matrix for rotation and translation
- [`rigid.transl()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/rigid.transl.md)
  : Creates a transformation matrix for translation with an offset of
  (h, k)
- [`rigid.refl()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/rigid.refl.md)
  : Creates a transformation matrix for reflection
- [`rigid.stretch()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/rigid.stretch.md)
  : Creates a transformation matrix for stretching
- [`combine.tr()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/combine.tr.md)
  : Combines rigid tranformation matrices

### Pre-Processing Helper Functions

Helper functions for generating pre-processing plots

- [`statPlot()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/statPlot.md)
  : Helper function for QC plots by generating intensity count data
