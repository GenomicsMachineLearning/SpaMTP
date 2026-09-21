# SpaMTP: Spatial Metabolite Annotation

## Overview

`SpaMTP` provides analysis and visualisation methods for spatial
metabolomics and paired spatial multi-omics experiments. It uses
standard Seurat, Cardinal, and Bioconductor infrastructure. This
standalone release bundles its annotation databases; example experiments
are available from the tutorial download links. It does not require
companion resource packages.

This vignette demonstrates a small, fully reproducible accurate-mass
annotation workflow. It deliberately uses an in-memory database so that
package checking does not require a network connection or a large data
download.

``` r

library(SpaMTP)
```

## Inspect available database resources

[`SpaMTPDatabaseInfo()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/SpaMTPDatabaseInfo.md)
reads the lightweight bundled registry without loading the large
annotation tables or accessing the network.

``` r

database_registry <- SpaMTPDatabaseInfo(version = "latest")
database_registry[
  database_registry$resource %in% c("chem_props", "smiles_features"),
  c("resource", "version", "category", "rows", "columns")
]
##           resource version  category   rows columns
## 1       chem_props   3.0.7      core 289754      11
## 16 smiles_features   3.0.7 structure 284818      36
```

## Create a small structure-aware database

The example contains D-2-hydroxyglutaric acid and glucose. SMILES
connectivity is decomposed into interpretable functional groups and
ionisation sites before candidate generation.

``` r

demo_database <- data.frame(
  ramp_id = c("RAMP_C_2HG", "RAMP_C_GLUCOSE"),
  chem_source_id = c("HMDB:HMDB0000606", "HMDB:HMDB0000122"),
  common_name = c("D-2-hydroxyglutaric acid", "D-glucose"),
  monoisotop_mass = c(148.037173, 180.063388),
  mol_formula = c("C5H8O5", "C6H12O6"),
  iso_smiles = c("O=C(O)CC(O)C(=O)O", "OCC1OC(O)C(O)C(O)C1O"),
  stringsAsFactors = FALSE
)

demo_database <- AnnotateSMILESStructure(demo_database)
demo_database[, c(
  "common_name", "carboxyl_sites", "hydroxyl_sites",
  "max_exchangeable_protons", "positive_mode_score",
  "negative_mode_score", "alkali_affinity_score"
)]
##                common_name carboxyl_sites hydroxyl_sites
## 1 D-2-hydroxyglutaric acid              2              3
## 2                D-glucose              0              5
##   max_exchangeable_protons positive_mode_score negative_mode_score
## 1                        2                0.52                0.94
## 2                        0                0.47                0.45
##   alkali_affinity_score
## 1                  0.99
## 2                  0.87
```

## Build and query an indexed search space

Candidate masses are generated once and sorted. Repeated observed peaks
can then be queried using a binary range search rather than scanning the
full database for every peak.

``` r

annotation_index <- BuildMZAnnotationIndex(
  db = demo_database,
  polarity = "negative",
  adducts = c("[M-H]-", "[2M-H]-"),
  infer_structure = "always"
)

observed_mz <- 147.02935
candidates <- QueryMZAnnotationIndex(
  observed_mz = observed_mz,
  index = annotation_index,
  ppm = 5
)

candidates[, c(
  "observed_mz", "expected_mz", "metabolite_names", "adduct",
  "ppm_error", "score"
)]
##   observed_mz expected_mz         metabolite_names adduct ppm_error      score
## 1    147.0293    147.0299 D-2-hydroxyglutaric acid    M-H  3.717158 0.07816111
```

## Reproducibility

For production analyses, record the database version returned by
[`SpaMTPDatabaseInfo()`](https://genomicsmachinelearning.github.io/SpaMTP/developmental/reference/SpaMTPDatabaseInfo.md)
together with ion mode, matrix profile, adduct rules, mass tolerance,
and annotation-score threshold.

``` r

sessionInfo()
## R version 4.6.1 (2026-06-24)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.5 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
## 
## time zone: UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] SpaMTP_0.99.0    BiocStyle_2.40.0
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.5        
##   [4] spatstat.utils_3.2-5   farver_2.1.2           rmarkdown_2.32        
##   [7] fs_2.1.0               ragg_1.5.2             vctrs_0.7.3           
##  [10] ROCR_1.0-12            spatstat.explore_3.8-2 htmltools_0.5.9       
##  [13] sass_0.4.10            sctransform_0.4.3      parallelly_1.48.0     
##  [16] KernSmooth_2.23-26     bslib_0.12.0           htmlwidgets_1.6.4     
##  [19] desc_1.4.3             ica_1.0-3              plyr_1.8.9            
##  [22] plotly_4.12.1          zoo_1.9-0              cachem_1.1.0          
##  [25] igraph_2.3.3           mime_0.13              lifecycle_1.0.5       
##  [28] pkgconfig_2.0.3        Matrix_1.7-5           R6_2.6.1              
##  [31] fastmap_1.2.0          fitdistrplus_1.2-6     future_1.75.0         
##  [34] shiny_1.14.0           digest_0.6.39          patchwork_1.3.2       
##  [37] S4Vectors_0.50.3       Seurat_5.5.1           tensor_1.5.1          
##  [40] RSpectra_0.16-2        irlba_2.3.7            textshaping_1.0.5     
##  [43] progressr_1.0.0        spatstat.sparse_3.2-0  httr_1.4.9            
##  [46] polyclip_1.10-7        abind_1.4-8            compiler_4.6.1        
##  [49] proxy_0.4-29           S7_0.2.2               BiocParallel_1.46.0   
##  [52] DBI_1.3.0              fastDummies_1.7.6      MASS_7.3-65           
##  [55] classInt_0.4-11        units_1.0-1            tools_4.6.1           
##  [58] lmtest_0.9-40          otel_0.2.0             httpuv_1.6.17         
##  [61] future.apply_1.20.2    goftest_1.2-3          glue_1.8.1            
##  [64] nlme_3.1-169           promises_1.5.0         sf_1.1-3              
##  [67] grid_4.6.1             Rtsne_0.17             cluster_2.1.8.2       
##  [70] reshape2_1.4.5         generics_0.1.4         gtable_0.3.6          
##  [73] spatstat.data_3.1-9    class_7.3-23           tidyr_1.3.2           
##  [76] data.table_1.18.6.1    sp_2.2-3               BiocGenerics_0.58.1   
##  [79] spatstat.geom_3.8-3    RcppAnnoy_0.0.23       ggrepel_0.9.8         
##  [82] RANN_2.6.3             pillar_1.11.1          stringr_1.6.0         
##  [85] spam_2.11-4            RcppHNSW_0.7.0         limma_3.68.5          
##  [88] later_1.4.8            splines_4.6.1          dplyr_1.2.1           
##  [91] lattice_0.22-9         survival_3.8-6         deldir_2.0-4          
##  [94] tidyselect_1.2.1       CardinalIO_1.10.0      miniUI_0.1.2          
##  [97] pbapply_1.7-5          knitr_1.52             gridExtra_2.3.1       
## [100] bookdown_0.48          ProtGenerics_1.44.0    matter_2.14.0         
## [103] scattermore_1.2        stats4_4.6.1           xfun_0.61             
## [106] Biobase_2.72.0         statmod_1.5.2          matrixStats_1.5.0     
## [109] stringi_1.8.9          yaml_2.3.12            evaluate_1.0.5        
## [112] codetools_0.2-20       tibble_3.3.1           BiocManager_1.30.27   
## [115] cli_3.6.6              ontologyIndex_2.12     uwot_0.2.5            
## [118] xtable_1.8-8           reticulate_1.47.0      systemfonts_1.3.2     
## [121] jquerylib_0.1.4        Rcpp_1.1.2             globals_0.19.1        
## [124] spatstat.random_3.5-1  zeallot_0.2.0          png_0.1-9             
## [127] spatstat.univar_3.2-0  parallel_4.6.1         pkgdown_2.2.1         
## [130] ggplot2_4.0.3          dotCall64_1.2          listenv_1.0.0         
## [133] viridisLite_0.4.3      e1071_1.7-17           scales_1.4.0          
## [136] ggridges_0.5.7         SeuratObject_5.4.0     Cardinal_3.14.0       
## [139] purrr_1.2.2            rlang_1.3.0            cowplot_1.2.0         
## [142] shinyjs_2.1.1
```
