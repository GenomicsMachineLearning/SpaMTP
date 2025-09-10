
#' Calculate annotation statistics for a single m/z value suggesting the most likely metabolite based on correlated pathway expression.
#'
#' This function evaluates the potential biological relevance of an annotation for a given m/z value by:
#' - Identifying pathways associated with each possible annotated metabolite
#' - Calculating colocalisation score between the m/z intensity and the expression of each corresponding pathway
#' - Ranking annotations by a combined z-score based on correlation strength and number of supporting significant pathways
#' - NOTE: this function requires `CreatePathwayAssay` and `CreatePathwayObject` to be run first
#'
#' @param mz Character or numeric. The target m/z feature. If numeric, the closest matching m/z in the dataset will be selected.
#' @param data A SpaMTP Seurat object containing both metabolite and pathway assays generated from `CreatePathwayObject`.
#' @param mz.assay Character string defining the name of the assay containing m/z features.
#' @param pathway.assay Character string matching the name of the assay containing pathway features (default = "pathway").
#' @param mz.slot Character string stating the slot to extract m/z values from (default = "scale.data").
#' @param pathway.slot Character string defining the slot to extract pathway features from (default = "scale.data").
#' @param corr_theshold Numeric value stating the correlation threshold to consider a pathway as significantly colocalized. If set to `0`, all pathways will be counted (default = 0).
#' @param corr_weight Numeric weight applied to correlation score in z-score calculation. If significance should be based more on the correlation, increase this value (default = 1).
#' @param n_weight Numeric weight applied to number of correlated pathways in z-score calculation (default = 1).
#'
#' @return A tibble with the ranked annotations for the m/z value, containing:
#' \describe{
#'   \item{metabolite}{Most common name of the metabolite associated with the RAMP ID}
#'   \item{ramp_id}{The RAMP ID corresponding to the annotation}
#'   \item{n_sig_path}{Number of correlated pathways above the threshold}
#'   \item{max_cor}{Maximum correlation value among significant pathways}
#'   \item{z_score}{Combined z-score used to rank annotations}
#'   \item{pval}{Unadjusted p-value}
#'   \item{pval_adj}{Adjusted p-value (BH method)}
#' }
#'
#' @importFrom Cardinal subset colocalized
#' @importFrom dplyr group_by summarise mutate slice_max pull
#' @importFrom tibble tibble
#' @importFrom SeuratObject CreateAssay5Object
#' @importFrom stats pnorm p.adjust
#' @importFrom utils head
#'
#' @examples
#' #data <- CreatePathwayObject(data,assay="SPT_pathway",slot = "scale.data")
#' #CalculateSingleAnnotationStatistics(mz = "mz-674.2805",data = data,mz.assay = "SPM",pathway.assay = "pathway",mz.slot = "scale.data")
#'
#' @export
CalculateSingleAnnotationStatistics <- function(mz, data, mz.assay, pathway.assay = "pathway", mz.slot= "scale.data", pathway.slot = "scale.data", corr_theshold = 0, corr_weight = 1, n_weight = 1){

  if(is.numeric(mz)){
    mz <- FindNearestMZ(data = data, target_mz = mz, assay = SM.assay)
  }

  DefaultAssay(data) <- pathway.assay
  pathway_list <- analytehaspathway %>%
    group_by(rampId) %>%
    summarise(pathways = list(unique(pathwayRampId)), .groups = "drop")

  named_list <- setNames(pathway_list$pathways, pathway_list$rampId)

  ## 2. Use source_df to match annotation sourceID to rampID

  row <- data[[mz.assay]]@meta.data[data[[mz.assay]]@meta.data$mz_names == mz,]
  source_df_copy <- source_df
  named_list_copy <- named_list

  message("Calculating top metabolites for provided m/z value. ")


  met_counts <- data[[mz.assay]][mz.slot][mz,,drop =FALSE]
  tran_counts <- data[["pathway"]]["counts"]

  gene_mappings <- data.frame(gene = rownames(tran_counts))

  rownames(tran_counts) <- unlist(lapply(1:length(rownames(tran_counts)), function (x) {
    paste0("mz-",(round(as.numeric(gsub("mz-", "", rownames(met_counts)[length(rownames(met_counts))],))) + 100), x)
  }))

  gene_mappings$mz <- rownames(tran_counts)
  gene_mappings$raw_mz <- gsub("mz-", "", gene_mappings$mz)

  data[["tmp"]] <- SeuratObject::CreateAssay5Object(counts = rbind(met_counts, tran_counts))
  data[["tmp"]][mz.slot] <- data[["tmp"]]["counts"]
  SM.assay <- "tmp"

  gc()

  main_cardinal <- ConvertSeuratToCardinal(data = data, assay = SM.assay, slot = mz.slot, verbose = FALSE)



  annotation_ids <- unlist(strsplit(row$all_Isomers_IDs, split = "; "))
  ramp_ids <- unlist(lapply(annotation_ids, function(analyte){
    source_df_copy[source_df_copy$sourceId == analyte, "rampId"]
  }))

  if (length(ramp_ids) > 1){
    key_pathways <- unique(unlist(named_list_copy[ramp_ids]))

    if (length(key_pathways) > 0){

      present_gene_mappings <- gene_mappings[gene_mappings$gene %in% key_pathways,]
      features <- c(row$raw_mz,present_gene_mappings$raw_mz)
      cardinal_subset <- Cardinal::subset(main_cardinal, mz %in% features)

      d_list <- list()
      d_list[["1"]] <- suppressWarnings(Cardinal::colocalized(cardinal_subset, mz=row$raw_mz, n = length(features)))
      for (i in names(d_list)){
        correlated_features <- d_list[[i]][order(d_list[[i]]$mz), ]
        correlated_features$features <- c(row$mz_names,present_gene_mappings$gene)
        correlated_features$modality <- c("metabolite", c(rep("gene", length(present_gene_mappings$gene))))
        correlated_features <- correlated_features[c("features", colnames(correlated_features)[!colnames(correlated_features) %in% c("mz", "features")])]
        d_list[[i]] <- correlated_features
      }

      for (i in names(d_list)){
        correlated_features <- d_list[[i]]
        correlated_features <- correlated_features[order(-abs(correlated_features$cor)), ]
        correlated_features$ident <- i
        correlated_features$rank <- c(1:length(correlated_features$ident))
        d_list[[i]] <- correlated_features

      }

      correlated_features <- data.frame(do.call(rbind, d_list))


      correlation_results <- list()
      for (annotation in ramp_ids){
        matched_pathways <- named_list_copy[[annotation]]
        if (!is.null(matched_pathways)){
          df <- correlated_features[correlated_features$features %in% matched_pathways, ]
          if (nrow(df) > 0){
            max_row <- df %>% slice_max(order_by = abs(cor), n = 1)
            df <- df[df$cor > corr_theshold,]
            correlation_results[[annotation]] <- list(
              "max_cor" = unique(pull(max_row, cor)),
              "max_path" = pull(max_row, features),
              "n_sig_path" = nrow(df)
            )
          }
        } else{
          correlation_results[[annotation]] <- list(
            "max_cor" = 0,
            "max_path" = "NA",
            "n_sig_path" = 0
          )
        }
      }
      metabolite_scores <- tibble::tibble(
        ramp_id = names(correlation_results),
        max_cor = sapply(correlation_results, function(x) x$max_cor),
        n_sig_path = sapply(correlation_results, function(x) x$n_sig_path))

      metabolite_scores <- metabolite_scores %>%
        mutate(
          abs_cor = abs(max_cor),
          z_cor = scale(abs_cor)[,1],
          z_path = scale(n_sig_path)[,1],
          z_score = (corr_weight *z_cor) + (n_weight* z_path)  # or use mean instead of sum
        )

      metabolite_scores <- metabolite_scores %>%
        mutate(
          pval = 1 - pnorm(z_score)
        )
      metabolite_scores <- metabolite_scores %>%
        mutate(
          pval_adj = p.adjust(pval, method = "BH")
        )

      metabolite_scores <- metabolite_scores[c("ramp_id","max_cor","n_sig_path", "z_score", "pval", "pval_adj")]
      metabolite_scores

    } else {
      stop("Provided m/z value contained annotations with no corresponding pathways ....")
    }
  } else {
    stop("Provided m/z value contained only 1 or less annotations ....")
  }
  message("Calculating Statistics")
  column_names <- c("ramp_id","max_cor","n_sig_path", "z_score", "pval", "pval_adj")

  if (is.null(metabolite_scores) || nrow(metabolite_scores) == 0) {
    return(tibble(ramp_id = NA, max_cor = NA, n_sig_path = NA, z_score = NA, pval = NA, pval_adj = NA, metabolite = NA))
  } else {
    get_most_common_name <- function(ramp_id) {
      if (is.na(ramp_id)) return(NA)

      matches <- source_df_copy[source_df_copy$rampId == ramp_id, ]
      if (nrow(matches) == 0) return(NA)

      name_counts <- sort(table(matches$commonName), decreasing = TRUE)
      return(names(name_counts)[1])
    }

    # Apply to df
    metabolite_scores$metabolite <- sapply(metabolite_scores$ramp_id, get_most_common_name)
    return(metabolite_scores[c("metabolite","ramp_id", "n_sig_path", "max_cor", "z_score", "pval", "pval_adj")])
  }

  return(metabolite_scores)

}


#' Calculate annotation statistics for all m/z value suggesting the most likely metabolite based on correlated pathway expression.
#'
#' This function evaluates the potential biological relevance of an annotation for a given m/z value by:
#' - Identifying pathways associated with each possible annotated metabolite
#' - Calculating colocalisation score between the m/z intensity and the expression of each corresponding pathway
#' - Ranking annotations by a combined z-score based on correlation strength and number of supporting significant pathways
#' - NOTE: this function requires `CreatePathwayAssay` and `CreatePathwayObject` to be run first
#'
#' @param mz Character or numeric. The target m/z feature. If numeric, the closest matching m/z in the dataset will be selected.
#' @param data A SpaMTP Seurat object containing both metabolite and RAMP_ID assays generated from `CreatePathwayAssay`.
#' @param mz.assay Character string defining the name of the assay containing m/z features.
#' @param pathway.assay Character string matching the name of the assay containing RAMP_ID features (default = "pathway").
#' @param mz.slot Character string stating the slot to extract m/z values from (default = "scale.data").
#' @param pathway.slot Character string defining the slot to extract pathway features from (default = "scale.data").
#' @param return.top Boolean indicating whether to return only the most likely metabolite with it's corresponding pval and score. If set to `FALSE`, a list will be returned with statistics for all possible metabolites per m/z (default = TRUE).
#' @param corr_theshold Numeric value stating the correlation threshold to consider a pathway as significantly colocalized. If set to `0`, all pathways will be counted (default = 0).
#' @param corr_weight Numeric weight applied to correlation score in z-score calculation. If significance should be based more on the correlation, increase this value (default = 1).
#' @param n_weight Numeric weight applied to number of correlated pathways in z-score calculation (default = 1).
#'
#' @return Either a data.frame containing the original annotations for all m/z values and their corresponding most likely metabolite, or a list contating statistics for each m/z value.
#'
#' @importFrom Cardinal subset colocalized
#' @importFrom dplyr group_by summarise mutate slice_max pull
#' @importFrom tibble tibble
#' @importFrom SeuratObject CreateAssay5Object
#' @importFrom stats pnorm p.adjust
#' @importFrom utils head
#'
#' @examples
#' #CalculateAnnotationStatistics(data = data,mz.assay = "SPM",pathway.assay = "merged",mz.slot = "scale.data")
#' }
#'
#' @export
CalculateAnnotationStatistics <- function(data, mz.assay, pathway.assay, mz.slot= "scale.data", pathway.slot = "scale.data", return.top = TRUE, corr_theshold = 0, corr_weight = 1, n_weight = 1){

  message("creating pathway expression object")
  y <- CreatePathwayObject(data,
                           assay=pathway.assay,
                           slot = pathway.slot)


  ## 1. get pathways for each ramp_id
  DefaultAssay(y) <- "pathway"
  pathway_list <- analytehaspathway %>%
    group_by(rampId) %>%
    summarise(pathways = list(unique(pathwayRampId)), .groups = "drop")

  named_list <- setNames(pathway_list$pathways, pathway_list$rampId)

  ## 2. Use source_df to match annotation sourceID to rampID

  meta_data <- y[[mz.assay]]@meta.data$raw_mz
  meta_rows <- y[[mz.assay]]@meta.data
  source_df_copy <- source_df
  named_list_copy <- named_list

  message("Calculating top metabolites for each m/z value. ")


  met_counts <- y[[mz.assay]][mz.slot] #[FindNearestMZ(data = y, target_mz = mz, assay = SM.assay),,drop =FALSE]
  tran_counts <- y[["pathway"]]["counts"]

  gene_mappings <- data.frame(gene = rownames(tran_counts))

  rownames(tran_counts) <- unlist(lapply(1:length(rownames(tran_counts)), function (x) {
    paste0("mz-",(round(as.numeric(gsub("mz-", "", rownames(met_counts)[length(rownames(met_counts))],))) + 100), x)
  }))

  gene_mappings$mz <- rownames(tran_counts)
  gene_mappings$raw_mz <- gsub("mz-", "", gene_mappings$mz)

  y[["tmp"]] <- SeuratObject::CreateAssay5Object(counts = rbind(met_counts, tran_counts))
  y[["tmp"]][mz.slot] <- y[["tmp"]]["counts"]
  SM.assay <- "tmp"

  gc()

  main_cardinal <- ConvertSeuratToCardinal(data = y, assay = SM.assay, slot = mz.slot, verbose = FALSE)


  mz_pathway_annotations <- lapply(1:length(meta_data), function(idx){
    row <- meta_rows[idx, ]
    annotation_ids <- unlist(strsplit(row$all_Isomers_IDs, split = "; "))
    ramp_ids <- unlist(lapply(annotation_ids, function(analyte){
      source_df_copy[source_df_copy$sourceId == analyte, "rampId"]
    }))

    if (length(ramp_ids) > 1){
      key_pathways <- unique(unlist(named_list_copy[ramp_ids]))

      if (length(key_pathways) > 0){

        present_gene_mappings <- gene_mappings[gene_mappings$gene %in% key_pathways,]
        features <- c(row$raw_mz,present_gene_mappings$raw_mz)
        cardinal_subset <- Cardinal::subset(main_cardinal, mz %in% features)

        d_list <- list()
        d_list[["1"]] <- suppressWarnings(Cardinal::colocalized(cardinal_subset, mz=row$raw_mz, n = length(features)))
        for (i in names(d_list)){
          correlated_features <- d_list[[i]][order(d_list[[i]]$mz), ]
          correlated_features$features <- c(row$mz_names,present_gene_mappings$gene)
          correlated_features$modality <- c("metabolite", c(rep("gene", length(present_gene_mappings$gene))))
          correlated_features <- correlated_features[c("features", colnames(correlated_features)[!colnames(correlated_features) %in% c("mz", "features")])]
          d_list[[i]] <- correlated_features
        }

        for (i in names(d_list)){
          correlated_features <- d_list[[i]]
          correlated_features <- correlated_features[order(-abs(correlated_features$cor)), ]
          correlated_features$ident <- i
          correlated_features$rank <- c(1:length(correlated_features$ident))
          d_list[[i]] <- correlated_features

        }

        correlated_features <- data.frame(do.call(rbind, d_list))


        correlation_results <- list()
        for (annotation in ramp_ids){
          matched_pathways <- named_list_copy[[annotation]]
          if (!is.null(matched_pathways)){
            df <- correlated_features[correlated_features$features %in% matched_pathways, ]
            if (nrow(df) > 0){
              max_row <- df %>% slice_max(order_by = abs(cor), n = 1)
              df <- df[df$cor > corr_theshold,]
              correlation_results[[annotation]] <- list(
                "max_cor" = unique(pull(max_row, cor)),
                "max_path" = pull(max_row, features),
                "n_sig_path" = nrow(df)
              )
            }
          } else{
            correlation_results[[annotation]] <- list(
              "max_cor" = 0,
              "max_path" = "NA",
              "n_sig_path" = 0
            )
          }
        }
        metabolite_scores <- tibble::tibble(
          ramp_id = names(correlation_results),
          max_cor = sapply(correlation_results, function(x) x$max_cor),
          n_sig_path = sapply(correlation_results, function(x) x$n_sig_path))

        metabolite_scores <- metabolite_scores %>%
          mutate(
            abs_cor = abs(max_cor),
            z_cor = scale(abs_cor)[,1],
            z_path = scale(n_sig_path)[,1],
            z_score = (corr_weight *z_cor) + (n_weight* z_path)  # or use mean instead of sum
          )

        metabolite_scores <- metabolite_scores %>%
          mutate(
            pval = 1 - pnorm(z_score)
          )
        metabolite_scores <- metabolite_scores %>%
          mutate(
            pval_adj = p.adjust(pval, method = "BH")
          )

        metabolite_scores <- metabolite_scores[c("ramp_id","max_cor","n_sig_path", "z_score", "pval", "pval_adj")]
        return(metabolite_scores)

      }
      return(NULL)
    }
    return(NULL)
  })

  message("Calculating Statistics")
  column_names <- c("ramp_id","max_cor","n_sig_path", "z_score", "pval", "pval_adj")

  if (return.top){

    df <- map_dfr(mz_pathway_annotations, function(df) {
      if (is.null(df) || nrow(df) == 0) {
        # Return a 1-row tibble of NAs if NULL or empty
        tibble(ramp_id = NA, max_cor = NA, n_sig_path = NA, z_score = NA,pval = NA, pval_adj = NA)
      } else {
        df %>%
          select(all_of(column_names)) %>%
          slice_min(pval_adj, with_ties = FALSE)
      }
    }, .id = "source")

    get_most_common_name <- function(ramp_id) {
      if (is.na(ramp_id)) return(NA)

      matches <- source_df_copy[source_df_copy$rampId == ramp_id, ]
      if (nrow(matches) == 0) return(NA)

      name_counts <- sort(table(matches$commonName), decreasing = TRUE)
      return(names(name_counts)[1])
    }

    # Apply to df
    df$metabolite <- sapply(df$ramp_id, get_most_common_name)
    df$ramp_id <- NULL
    df <- df[c("metabolite","n_sig_path", "max_cor", "z_score", "pval", "pval_adj")]

    combined_df <- bind_cols(meta_rows[c(1,2)], df)
    combined_df$original_annotations <- meta_rows$all_IsomerNames

    return(combined_df)
  } else {
    warning("Returning results for all annotations per m/z. For returning a data.frame with only the most significant metabolite per m/z set `return.type` == 'top' ")

    result_list <- lapply(mz_pathway_annotations, function(x){
      if (is.null(x) || nrow(x) == 0) {
        return(tibble(ramp_id = NA, max_cor = NA, n_sig_path = NA, z_score = NA, pval = NA, pval_adj = NA, metabolite = NA))
      } else {
        get_most_common_name <- function(ramp_id) {
          if (is.na(ramp_id)) return(NA)

          matches <- source_df_copy[source_df_copy$rampId == ramp_id, ]
          if (nrow(matches) == 0) return(NA)

          name_counts <- sort(table(matches$commonName), decreasing = TRUE)
          return(names(name_counts)[1])
        }

        # Apply to df
        x$metabolite <- sapply(x$ramp_id, get_most_common_name)
        return(x[c("metabolite","ramp_id", "n_sig_path", "max_cor", "z_score", "pval", "pval_adj")])
      }
    })

    names(result_list) <- meta_rows$mz_names
    return(result_list)

  }
}
