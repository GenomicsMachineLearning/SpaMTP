# Calculate annotation statistics for a single m/z value suggesting the most likely metabolite based on correlated pathway expression.

This function evaluates the potential biological relevance of an
annotation for a given m/z value by:

- Identifying pathways associated with each possible annotated
  metabolite

- Calculating colocalisation score between the m/z intensity and the
  expression of each corresponding pathway

- Ranking annotations by a combined z-score based on correlation
  strength and number of supporting significant pathways

- NOTE: this function requires `CreatePathwayAssay` and
  `CreatePathwayObject` to be run first

## Usage

``` r
CalculateSingleAnnotationStatistics(
  mz,
  data,
  mz.assay,
  pathway.assay = "pathway",
  mz.slot = "scale.data",
  pathway.slot = "scale.data",
  corr_theshold = 0,
  corr_weight = 1,
  n_weight = 1
)
```

## Arguments

- mz:

  Character or numeric. The target m/z feature. If numeric, the closest
  matching m/z in the dataset will be selected.

- data:

  A SpaMTP Seurat object containing both metabolite and pathway assays
  generated from `CreatePathwayObject`.

- mz.assay:

  Character string defining the name of the assay containing m/z
  features.

- pathway.assay:

  Character string matching the name of the assay containing pathway
  features (default = "pathway").

- mz.slot:

  Character string stating the slot to extract m/z values from (default
  = "scale.data").

- pathway.slot:

  Character string defining the slot to extract pathway features from
  (default = "scale.data").

- corr_theshold:

  Numeric value stating the correlation threshold to consider a pathway
  as significantly colocalized. If set to `0`, all pathways will be
  counted (default = 0).

- corr_weight:

  Numeric weight applied to correlation score in z-score calculation. If
  significance should be based more on the correlation, increase this
  value (default = 1).

- n_weight:

  Numeric weight applied to number of correlated pathways in z-score
  calculation (default = 1).

## Value

A tibble with the ranked annotations for the m/z value, containing:

- metabolite:

  Most common name of the metabolite associated with the RAMP ID

- ramp_id:

  The RAMP ID corresponding to the annotation

- n_sig_path:

  Number of correlated pathways above the threshold

- max_cor:

  Maximum correlation value among significant pathways

- z_score:

  Combined z-score used to rank annotations

- pval:

  Unadjusted p-value

- pval_adj:

  Adjusted p-value (BH method)

## Examples

``` r
#data <- CreatePathwayObject(data,assay="SPT_pathway",slot = "scale.data")
#CalculateSingleAnnotationStatistics(mz = "mz-674.2805",data = data,mz.assay = "SPM",pathway.assay = "pathway",mz.slot = "scale.data")
```
