# Calculate annotation statistics for all m/z value suggesting the most likely metabolite based on correlated pathway expression.

This function evaluates the potential biological relevance of an
annotation for a given m/z value by:

- Identifying pathways associated with each possible annotated
  metabolite

- Calculating colocalisation score between the m/z intensity and the
  expression of each corresponding pathway

- Ranking annotations by a combined z-score based on correlation
  strength and number of supporting significant pathways

- NOTE: this function requires `CreatePathwayAssay` to be run first

## Usage

``` r
CalculateAnnotationStatistics(
  data,
  mz.assay,
  pathway.assay,
  mz.slot = "scale.data",
  pathway.slot = "scale.data",
  return.top = TRUE,
  corr_theshold = 0,
  corr_weight = 1,
  n_weight = 1
)
```

## Arguments

- data:

  A SpaMTP Seurat object containing both metabolite and RAMP_ID assays
  generated from `CreatePathwayAssay`.

- mz.assay:

  Character string defining the name of the assay containing m/z
  features.

- pathway.assay:

  Character string matching the name of the assay containing RAMP_ID
  features (default = "pathway").

- mz.slot:

  Character string stating the slot to extract m/z values from (default
  = "scale.data").

- pathway.slot:

  Character string defining the slot to extract pathway features from
  (default = "scale.data").

- return.top:

  Boolean indicating whether to return only the most likely metabolite
  with it's corresponding pval and score. If set to `FALSE`, a list will
  be returned with statistics for all possible metabolites per m/z
  (default = TRUE).

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

- mz:

  Character or numeric. The target m/z feature. If numeric, the closest
  matching m/z in the dataset will be selected.

## Value

Either a data.frame containing the original annotations for all m/z
values and their corresponding most likely metabolite, or a list
contating statistics for each m/z value.

## Examples

``` r
#CalculateAnnotationStatistics(data = data,mz.assay = "SPM",pathway.assay = "merged",mz.slot = "scale.data")
```
