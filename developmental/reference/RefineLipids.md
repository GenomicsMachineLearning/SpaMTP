# Uses common lipid nomenclature to simplify lipid annotations

Based on rgoslin, this function has the ability to simplify lipid
annotations into more general terms, groups and classes. This function
can input a data frame containing a column storing the annotation and
refine only the entries which are annotated as known lipid structures.

## Usage

``` r
RefineLipids(
  data,
  annotation.column = "annotations",
  database = "HMDB",
  lipid_info = "simple"
)
```

## Arguments

- data:

  data.frame containing a column with lipid annotations

- annotation.column:

  Character string matching the column name containing the lipid
  annotations in the provided df (default = "annotations").

- database:

  Character string defining the database matching the lipid annotations.
  Possible entries include
  c('Shorthand2020','Goslin','FattyAcids','LipidMaps','SwissLipids','HMDB')
  (default = "HMDB").

- lipid_info:

  Character string defining what level of information to return for each
  lipid. Options are either "all" or "simple". "Simple" returns a
  smaller, simplified list of annotations for each lipid. Please visit
  https://bioconductor.org/packages/release/bioc/vignettes/rgoslin/inst/doc/introduction.html
  to see all possible outputs (default = "simple")).

## Value

Data.frame containing additional columns with simplified lipid names

## Examples

``` r
# RefineLipids(DEPs_df)
```
