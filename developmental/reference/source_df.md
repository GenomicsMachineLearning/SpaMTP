# source_df: A dataframe containing source information about RAMP_ID analyte used for analysis

Pruned from RaMP-DB 3.0.7. This object maps RaMP analytes to source
identifiers, names, database provenance, and pathway counts.

## Usage

``` r
source_df
```

## Format

### A data frame with 1,051,927 rows and 8 variables:

- sourceId:

  Source database analyte ID (character)

- rampId:

  RaMP analyte ID (character)

- IDtype:

  Identifier type (character)

- geneOrCompound:

  Whether the analyte is a gene or compound (character)

- commonName:

  Common name of analyte (character)

- priorityHMDBStatus:

  Priority level of the analyte (character)

- dataSource:

  Relative source databases where analyte is represented (character)

- pathwayCount:

  Number of pathways analyte is present in (integer)

## Source

<https://github.com/ncats/RaMP-DB>
