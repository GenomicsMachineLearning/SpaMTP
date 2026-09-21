# Inspect SpaMTP database resources

Inspect SpaMTP database resources

## Usage

``` r
SpaMTPDatabaseInfo(version = NULL)
```

## Arguments

- version:

  Optional bundled database snapshot version. `NULL` and `"latest"` list
  the snapshot shipped with this SpaMTP release.

## Value

A data frame describing bundled resources, including their version,
dimensions, dataset names, and source-package file sizes and MD5
checksums. An unavailable version returns a zero-row data frame.

## Examples

``` r
utils::str(formals(SpaMTPDatabaseInfo))
#> Dotted pair list of 1
#>  $ version: NULL
SpaMTPDatabaseInfo()
#>             resource version  source  category     source_object
#> 1         chem_props   3.0.7 bundled      core        chem_props
#> 2          source_df   3.0.7 bundled      core         source_df
#> 3            analyte   3.0.7 bundled      core           analyte
#> 4  analytehaspathway   3.0.7 bundled      core analytehaspathway
#> 5            pathway   3.0.7 bundled      core           pathway
#> 6   ramp_db_metadata   3.0.7 bundled      core  ramp_db_metadata
#> 7          ramp_hmdb   3.0.7 bundled  topology         RAMP_hmdb
#> 8          ramp_kegg   3.0.7 bundled  topology         RAMP_kegg
#> 9      ramp_reactome   3.0.7 bundled  topology     RAMP_Reactome
#> 10  ramp_wikipathway   3.0.7 bundled  topology  RAMP_wikipathway
#> 11           hmdb_db   3.0.7 bundled    legacy           HMDB_db
#> 12          chebi_db   3.0.7 bundled    legacy          Chebi_db
#> 13      lipidmaps_db   3.0.7 bundled    legacy      Lipidmaps_db
#> 14           gnps_db   3.0.7 bundled    legacy           GNPS_db
#> 15    filtered_fmp10   3.0.7 bundled    legacy    filtered_fmp10
#> 16   smiles_features   3.0.7 bundled structure   smiles_features
#>                    rdata_path    rows columns serialized_bytes
#> 1         data/chem_props.rda  289754      11         12012876
#> 2          data/source_df.rda 1051927       8          9122600
#> 3            data/analyte.rda  463257       2          1204836
#> 4  data/analytehaspathway.rda 1355476       3          1597020
#> 5            data/pathway.rda  122936       5          1930672
#> 6   data/ramp_db_metadata.rda       9      NA             1368
#> 7          data/RAMP_hmdb.rda   48671      NA          1294408
#> 8          data/RAMP_kegg.rda     325      NA           314904
#> 9      data/RAMP_Reactome.rda    2460      NA          2972092
#> 10  data/RAMP_wikipathway.rda     685      NA           236856
#> 11           data/HMDB_db.rda   26538      37          9204653
#> 12          data/Chebi_db.rda   46297      37         10410280
#> 13      data/Lipidmaps_db.rda    9493      37          2588253
#> 14           data/GNPS_db.rda     489      37            96626
#> 15    data/filtered_fmp10.rda    4166       9            71479
#> 16   data/smiles_features.rda  284818      36          4027336
#>                                 md5
#> 1  a87e3ac024a4b526a1492e2a786a3fd1
#> 2  916e2d41d4b265a8b48232c8efb8a663
#> 3  3bbfa9f60ce096a35c03ef95997fd2e9
#> 4  d52f64d24d5ca2c2ed3e5c7ec3181610
#> 5  90207c7b5b389b4de48a3f69f49ce00d
#> 6  a2739146ecd68512acbc27b2ff9b95d2
#> 7  00530e1d781c0cb348c1a7ce28060b58
#> 8  7c0a75a53d1c8a542dcb632157d604f0
#> 9  44467aa8f867a7226df8150556a8ad8c
#> 10 c97ecbd5d1330f5580d38f4a0cff6a06
#> 11 9cb1447f697db32b5f7f3f859e94e6a9
#> 12 4254dd1ec2bd9d434d0daab4f4867baf
#> 13 9984f1e6d1ca4eef858378a30a352bde
#> 14 3b5ff60e6b5c493eb5b73e0216d2ed62
#> 15 39af800b4d76fb707995e9d787e9f441
#> 16 82678726add415f9ee0d5b612b89c404
```
