# chem_props: A database containing the chemical properties and metadata of each RAMP_DB analyte

Pruned from RaMP-DB 3.0.7. This object contains chemical structures,
monoisotopic masses, names, identifiers, and formulae. It can be passed
directly to
[`BuildMZAnnotationIndex()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/BuildMZAnnotationIndex.md)
or
[`AnnotateMZ()`](https://genomicsmachinelearning.github.io/SpaMTP/reference/AnnotateMZ.md).

## Usage

``` r
chem_props
```

## Format

### A data frame with 289,754 rows and 11 variables:

- ramp_id:

  RaMP analyte ID (character)

- chem_data_source:

  Relative source database for the respective analyte (character)

- chem_source_id:

  Relative source database ID for analyste (character)

- iso_smiles:

  Smile structure for relative analyte (character)

- inchi_key_prefix:

  InChlKey prefix of Internation Chemical Identifier (InChl) (character)

- inchi_key:

  Full InChlKey for corresponding analyte (character)

- inchi:

  Full InChl identifier structure for analyte (character)

- mw:

  Molecular weight for analyte (double)

- monoisotop_mass:

  Relative monoisotopic mass for analyate (double)

- common_name:

  Analytes common name (character)

- mol_formula:

  Analytes simplified molecular fomula (character)

## Source

<https://github.com/ncats/RaMP-DB>
