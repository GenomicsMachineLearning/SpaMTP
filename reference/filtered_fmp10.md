# filtered_fmp10: data.frame containing FMP10+ metabolite mappings

This dataset contains metabolite information for various m/z values
corresponding to metabolite ID's from the HMDB, LipidMaps and ChEMI
databases.

## Usage

``` r
filtered_fmp10
```

## Format

A data frame with the following columns:

- `mass`:

  (numeric) The mass-to-charge ratio (*m/z*) of the metabolite.

- `annotation`:

  (character) The metabolite name or structural description.

- `Adduct`:

  (character) The ionisation adduct associated with the metabolite,
  e.g., `[M+K]`.

- `Formula`:

  (character) The chemical formula of the metabolite, e.g., `C9H16O`.

- `Isomers`:

  (character) A unique identifier for the isomers of the metabolite,
  often linked to external databases.

- `Isomers_IDs`:

  (character) The database-specific IDs for the isomers, such as
  `LIPIDMAPS:LMFA05000118`.

- `Error`:

  (numeric) The mass error or difference between the observed and
  theoretical *m/z*, typically in parts-per-million (ppm).

- `IsomerNames`:

  (character) Names of the isomers for the metabolite.

- `Reference_mz`:

  (numeric) The reference mass-to-charge ratio (*m/z*) used for
  comparison or alignment.
