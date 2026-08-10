# RAMP_kegg: A list containing network plot information about pathways from the KEGG database

This object contains a collection of information for each RAMP KEGG
network, including their source, destination, direction, and reaction
type for both proteins and metabolites. Topology sources were retrieved
on 2022-11-01; node identifiers were re-harmonised and integrity-checked
against the bundled RaMP-DB 3.0.7 snapshot.

## Usage

``` r
RAMP_kegg
```

## Format

### A list with 10 elements:

- id:

  KEGG pathway identifier (character)

- title:

  Pathway title (character)

- database:

  Source database (character)

- species:

  Species (character)

- protEdges:

  A data frame with 1 row and 4 variables for protein edges, including
  `src`, `dest`, `directed`, and `reaction_type` (various types)

- protPropEdges:

  A data frame with 747 rows and 4 variables for protein-protein
  interactions: `src` (character), `dest` (character), `directed`
  (integer), and `reaction_type` (integer)

- metabolEdges:

  A data frame with 1 row and 4 variables for metabolite edges,
  including `src`, `dest`, `directed`, and `reaction_type` (various
  types)

- metabolPropEdges:

  A data frame with 88 rows and 4 variables for metabolite interactions:
  `src` (character), `dest` (character), `directed` (integer), and
  `reaction_type` (integer)

- mixedEdges:

  A data frame with 311 rows and 4 variables for mixed interactions
  between proteins and metabolites: `src` (character), `dest`
  (character), `directed` (integer), and `reaction_type` (integer)

- timestamp:

  The date of data extraction (Date)
