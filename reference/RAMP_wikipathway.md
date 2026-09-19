# RAMP_wikipathway: A list containing network plot information about pathways from the Wiki database

This object contains a collection of information for each RAMP Wiki
network, including their source, destination, direction, and reaction
type for both proteins and metabolites. Topology sources were retrieved
on 2022-11-01; node identifiers were re-harmonised and integrity-checked
against the bundled RaMP-DB 3.0.7 snapshot.

## Usage

``` r
RAMP_wikipathway
```

## Format

### A list with 10 elements:

- id:

  Wiki pathway identifier (character)

- title:

  Pathway title (character)

- database:

  Source database (character)

- species:

  Species (character)

- protEdges:

  Protein interactions: a data frame with `src`, `dest`, `directed`,
  `reaction_type`, and `source_reaction_type`. Row counts vary by
  pathway.

- protPropEdges:

  Propagated protein interactions: a data frame with `src`, `dest`,
  `directed`, `reaction_type`, and `source_reaction_type`. Row counts
  vary by pathway.

- metabolEdges:

  Metabolite interactions: a data frame with `src`, `dest`, `directed`,
  `reaction_type`, and `source_reaction_type`. Row counts vary by
  pathway.

- metabolPropEdges:

  Propagated metabolite interactions: a data frame with `src`, `dest`,
  `directed`, `reaction_type`, and `source_reaction_type`. Row counts
  vary by pathway.

- mixedEdges:

  Protein-metabolite interactions: a data frame with `src`, `dest`,
  `directed`, `reaction_type`, and `source_reaction_type`. Row counts
  vary by pathway.

- timestamp:

  The date of data extraction (Date)

## Details

Interaction labels and directions were restored from the pinned graphite
archive. Edge tables retain the exact label in `source_reaction_type`;
`spamtp_interaction_repair` records the repair provenance on the
collection. See
[`vignette("Pathway_Database_Integration", package = "SpaMTP")`](https://genomicsmachinelearning.github.io/SpaMTP/articles/Pathway_Database_Integration.md).
