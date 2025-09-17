# processed/

Layout for processed SABR outputs intended for quick and predictable loading in R.

## Structure

  - `sequences/` — FASTA/FASTA-like sequence artifacts
  - `tables/`
    - `raw/` — untransformed exported tables (e.g., counts)
    - `derived/` — transformed tables (e.g., relative abundance)
  - `rdata/`
    - `dataframes/` - R objects for dataframes used throughout the project. Contains raw and derived versions.
    - `phyloseq/` — R objects for phyloseq (single or lists)
  - `asv_tables/` — R objects containing tabular artifacts for ASV counts.
    - `raw/` — untransformed exported tables (e.g., counts)
    - `derived/` — transformed tables
  - `metadata/` — metadata objects and dictionaries (taxonomy, cleaned metadata)

## Notes

- Keep raw vs derived separate to avoid confusion in analysis.
- Prefer loading `.rda` from `rdata/` or directly with `load()` for speed inside R.
- CSVs are retained for portability and inspection.
