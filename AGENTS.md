# HOSHI
HOSHI is a report generator that builds reports from 16S rRNA taxonomic
classifier output (relative abundance and/or counts).

## Scope
- Read output from a taxonomic classifier (taxonomic abundance — see
  `test_data/**_output/`, usually a tabular format) and populate it with
  metadata such as patient name.
- Attempt to annotate bacteria based on its `tax_id` when available.
- Convert taxonomic output to Kraken2 report format.

### Medical report
Clinical-style 16S bacterial detection report with patient/QC/interpretation
metadata and pathogen classification.

#### Pathogen sheet
`src/hoshi/assets/pathogen_sheet.csv` is the source of truth for classifying
bacteria into pathogen classes (commensal / potential / opportunistic / primary).

### Environmental report
TODO: Being draft

## Coding style
- Avoid line-by-line comments unless the code is genuinely technical.
  Prefer comments that explain overall functionality over technical detail.
- Short, concise comments on what a function/module is for and how it works.

### Testing
- Use pytest for testing.
