# HOSHI
HOSHI is a report generator that builds reports from 16S rRNA taxonomic
classifier output (relative abundance and/or counts). It works by
reading output from a taxonomic classifier (taxonomic abundance — see
`test_data/*_output/`, usually a tabular format) and combining it with
metadata such as patient name into a report.

## Scope

### Overall
- Classifier output (EMU, Savont) is read by an `AbundanceReader`
  (`src/hoshi/lib/reader.py`) into a `SummarizedExperiment`
  (`src/hoshi/lib/experiment.py`). Reports render from that shape only, so
  adding a classifier means adding a reader, not touching the reports.
- Convert taxonomic output to Kraken2 report format or a species count table
  (`hoshi convert`).
- Enrich a `tax_id` table with taxonomy lineage from an NCBI taxdb
  (`hoshi enrich`).

### Medical report
Clinical-style 16S bacterial detection report (`hoshi report-medical`) with
patient/QC/interpretation metadata and pathogen classification. This report
annotates bacteria based on their `tax_id` using the
[pathogen sheet](#pathogen-sheet).

#### Pathogen sheet
`src/hoshi/assets/pathogen_sheet.csv` is the source of truth for classifying
organisms into pathogen classes (commensal / potential / opportunistic /
primary), matched on its `taxid` column. Most entries are bacteria; the sheet
also lists some viruses, fungi and protists.

### Microbiome report
Single-sample (`hoshi report-single`) and multi-sample (`hoshi report-multi`,
one tab per sample) reports with a species table and a taxonomy Sankey. The
single-sample report also shows alpha diversity.

### Environmental report
TODO: Being drafted

## Coding style
- PEP 8, 4-space indents, type hints where practical.
- Avoid line-by-line comments unless the code is genuinely technical.
- Prefer a short comment on what a function/module is for and how it works
  over comments that narrate technical detail.
- Keep Jinja templates self-documenting (`{# ... #}`) and align context keys
  with DataFrame column names.
- Conventional commits, e.g. `feat: add species coverage table`.

## Testing
- Use pytest for testing (`uv run pytest`).
- Tests mirror module names: `src/hoshi/lib/egress.py` → `tests/test_egress.py`.
- Reference inputs live under `test_data/`; keep new samples lightweight.
