# Hoshi — 16S rRNA Sequencing Report Generator

Hoshi turns the output of 16S rRNA taxonomic classifiers into readable HTML/PDF
reports. Heterogeneous tool outputs (EMU, Savont, …) are read into a single
internal representation — a Bioconductor-style `SummarizedExperiment` — and every
report renders from that common shape, so the reporting layer never depends on
which classifier produced the data.

## Features

- **Multiple input formats** — Savont output directories (default) and EMU
  `*_rel-abundance.tsv` files, behind one `AbundanceReader` interface.
- **Single-sample report** — responsive HTML with a species table, alpha
  diversity, and an interactive taxonomy Sankey; optional A4 PDF.
- **Multi-sample report** — one combined HTML with a tab per sample.
- **Medical detection report** — clinical-style HTML/PDF for 16S bacterial
  detection, with optional patient/QC/interpretation metadata and pathogen
  classification.
- **Format conversion** — export a sample to Kraken2 report format.
- **Taxonomy enrichment** — fill lineage columns for a `tax_id` table from an
  NCBI taxdb.

## Requirements

- Python **3.12+**
- [`uv`](https://docs.astral.sh/uv/) for dependency management (recommended)
- Runtime deps: `jinja2`, `pandas`, `plotly`, `kaleido`, `taxopy`, `weasyprint`
  (installed automatically)

> PDF output uses **WeasyPrint**, which relies on native libraries
> (`libpango`, `libcairo`, `libgdk-pixbuf`). On Debian/Ubuntu:
> `apt-get install libpango-1.0-0 libpangocairo-1.0-0 libgdk-pixbuf-2.0-0`.
> A `Dockerfile` is provided if you'd rather not install these system-wide.

## Setup

```bash
# From the project root:
uv sync                       # create the venv and install deps (incl. dev group)
uv run hoshi --help           # verify the CLI

# Or, without uv:
pip install -e .
hoshi --help
```

## Quickstart

### Single-sample report

```bash
# Savont output directory (default input format)
uv run hoshi report-single test_data/savont_output/test_ind/savont-out-sample01 \
    -n sample01 -o dist/sample01_report.html

# EMU output directory, with a PDF alongside the HTML
uv run hoshi report-single --input-format emu \
    test_data/emu_output/test_ind/sample01 \
    -n sample01 -o dist/sample01_report.html --pdf
```

### Multi-sample report

```bash
uv run hoshi report-multi \
    test_data/savont_output/test_ind/savont-out-sample01 \
    test_data/savont_output/test_ind/savont-out-soil \
    -o dist/multi_report.html
```

### Medical detection report

```bash
# Savont directory + optional clinical metadata JSON
uv run hoshi report-medical \
    test_data/savont_output/test_ind/savont-out-sample01 \
    -m test_data/medical/metadata_only.json \
    -o dist/sample01_medical.html --pdf

# EMU TSV file as input
uv run hoshi report-medical --input-format emu \
    test_data/emu_output/emu-mock01.tsv \
    -o dist/mock01_medical.html
```

### Convert to Kraken2 format

```bash
uv run hoshi convert test_data/emu_output/emu-mock01.tsv -o dist/mock01_kraken2.txt
```

### Convert Savont output to a species count table

Savont's native species output lists relative abundance and taxonomy but drops
both the NCBI `tax_id` and the estimated read count. This re-emits the species
table (one row per species) with those two columns restored.

```bash
uv run hoshi convert --input-format savont --output-format table \
    test_data/savont_output/test_ind/savont-out-sample01 \
    -o dist/sample01_species_counts.tsv
```

Output columns: `relative_abundance`, `estimated_count`, `tax_id`, `species`,
`genus`, `family`, `order`, `class`, `phylum`, `superkingdom`.

### Enrich a tax_id table with lineage

```bash
uv run hoshi enrich my_table.tsv -o my_table_enriched.tsv
```

## Commands

| Command          | Purpose                                          | Input                                 | Output      |
|------------------|--------------------------------------------------|---------------------------------------|-------------|
| `report-single`  | One responsive HTML report for a single sample   | Savont dir (default) or EMU dir       | HTML (+PDF) |
| `report-multi`   | Combined report, one tab per sample              | Multiple Savont/EMU dirs              | HTML (+PDF) |
| `report-medical` | Clinical 16S bacterial detection report          | Savont dir or EMU TSV + optional JSON | HTML (+PDF) |
| `convert`        | Convert a sample between formats                 | EMU TSV or Savont dir                 | Kraken2 or species count table |
| `enrich`         | Add taxonomy lineage columns to a `tax_id` table | Table with a `tax_id` column          | TSV         |

Run `uv run hoshi <command> --help` for the full flag list of any command.

### Common flags

- `--input-format {savont,emu}` — classifier input format (`report-single`,
  `report-multi`, `report-medical`). Default: `savont`.
- `-o, --output` — output path. `report-multi` prints HTML to stdout when omitted.
- `--pdf` / `--output-pdf PATH` — also render a PDF via WeasyPrint. For
  `report-multi`, `--pdf` requires `-o/--output`.
- `-n, --name` / `--title` — sample name / page title (`report-single`).

### `report-medical` specifics

- `-m, --metadata PATH` — clinical metadata JSON (patient/specimen IDs, QC
  items, conclusion, lab identity, …). Any omitted field renders as `N/A`; an
  omitted `authorized_by` renders as a blank signature line.
- `--pathogen-sheet PATH` — CSV mapping NCBI `tax_id` → pathogen class
  (commensal / potential / opportunistic / primary). Default: the pathogen
  sheet bundled with hoshi (`src/hoshi/assets/pathogen_sheet.csv`). Pass an
  empty value to disable the column.
- `--top N` — number of top organisms to include (default: 5).

## Architecture at a glance

- **`lib/experiment.py`** — `SummarizedExperiment`: aligned assay matrices
  (`abundance`, `counts`), `row_data` (taxonomy), `col_data` (sample metadata),
  and free-form `metadata`. Validated on construction; transformations return new
  instances.
- **`lib/reader.py`** — `AbundanceReader` protocol + `build_reader(...)` factory.
  Concrete readers: `EmuReader`, `SavontReader`. Adding a classifier means adding
  a reader; reports are untouched.
- **`lib/report.py`** — `Report` composite: a `SummarizedExperiment` plus
  report-time data (e.g. per-species confidence, clinical metadata).
- **`command/`** — thin CLI glue + Jinja2 report generation.

```
src/hoshi/
├── cli.py                     # argparse entry point; registers sub-commands
├── lib/
│   ├── experiment.py          # SummarizedExperiment (internal standard)
│   ├── reader.py              # AbundanceReader protocol + build_reader
│   ├── emu_reader.py          # EMU reader
│   ├── savont_reader.py       # Savont reader
│   ├── ingress*.py            # abundance ingestion helpers
│   ├── egress.py              # experiment_to_kraken2 / experiment_to_count_table
│   ├── report.py              # Report composite (SE + report-time data)
│   ├── medical.py             # medical report builder / view model
│   ├── pathogen.py            # pathogen sheet lookup
│   ├── sankey.py              # Sankey data + figure
│   ├── diversity.py           # alpha diversity
│   └── taxdb.py               # NCBI lineage lookup (enrich)
└── command/                   # thin CLI glue + report generation
    ├── convert.py  enrich.py
    ├── report_single.py  report_multi.py  report_medical.py
    └── templates/
        ├── microbiome/        # single / multi templates
        └── medical/           # medical report template
```

## Development

```bash
uv run pytest                          # run the test suite
uv run pytest tests/test_ingress.py    # focused run
uv run ruff check src tests            # lint
uv run ruff format                     # auto-format
```

Conventions (see [`AGENTS.md`](AGENTS.md) for the full guide):

- PEP 8, 4-space indents, type hints where practical.
- Tests mirror module names: `src/hoshi/lib/ingress.py` → `tests/test_ingress.py`.
- Reference inputs live under `test_data/`; keep new samples lightweight.
- Keep Jinja templates self-documenting (`{# ... #}`) and align context keys with
  DataFrame column names.
- Conventional commits: `feat: add species coverage table`.

## Known gaps

- `pyproject.toml` declares `hoshi-report` and `hoshi-summary` console scripts
  pointing at `hoshi.command.report:main` / `hoshi.command.summary:main`, which
  do not currently exist. Use the `hoshi` entry point; these two are stale and
  should be removed or implemented.
- Report templates are functional but minimal ("barebone"). Enriching them with
  more information — QC metrics, run/read stats from other tools — is the "complete
  report" goal tracked in [`SPECS.md`](SPECS.md).
- Bracken / QIIME2 / BIOM readers and the comparison report are planned, not
  implemented.

## License

See [`LICENSE`](LICENSE).
