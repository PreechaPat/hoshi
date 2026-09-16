"""
hoshi convert — Convert between metagenomics formats.

Reads a single sample input file/directory in one format (Emu TSV or a Savont
output directory) via the SummarizedExperiment intermediate representation, then
exports to a target format (Kraken2 report or a flat species count table).

The ``table`` output addresses the gap in Savont's native species output, which
lists relative abundance and taxonomy but drops both the NCBI ``tax_id`` and the
estimated read count. ``convert --output-format table`` re-emits the species
table with those two columns restored.

Example usage:
    hoshi convert test_data/emu_output/test_batch/output/sample01_rel-abundance.tsv \
        -o sample01_kraken2.txt

    hoshi convert --input-format savont --output-format table \
        test_data/savont_output/test_ind/savont-out-sample01 \
        -o sample01_species_counts.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from hoshi.lib.egress import count_table_to_tsv, experiment_to_kraken2
from hoshi.lib.ingress import (
    read_emu_abundance_into_summarizedexperiment,
    read_savont_abundance_into_summarizedexperiment,
)

_SUPPORTED_INPUT_FORMATS = ("emu", "savont")
_SUPPORTED_OUTPUT_FORMATS = ("kraken2", "table")


def run(args: argparse.Namespace) -> int:
    input_format = args.input_format
    output_format = args.output_format
    input_path = Path(args.input_file)
    output = Path(args.output) if args.output else None

    # Validate input path exists (a file for emu, a directory for savont)
    if not input_path.exists():
        print(f"Error: input path not found: {input_path}", file=sys.stderr)
        return 1

    # Load into SummarizedExperiment (single sample)
    if input_format == "emu":
        experiment = read_emu_abundance_into_summarizedexperiment(input_path)
    elif input_format == "savont":
        sample_name = args.name or input_path.name
        experiment = read_savont_abundance_into_summarizedexperiment(
            input_path, sample_names=[sample_name]
        )
    else:
        print(
            f"Error: unsupported input format '{input_format}'. "
            f"Supported: {', '.join(_SUPPORTED_INPUT_FORMATS)}",
            file=sys.stderr,
        )
        return 1

    # Convert and output
    if output_format == "kraken2":
        rendered = experiment_to_kraken2(experiment)
    elif output_format == "table":
        rendered = count_table_to_tsv(experiment)
    else:
        print(
            f"Error: unsupported output format '{output_format}'. "
            f"Supported: {', '.join(_SUPPORTED_OUTPUT_FORMATS)}",
            file=sys.stderr,
        )
        return 1

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
        print(f"Written: {output}")
    else:
        sys.stdout.write(rendered)

    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "convert",
        help="Convert a single sample between metagenomics output formats.",
        description=(
            "Convert a single-sample metagenomics classification output "
            "between formats. Reads input via a SummarizedExperiment "
            "intermediate, enabling conversion from any supported input "
            "to any supported output format. Use '--output-format table' to "
            "export a flat species table with tax_id and estimated counts "
            "(the columns Savont's native species output omits)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        help=(
            "Input to convert (single sample): an Emu rel-abundance TSV, or a "
            "Savont output directory when '--input-format savont'."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file path. If omitted, writes to stdout.",
    )
    parser.add_argument(
        "-n", "--name",
        help=(
            "Sample name (savont input only). Defaults to the directory name."
        ),
    )
    parser.add_argument(
        "--input-format",
        choices=_SUPPORTED_INPUT_FORMATS,
        default="emu",
        help="Input format (default: emu).",
    )
    parser.add_argument(
        "--output-format",
        choices=_SUPPORTED_OUTPUT_FORMATS,
        default="kraken2",
        help="Output format (default: kraken2).",
    )

    parser.set_defaults(func=run)
    return parser
