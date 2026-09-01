"""
hoshi convert — Convert between metagenomics formats.

Reads a single sample input file in one format (e.g., Emu) via the
SummarizedExperiment intermediate representation, then exports to a
target format (e.g., Kraken2 report).

Example usage:
    hoshi convert test_data/emu_output/test_batch/output/sample01_rel-abundance.tsv \
        -o sample01_kraken2.txt

    hoshi convert --input-format emu --output-format kraken2 \
        test_data/emu_output/test_batch/output/sample01_rel-abundance.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from hoshi.lib.egress import experiment_to_kraken2
from hoshi.lib.ingress import read_emu_abundance_into_summarizedexperiment

_SUPPORTED_INPUT_FORMATS = ("emu",)
_SUPPORTED_OUTPUT_FORMATS = ("kraken2",)


def run(args: argparse.Namespace) -> int:
    input_format = args.input_format
    output_format = args.output_format
    input_path = Path(args.input_file)
    output = Path(args.output) if args.output else None

    # Validate input file exists
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        return 1

    # Load into SummarizedExperiment (single sample)
    if input_format == "emu":
        experiment = read_emu_abundance_into_summarizedexperiment(input_path)
    else:
        print(
            f"Error: unsupported input format '{input_format}'. "
            f"Supported: {', '.join(_SUPPORTED_INPUT_FORMATS)}",
            file=sys.stderr,
        )
        return 1

    # Convert and output
    if output_format == "kraken2":
        report = experiment_to_kraken2(experiment)
        if output:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(report)
            print(f"Written: {output}")
        else:
            sys.stdout.write(report)
    else:
        print(
            f"Error: unsupported output format '{output_format}'. "
            f"Supported: {', '.join(_SUPPORTED_OUTPUT_FORMATS)}",
            file=sys.stderr,
        )
        return 1

    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "convert",
        help="Convert a single sample between metagenomics output formats.",
        description=(
            "Convert a single-sample metagenomics classification output "
            "between formats. Reads input via a SummarizedExperiment "
            "intermediate, enabling conversion from any supported input "
            "to any supported output format."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        help="Input file to convert (single sample).",
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file path. If omitted, writes to stdout.",
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
