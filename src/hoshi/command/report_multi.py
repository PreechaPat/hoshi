"""Generate a combined multi-sample HTML report.

Each input is a classifier output *directory* (Savont by default, or EMU via
``--input-format emu``). Every sample is read through the unified reader
interface into a :class:`SummarizedExperiment`, wrapped in a :class:`Report`
(so confidence and other report-time data are layered consistently with the
single-sample and medical reports), and rendered as one tab in the multi-sample
template.

    dirs -> build_reader -> SummarizedExperiment -> Report -> per-sample tab
"""

from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.reader import (
    DEFAULT_INPUT_FORMAT,
    SUPPORTED_INPUT_FORMATS,
    build_reader,
)
from hoshi.lib.report import Report
from hoshi.lib.sankey import get_sankey_data, render_sankey_figure

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_REPORT_TEMPLATE = "microbiome/multisample.html.j2"

# tax_id values used for control/meta rows that must not appear as organisms.
_META_TAX_IDS = {"unmapped", "mapped_filtered", "mapped_unclassified"}


def _dom_id(name: str, index: int) -> str:
    """Build a stable, unique DOM id for a sample tab."""
    slug = re.sub(r"[^a-zA-Z0-9_-]+", "-", name).strip("-").lower()
    return f"tab-{index}-{slug or 'sample'}"


def _build_table(report: Report, name: str, index: int) -> dict:
    """Build one tab's render context from a single-sample :class:`Report`."""
    df = report.experiment.to_dataframe()

    # Drop control/meta rows and sort by abundance for display.
    display = df[~df["tax_id"].astype(str).isin(_META_TAX_IDS)].copy()
    display = display.sort_values("abundance", ascending=False).reset_index(drop=True)

    # Surface per-species confidence (Savont) as a column when present.
    confidence = report.confidence
    cols = ["species", "tax_id", "abundance", "estimated counts"]
    if confidence:
        display["confidence"] = display["tax_id"].astype(str).map(confidence).round(1)
        cols.append("confidence")
    table_html = display[cols].to_html(index=False, border=0, classes="data-table")

    # Per-sample Sankey (interactive HTML).
    try:
        sankey_data = get_sankey_data(df)
        fig = render_sankey_figure(sankey_data, title=f"Taxonomy Flow: {name}")
        sankey_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    except Exception as exc:  # noqa: BLE001 - render best effort
        sankey_html = f"<p>Sankey plot could not be generated: {exc!s}</p>"

    return {
        "dom_id": _dom_id(name, index),
        "name": name,
        "html": table_html,
        "sankey_html": sankey_html,
    }


def run(args: argparse.Namespace) -> int:
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2"), default_for_string=True),
    )

    tables = []
    for index, path_str in enumerate(args.inputs):
        path = Path(path_str)
        if not path.is_dir():
            raise NotADirectoryError(
                f"Expected a sample directory, got: {path}. "
                f"report-multi consumes {args.input_format} output directories."
            )
        name = path.name
        reader = build_reader(args.input_format, path, sample_name=name)
        report = Report.from_experiment(reader.to_summarized_experiment())
        tables.append(_build_table(report, name, index))

    template = env.get_template(_REPORT_TEMPLATE)
    html = template.render(
        page_title=args.title or "Multi-Sample Report",
        tables=tables,
        report_date=datetime.now().strftime("%Y-%m-%d"),
    )

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(html, encoding="utf-8")
        print(f"HTML report: {output_path}")
    else:
        print(html)

    # Optionally also render a PDF via WeasyPrint (mirrors report-single /
    # report-medical). Needs a concrete HTML output path to anchor the PDF path.
    if args.pdf or args.output_pdf:
        if not args.output:
            print(
                "Error: --pdf requires -o/--output so the PDF path can be derived.",
                file=sys.stderr,
            )
            return 1
        pdf_path = _default_pdf_path(Path(args.output), args.output_pdf)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            _generate_pdf(html, pdf_path)
        except RuntimeError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        print(f"PDF report:  {pdf_path}")

    return 0


def _generate_pdf(html_content: str, output_path: Path) -> None:
    """Convert HTML to PDF using WeasyPrint.

    Note: Plotly Sankey (JavaScript) cannot render in WeasyPrint, so the Sankey
    panels show as placeholder text in the PDF. Use the HTML output for the
    interactive diagrams.
    """
    try:
        from weasyprint import HTML  # noqa: PLC0415
    except ImportError as exc:
        raise RuntimeError(
            "WeasyPrint is required for PDF generation. "
            "Install it with: pip install weasyprint"
        ) from exc

    HTML(string=html_content).write_pdf(str(output_path))


def _default_pdf_path(html_path: Path, pdf_arg: str | None) -> Path:
    """Resolve the PDF output path: explicit arg, else HTML path with .pdf suffix."""
    if pdf_arg:
        return Path(pdf_arg)
    return html_path.with_suffix(".pdf")


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "report-multi",
        help="Generate a combined HTML report for multiple samples.",
        description=(
            "Generate a combined multi-sample HTML report from classifier output "
            "directories. Defaults to Savont input; pass --input-format emu for "
            "EMU output folders. Each sample becomes a tab with its species table "
            "and taxonomy Sankey."
        ),
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help=(
            "Sample output directories to include. For --input-format savont: "
            "Savont output directories. For --input-format emu: EMU output "
            "directories (each containing *_rel-abundance.tsv)."
        ),
    )
    parser.add_argument(
        "--input-format",
        choices=SUPPORTED_INPUT_FORMATS,
        default=DEFAULT_INPUT_FORMAT,
        help=f"Classifier input format (default: {DEFAULT_INPUT_FORMAT}).",
    )
    parser.add_argument("-o", "--output", help="Output HTML file path.")
    parser.add_argument("--title", help="Title for the summary report.")
    parser.add_argument(
        "--pdf",
        action="store_true",
        help="Also render a PDF (via WeasyPrint) next to the HTML output. "
        "Requires -o/--output.",
    )
    parser.add_argument(
        "--output-pdf",
        help=(
            "Path to write the PDF report. Implies --pdf. Defaults to the HTML "
            "output path with a .pdf suffix. Requires -o/--output."
        ),
    )

    parser.set_defaults(func=run)
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    build_parser(subparsers)
    raise SystemExit(run(parser.parse_args()))
