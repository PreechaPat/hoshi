"""Generate a single-sample HTML report, optionally with a paged PDF.

The HTML output is responsive (single-page, navbar, mobile-friendly).
The PDF output is a formal 5-page A4 report (Cover, Raw Sequences, Species,
Taxonomy Sankey, Appendix/Citation/Disclaimer) generated via WeasyPrint.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.diversity import compute_diversity
from hoshi.lib.reader import (
    DEFAULT_INPUT_FORMAT,
    SUPPORTED_INPUT_FORMATS,
    build_reader,
)
from hoshi.lib.report import Report
from hoshi.lib.sankey import get_sankey_data, render_sankey_figure

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_REPORT_TEMPLATE = "microbiome/singlesample.html.j2"


def _prepare_report_data(
    report: Report,
    sample_name: str,
) -> tuple[dict, object]:
    """Prepare shared report data from a :class:`Report`.

    Returns (report_data dict, DiversityStats). Confidence carried on the
    ``Report`` is surfaced as a headline percentage and per-species column.
    """
    se = report.experiment
    df = se.to_dataframe()
    stats = compute_diversity(df)

    # Filter out meta rows for display
    meta_ids = {"unmapped", "mapped_filtered", "mapped_unclassified"}
    df_display = df[~df["tax_id"].astype(str).isin(meta_ids)].copy()
    df_display = df_display.sort_values("abundance", ascending=False).reset_index(drop=True)

    # Surface per-species confidence (Savont) as an extra column when present.
    confidence = report.confidence
    if confidence:
        df_display["confidence"] = (
            df_display["tax_id"].astype(str).map(confidence).round(1)
        )
        display_cols = ["species", "tax_id", "abundance", "estimated counts", "confidence"]
    else:
        display_cols = ["species", "tax_id", "abundance", "estimated counts"]
    df_subset = df_display[display_cols]

    top5_html = df_subset.head(5).to_html(index=False, border=0, classes="data-table")
    top10_html = df_subset.head(10).to_html(index=False, border=0, classes="data-table")

    # Sankey diagram (interactive for HTML)
    sankey_html = ""
    sankey_svg = ""
    try:
        sankey_data = get_sankey_data(df)
        fig = render_sankey_figure(sankey_data, title=f"Taxonomy Flow: {sample_name}")
        sankey_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
        # Static SVG for PDF (kaleido required)
        try:
            sankey_svg = fig.to_image(format="svg", width=800, height=600).decode("utf-8")
        except Exception:
            sankey_svg = ""
    except Exception as e:
        sankey_html = f"<p>Sankey plot could not be generated: {e!s}</p>"

    # Source file info
    source_file = ""
    sample_dir = ""
    if not se.col_data.empty and "source_file" in se.col_data.columns:
        source_file = se.col_data.iloc[0]["source_file"]
        sample_dir = str(Path(source_file).parent.resolve())

    confidence_pct = report.confidence_pct

    report_data = {
        "name": sample_name,
        "source_path": source_file,
        "sample_dir": sample_dir,
        "source": se.metadata.get("source", ""),
        "confidence_pct": round(confidence_pct, 1) if confidence_pct is not None else None,
        "top5_html": top5_html,
        "top10_html": top10_html,
        "sankey_html": sankey_html,
        "sankey_svg": sankey_svg,
    }

    return report_data, stats


def generate_single_html_report(
    report: Report,
    *,
    page_title: str | None = None,
    sample_name: str | None = None,
) -> str:
    """Generate a report from the unified template.

    On screen: renders as a sidebar-navigated single-page app.
    In print/PDF: renders as paged A4 via CSS @media print rules.

    Parameters
    ----------
    report : Report
        A single-sample report composed over a SummarizedExperiment.
    page_title : str, optional
        Override the full page title.
    sample_name : str, optional
        Sample name for the report header. Falls back to the experiment's
        sample_id.
    """
    name = sample_name or str(report.experiment.sample_ids[0])
    title = page_title or name

    report_data, stats = _prepare_report_data(report, name)

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2"), default_for_string=True),
    )
    template = env.get_template(_REPORT_TEMPLATE)
    return template.render(
        page_title=title,
        table=report_data,
        stats=stats,
        report_date=datetime.now().strftime("%Y-%m-%d"),
    )


def generate_paged_html_report(
    report: Report,
    *,
    page_title: str | None = None,
    sample_name: str | None = None,
) -> str:
    """Generate HTML for PDF conversion (same template, WeasyPrint uses @media print).

    The unified template renders as paged A4 when processed by WeasyPrint
    because WeasyPrint applies @media print and @page CSS rules.
    """
    name = sample_name or str(report.experiment.sample_ids[0])
    title = page_title or name

    report_data, stats = _prepare_report_data(report, name)

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2"), default_for_string=True),
    )
    template = env.get_template(_REPORT_TEMPLATE)
    return template.render(
        page_title=title,
        table=report_data,
        stats=stats,
        report_date=datetime.now().strftime("%Y-%m-%d"),
    )


def _generate_pdf(html_content: str, output_path: Path) -> None:
    """Convert HTML to PDF using WeasyPrint.

    Note: Plotly Sankey (JavaScript) cannot render in WeasyPrint.
    The Sankey page will show as placeholder text in the PDF.
    For interactive diagrams, use the HTML output.
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


def run(args: argparse.Namespace) -> int:
    sample_dir = Path(args.sample_dir)

    if not sample_dir.is_dir():
        raise NotADirectoryError(f"Expected a sample directory, got: {sample_dir}")

    # Pipeline: reader (savont default / emu) → SummarizedExperiment → Report
    # → reports. The reader owns its own folder layout and file discovery; the
    # Report composite layers report-time data (confidence) on the experiment.
    name = args.name or sample_dir.name
    reader = build_reader(args.input_format, sample_dir, sample_name=name)
    se = reader.to_summarized_experiment()
    report = Report.from_experiment(se)

    # Generate responsive HTML report
    html = generate_single_html_report(report, page_title=args.title, sample_name=name)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"HTML report: {output_path}")

    # Optionally also render a PDF via WeasyPrint.
    if args.pdf or args.output_pdf:
        paged_html = generate_paged_html_report(
            report, page_title=args.title, sample_name=name
        )
        pdf_path = _default_pdf_path(output_path, args.output_pdf)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            _generate_pdf(paged_html, pdf_path)
        except RuntimeError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        print(f"PDF report:  {pdf_path}")

    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "report-single",
        help="Generate a single-sample HTML report from a classifier output folder.",
        description=(
            "Generate an HTML microbiome report from a classifier output directory. "
            "Defaults to Savont input; pass --input-format emu for an EMU output "
            "folder (the reader discovers *_rel-abundance.tsv and other files). "
            "Renders a responsive HTML report; pass --pdf for a formal A4 PDF too."
        ),
    )
    parser.add_argument(
        "sample_dir",
        help=(
            "Path to the sample output directory. For --input-format savont: a "
            "Savont output directory. For --input-format emu: an EMU output "
            "directory containing *_rel-abundance.tsv."
        ),
    )
    parser.add_argument(
        "--input-format",
        choices=SUPPORTED_INPUT_FORMATS,
        default=DEFAULT_INPUT_FORMAT,
        help=f"Classifier input format (default: {DEFAULT_INPUT_FORMAT}).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="report_single.html",
        help="Path to write the generated HTML report (default: report_single.html).",
    )
    parser.add_argument(
        "-n",
        "--name",
        help="Sample name displayed in the report header (e.g. 'sample01'). "
        "Defaults to the folder name.",
    )
    parser.add_argument(
        "--title",
        help="Override the full page title (takes precedence over --name).",
    )
    parser.add_argument(
        "--pdf",
        action="store_true",
        help="Also render a PDF (via WeasyPrint) next to the HTML output.",
    )
    parser.add_argument(
        "--output-pdf",
        help=(
            "Path to write the PDF report. Implies --pdf. Defaults to the HTML "
            "output path with a .pdf suffix."
        ),
    )

    parser.set_defaults(func=run)
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    build_parser(subparsers)
    raise SystemExit(run(parser.parse_args()))
