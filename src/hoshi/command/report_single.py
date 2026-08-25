"""Generate a single-sample HTML report, optionally with a paged PDF.

The HTML output is responsive (single-page, navbar, mobile-friendly).
The PDF output is a formal 5-page A4 report (Cover, Raw Sequences, Species,
Taxonomy Sankey, Appendix/Citation/Disclaimer) generated via WeasyPrint.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.diversity import compute_diversity
from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.ingress import emu_to_experiment
from hoshi.lib.sankey import get_sankey_data, render_sankey_figure

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_REPORT_TEMPLATE = "single_report.html.j2"


def _find_emu_files(sample_dir: Path) -> dict[str, Path | None]:
    """Discover EMU output files in a sample directory.

    Returns a dict with keys: 'abundance', 'read_assignments', 'unclassified', 'unmapped'.
    Values are Path objects or None if not found.
    """
    files: dict[str, Path | None] = {
        "abundance": None,
        "read_assignments": None,
        "unclassified": None,
        "unmapped": None,
    }

    for f in sample_dir.iterdir():
        if not f.is_file():
            continue
        name = f.name
        if name.endswith("_rel-abundance.tsv"):
            files["abundance"] = f
        elif name.endswith("_read-assignment-distributions.tsv"):
            files["read_assignments"] = f
        elif name.endswith("_unclassified_mapped.fastq.gz"):
            files["unclassified"] = f
        elif name.endswith("_unmapped.fastq.gz"):
            files["unmapped"] = f

    return files


def _prepare_report_data(
    se: SummarizedExperiment,
    sample_name: str,
) -> tuple[dict, object]:
    """Prepare shared report data from a SummarizedExperiment.

    Returns (report_data dict, DiversityStats).
    """
    df = se.to_dataframe()
    stats = compute_diversity(df)

    # Filter out meta rows for display
    meta_ids = {"unmapped", "mapped_filtered", "mapped_unclassified"}
    df_display = df[~df["tax_id"].astype(str).isin(meta_ids)].copy()
    df_display = df_display.sort_values("abundance", ascending=False).reset_index(drop=True)
    df_subset = df_display[["species", "tax_id", "abundance", "estimated counts"]]

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

    report_data = {
        "name": sample_name,
        "source_path": source_file,
        "sample_dir": sample_dir,
        "top5_html": top5_html,
        "top10_html": top10_html,
        "sankey_html": sankey_html,
        "sankey_svg": sankey_svg,
    }

    return report_data, stats


def generate_single_html_report(
    se: SummarizedExperiment,
    *,
    page_title: str | None = None,
    sample_name: str | None = None,
) -> str:
    """Generate a report from the unified template.

    On screen: renders as a sidebar-navigated single-page app.
    In print/PDF: renders as paged A4 via CSS @media print rules.

    Parameters
    ----------
    se : SummarizedExperiment
        A single-sample experiment container with assays and row_data.
    page_title : str, optional
        Override the full page title.
    sample_name : str, optional
        Sample name for the report header. Falls back to SE's sample_id.
    """
    name = sample_name or str(se.sample_ids[0])
    title = page_title or name

    report_data, stats = _prepare_report_data(se, name)

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
    se: SummarizedExperiment,
    *,
    page_title: str | None = None,
    sample_name: str | None = None,
) -> str:
    """Generate HTML for PDF conversion (same template, WeasyPrint uses @media print).

    The unified template renders as paged A4 when processed by WeasyPrint
    because WeasyPrint applies @media print and @page CSS rules.
    """
    name = sample_name or str(se.sample_ids[0])
    title = page_title or name

    report_data, stats = _prepare_report_data(se, name)

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


def run(args: argparse.Namespace) -> int:
    sample_dir = Path(args.sample_dir)

    if not sample_dir.is_dir():
        raise NotADirectoryError(f"Expected a sample directory, got: {sample_dir}")

    # Discover Emu files
    emu_files = _find_emu_files(sample_dir)
    abundance_file = emu_files["abundance"]

    if abundance_file is None:
        raise FileNotFoundError(
            f"No *_rel-abundance.tsv file found in {sample_dir}. "
            "Is this an EMU output directory?"
        )

    # Pipeline: Emu TSV → SummarizedExperiment → generate reports
    name = args.name or sample_dir.name
    se = emu_to_experiment(abundance_file, sample_names=[name])

    # Generate responsive HTML report
    html = generate_single_html_report(se, page_title=args.title, sample_name=name)
    html_path = Path(args.output_html)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(html, encoding="utf-8")
    print(f"HTML report: {html_path}")

    # Generate PDF report
    paged_html = generate_paged_html_report(se, page_title=args.title, sample_name=name)
    pdf_path = Path(args.output_pdf)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    _generate_pdf(paged_html, pdf_path)
    print(f"PDF report:  {pdf_path}")

    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "report-single",
        help="Generate a single-sample HTML report from an EMU output folder.",
        description=(
            "Generate an HTML microbiome report from an EMU output directory. "
            "Automatically discovers *_rel-abundance.tsv and other EMU files. "
            "Produces both a responsive HTML report and a formal 5-page A4 PDF."
        ),
    )
    parser.add_argument(
        "sample_dir",
        help="Path to the EMU sample output directory (contains *_rel-abundance.tsv).",
    )
    parser.add_argument(
        "--output-html",
        default="output.html",
        help="Path to write the generated HTML report (default: output.html).",
    )
    parser.add_argument(
        "--output-pdf",
        default="output.pdf",
        help="Path to write the generated PDF report (default: output.pdf).",
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

    parser.set_defaults(func=run)
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    build_parser(subparsers)
    raise SystemExit(run(parser.parse_args()))
