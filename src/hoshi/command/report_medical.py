"""Generate a 16S rRNA medical detection report as HTML — end to end.

Consumes a single-sample classifier output (EMU ``*_rel-abundance.tsv`` file or
a Savont output directory) plus an optional clinical metadata JSON, builds a
:class:`Report` intermediate (organisms derived from the experiment; clinical
fields carried in ``Report.metadata``), and renders a clinical-style HTML
report via Jinja2 in a single step.

    inputs -> _load_experiment -> build_medical_report -> Report
           -> report_to_medical_data -> Jinja2 -> HTML

The metadata JSON is optional; any omitted clinical field falls back to
``"N/A"``.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.medical import build_medical_report, report_to_medical_data
from hoshi.lib.report import Report

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_REPORT_TEMPLATE = "medical/medical_report.html.j2"

_SUPPORTED_INPUT_FORMATS = ("emu", "savont")

# Clinical fields sourced from the metadata JSON; missing ones fall back to N/A.
_NA = "N/A"
_CLINICAL_KEYS = (
    "report_id",
    "patient_id",
    "specimen_id",
    "specimen_type",
    "collection_date",
    "status",
    "conclusion",
    "lab_name",
    "reference_db",
    "method",
    "authorized_by",
)


def generate_medical_report(data: dict, report_date: str | None = None) -> str:
    """Render the medical report HTML from a data dict.

    Parameters
    ----------
    data : dict
        Report data containing patient info, organisms, QC, etc.
    report_date : str, optional
        Override the report date. Defaults to current date/time.
    """
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2"), default_for_string=True),
    )
    template = env.get_template(_REPORT_TEMPLATE)

    if report_date is None:
        report_date = datetime.now().strftime("%d %b %Y %H:%M")

    # Build QC items list; default to a pass/pass pair when none supplied.
    qc_items = data.get("qc_items")
    if qc_items is None:
        qc_items = [
            {"name": "Read quality", "status": "pass"},
            {"name": "Read support", "status": "pass"},
        ]

    return template.render(
        page_title="Full-Length 16S rRNA Bacterial Detection Report",
        report_id=data["report_id"],
        report_date=report_date,
        status=data.get("status", "Final"),
        lab_name=data.get("lab_name"),
        patient_id=data["patient_id"],
        specimen_id=data["specimen_id"],
        specimen_type=data["specimen_type"],
        collection_date=data["collection_date"],
        conclusion=data.get("conclusion", "not_detected"),
        organisms=data.get("organisms", []),
        qc_items=qc_items,
        reference_db=data.get("reference_db"),
        method=data.get("method"),
        authorized_by=data.get("authorized_by"),
    )


def generate_medical_report_from_report(
    report: Report,
    *,
    report_date: str | None = None,
    top: int | None = None,
) -> str:
    """Render the medical report HTML from a :class:`Report` intermediate.

    This is the ``inputs -> Report -> Jinja2`` path: the ``Report`` is flattened
    into the render dict (organisms derived from the experiment, clinical fields
    and pathogenicity combined from ``metadata``), then handed to
    :func:`generate_medical_report`.
    """
    data = report_to_medical_data(report, top=top)
    return generate_medical_report(data, report_date=report_date)


def _load_experiment(input_format: str, input_path: Path) -> SummarizedExperiment:
    """Load a single-sample classifier output into a SummarizedExperiment."""
    # Imported lazily to keep the heavy ingress deps off the JSON-only path.
    from hoshi.lib.ingress import (  # noqa: PLC0415
        read_emu_abundance_into_summarizedexperiment,
        read_savont_abundance_into_summarizedexperiment,
    )

    if input_format == "emu":
        if not input_path.is_file():
            raise ValueError(f"EMU input must be a TSV file: {input_path}")
        return read_emu_abundance_into_summarizedexperiment(input_path)

    if input_format == "savont":
        if not input_path.is_dir():
            raise ValueError(f"Savont input must be a directory: {input_path}")
        return read_savont_abundance_into_summarizedexperiment(
            input_path, sample_names=[input_path.name]
        )

    raise ValueError(
        f"Unsupported input format '{input_format}'. "
        f"Supported: {', '.join(_SUPPORTED_INPUT_FORMATS)}"
    )


def _extract_confidence(experiment: SummarizedExperiment) -> dict[str, float]:
    """Pull per-``tax_id`` confidence for the single sample, if the source has it.

    Savont stores ``metadata["species_confidence"][sample] = {tax_id: pct}``.
    Sources without a per-call identity signal (e.g. EMU) yield ``{}``.
    """
    per_sample = experiment.metadata.get("species_confidence")
    if not per_sample:
        return {}
    if experiment.n_samples != 1:
        return {}
    sample = str(experiment.sample_ids[0])
    return dict(per_sample.get(sample, {}))


def _load_metadata(metadata_path: Path | None) -> dict:
    """Load the optional clinical metadata JSON.

    Returns an empty dict when no path is given. Only recognized clinical keys
    plus ``qc_items`` are carried through; everything else is ignored.
    """
    if metadata_path is None:
        return {}
    if not metadata_path.is_file():
        raise ValueError(f"Metadata file not found: {metadata_path}")
    with open(metadata_path, encoding="utf-8") as f:
        raw = json.load(f)
    if not isinstance(raw, dict):
        raise ValueError("Metadata JSON must be a JSON object.")
    return raw


def _build_clinical(meta: dict) -> dict:
    """Build the clinical envelope from metadata, defaulting missing keys to N/A."""
    clinical: dict = {}
    for key in _CLINICAL_KEYS:
        value = meta.get(key)
        clinical[key] = value if value is not None else _NA
    return clinical


def _determine_conclusion(organisms: list[dict]) -> str:
    """Default conclusion when the metadata JSON does not specify one."""
    return "not_detected" if not organisms else "organism_detected"


def _generate_pdf(html_content: str, output_path: Path) -> None:
    """Convert HTML to PDF using WeasyPrint.

    The medical template carries print/``@page`` CSS (margins, non-splitting
    tables, repeating footer) that WeasyPrint applies during conversion.
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
    # Classifier output + optional metadata JSON -> Report -> HTML.
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input not found: {input_path}", file=sys.stderr)
        return 1

    try:
        experiment = _load_experiment(args.input_format, input_path)
        meta = _load_metadata(Path(args.metadata) if args.metadata else None)
    except (ValueError, KeyError, FileNotFoundError, json.JSONDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    confidence = _extract_confidence(experiment)
    clinical = _build_clinical(meta)
    qc_items = meta.get("qc_items")

    report = build_medical_report(
        experiment,
        clinical=clinical,
        qc_items=qc_items,
        confidence=confidence,
    )

    # Conclusion comes from the metadata JSON; auto-derive only when absent.
    data = report_to_medical_data(report, top=args.top)
    if meta.get("conclusion") is None:
        report = report.with_metadata(conclusion=_determine_conclusion(data["organisms"]))

    try:
        html = generate_medical_report_from_report(report, top=args.top)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"HTML report: {output_path}")

    # Optionally also render a PDF via WeasyPrint.
    if args.pdf or args.output_pdf:
        pdf_path = _default_pdf_path(output_path, args.output_pdf)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            _generate_pdf(html, pdf_path)
        except RuntimeError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        print(f"PDF report:  {pdf_path}")

    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "report-medical",
        help="Generate a 16S rRNA medical detection report as HTML (end to end).",
        description=(
            "Generate a clinical-style HTML report for 16S rRNA bacterial detection. "
            "Consumes a single-sample classifier output (EMU *_rel-abundance.tsv file "
            "or a Savont output directory) plus an optional clinical metadata JSON, "
            "builds a Report intermediate, and renders HTML in one step. Missing "
            "clinical fields fall back to 'N/A'."
        ),
    )
    parser.add_argument(
        "input",
        help=(
            "Classifier output path. For --input-format emu: an EMU "
            "*_rel-abundance.tsv file. For --input-format savont: a Savont output "
            "directory."
        ),
    )
    parser.add_argument(
        "--input-format",
        choices=_SUPPORTED_INPUT_FORMATS,
        default="emu",
        help="Classifier input format (default: emu).",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        help=(
            "Path to an optional clinical metadata JSON (report_id, patient_id, "
            "specimen_*, conclusion, authorized_by, qc_items, ...). Omitted fields "
            "fall back to 'N/A'."
        ),
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="Number of top organisms to include (default: 5).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="report_medical.html",
        help="Path to write the generated HTML report (default: report_medical.html).",
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
