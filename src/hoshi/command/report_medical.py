"""Generate a 16S rRNA medical detection report as HTML — end to end.

Consumes a single-sample classifier output (EMU ``*_rel-abundance.tsv`` file or
a Savont output directory) plus an optional clinical metadata JSON, builds a
:class:`Report` intermediate (organisms derived from the experiment; clinical
fields carried in ``Report.metadata``), and renders a clinical-style HTML
report via Jinja2 in a single step.

    inputs -> _load_experiment -> build_medical_report -> Report
           -> report_to_medical_data -> Jinja2 -> HTML

The metadata JSON is optional; any omitted clinical field is passed to the
template as ``None`` and the template renders its own fallback (``"N/A"`` for
identifying fields, a blank signature line for ``authorized_by``).
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
from hoshi.lib.pathogen import DEFAULT_PATHOGEN_SHEET
from hoshi.lib.report import Report

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_REPORT_TEMPLATE = "medical/medical_report.html.j2"

_SUPPORTED_INPUT_FORMATS = ("savont", "emu")

# Clinical fields sourced from the metadata JSON. Missing ones are passed to the
# template as ``None``; the template owns all presentation defaults (``N/A`` for
# identifying fields, a blank signature line for ``authorized_by``).
_CLINICAL_KEYS = (
    # Report identity
    "report_id",
    "report_date",
    # Patient
    "patient_id",
    "patient_name",
    "dob",
    "age",
    "gender",
    # Specimen
    "specimen_id",
    "specimen_type",
    "collection_date",
    "received_date",
    # Ordering / provider
    "ordering_physician",
    "healthcare_provider",
    "reason_for_testing",
    "test_performed",
    # Report-level result / method
    "status",
    "conclusion",
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
        report_date = data.get("report_date") or datetime.now().strftime("%d %b %Y %H:%M")

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
        patient_id=data["patient_id"],
        patient_name=data.get("patient_name"),
        dob=data.get("dob"),
        age=data.get("age"),
        gender=data.get("gender"),
        specimen_id=data["specimen_id"],
        specimen_type=data["specimen_type"],
        collection_date=data["collection_date"],
        received_date=data.get("received_date"),
        ordering_physician=data.get("ordering_physician"),
        healthcare_provider=data.get("healthcare_provider"),
        reason_for_testing=data.get("reason_for_testing"),
        test_performed=data.get("test_performed"),
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
            input_path, sample_name=input_path.name
        )

    raise ValueError(
        f"Unsupported input format '{input_format}'. "
        f"Supported: {', '.join(_SUPPORTED_INPUT_FORMATS)}"
    )


# Top-level objects in the metadata JSON whose keys are flattened up to the
# clinical envelope. ``report_metadata`` carries the report/patient/specimen/
# provider fields (incl. authorized_by); ``method`` carries reference_db /
# method. Everything else (notably ``qc_items``) is left at the top level.
_METADATA_GROUPS = ("report_metadata", "method")


def _load_metadata(metadata_path: Path | None) -> dict:
    """Load the optional clinical metadata JSON and flatten it.

    The JSON groups clinical fields under ``report_metadata`` (report identity,
    patient, specimen, provider, authorized_by) and ``method`` (reference_db,
    method), with ``qc_items`` at the top level. This flattens those groups back
    into the flat clinical envelope the rest of the pipeline (and the template)
    consumes. A legacy flat JSON — no group objects — is accepted unchanged.
    Returns an empty dict when no path is given.
    """
    if metadata_path is None:
        return {}
    if not metadata_path.is_file():
        raise ValueError(f"Metadata file not found: {metadata_path}")
    with open(metadata_path, encoding="utf-8") as f:
        raw = json.load(f)
    if not isinstance(raw, dict):
        raise ValueError("Metadata JSON must be a JSON object.")
    return _flatten_metadata(raw)


def _flatten_metadata(raw: dict) -> dict:
    """Flatten the grouped metadata objects into a flat clinical dict.

    Keys inside ``report_metadata`` / ``method`` are lifted to the top level;
    any keys already at the top level (e.g. ``qc_items``, or fields from a
    legacy flat JSON) are preserved. Grouped keys take precedence on conflict.
    """
    flat = {key: value for key, value in raw.items() if key not in _METADATA_GROUPS}
    for group in _METADATA_GROUPS:
        section = raw.get(group)
        if section is None:
            continue
        if not isinstance(section, dict):
            raise ValueError(f"Metadata '{group}' must be a JSON object.")
        flat.update(section)
    return flat


def _load_pathogens(sheet_arg: str | None) -> dict[str, str] | None:
    """Load per-tax_id pathogen classification from the sheet.

    Returns ``None`` when no sheet is configured (disables the column). A
    configured-but-missing file is an error.
    """
    if not sheet_arg:
        return None
    from hoshi.lib.pathogen import PathogenDB  # noqa: PLC0415

    return PathogenDB.from_csv(sheet_arg).describe_lookup()


def _build_clinical(meta: dict) -> dict:
    """Build the clinical envelope from metadata.

    Missing keys are passed through as ``None`` (not the literal ``"N/A"``).
    Presentation defaults are owned entirely by the Jinja2 template, which
    renders ``"N/A"`` for absent identifying fields and a blank signature line
    (``____``) for an absent ``authorized_by``.
    """
    return {key: meta.get(key) for key in _CLINICAL_KEYS}


def _is_commensal(pathogenic: str | None) -> bool:
    """True when an organism's pathogen classification is commensal.

    The ``pathogenic`` value is the describe() string from the pathogen sheet
    (e.g. ``"primary"``, ``"opportunistic"``, ``"commensal (gut, stool)"``).
    Only strings beginning with ``"commensal"`` count as commensal; unknown
    (``None``/empty) classifications are treated as non-commensal.
    """
    return bool(pathogenic) and str(pathogenic).strip().lower().startswith("commensal")


def _determine_conclusion(organisms: list[dict]) -> str:
    """Default conclusion when the metadata JSON does not specify one.

    - No organisms detected            -> ``"not_detected"``
    - Any non-commensal organism found -> ``"pathogen_detected"``
    - All detected organisms commensal -> ``"organism_detected"``
    """
    if not organisms:
        return "not_detected"
    if any(not _is_commensal(org.get("pathogenic")) for org in organisms):
        return "pathogen_detected"
    return "organism_detected"


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

    clinical = _build_clinical(meta)
    qc_items = meta.get("qc_items")

    # Per-tax_id pathogen classification from the sheet (commensal / potential /
    # opportunistic / primary; commensal carries its body sites). Absent tax_ids
    # render as N/A. A missing sheet path simply disables the column.
    try:
        pathogens = _load_pathogens(getattr(args, "pathogen_sheet", None))
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    # Sequence identity is auto-extracted from the experiment inside build_medical_report
    # (Savont populates it; EMU does not), so we don't pull it out here.
    report = build_medical_report(
        experiment,
        clinical=clinical,
        pathogens=pathogens,
        qc_items=qc_items,
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
            "clinical fields fall back to 'N/A' (an omitted authorized_by renders "
            "as a blank signature line)."
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
        default="savont",
        help="Classifier input format (default: savont).",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        help=(
            "Path to an optional clinical metadata JSON. Fields are grouped "
            "under 'report_metadata' (report_id, patient_*, specimen_*, "
            "ordering_physician, healthcare_provider, authorized_by, ...) and "
            "'method' (reference_db, method), with 'qc_items' at the top level. "
            "Omitted fields fall back to 'N/A'."
        ),
    )
    parser.add_argument(
        "--pathogen-sheet",
        default=str(DEFAULT_PATHOGEN_SHEET),
        help=(
            "Path to the pathogen sheet CSV used to classify organisms "
            "(commensal / potential / opportunistic / primary), keyed by NCBI "
            "tax_id. Pass an empty value to disable the classification column "
            "(default: the sheet bundled with hoshi)."
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
