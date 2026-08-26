"""Generate a 16S rRNA medical detection report as HTML.

Accepts a JSON input file describing patient information, specimen details,
organisms detected, and QC status. Produces a clinical-style HTML report
suitable for printing or PDF conversion.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_REPORT_TEMPLATE = "medical/medical_report.html.j2"


def _load_report_data(json_path: Path) -> dict:
    """Load and validate report data from a JSON file.

    Expected JSON structure:
    {
        "report_id": "16S-2026-000184",
        "patient_id": "HN-XXXXXX",
        "specimen_id": "SP-26-001842",
        "specimen_type": "Synovial fluid",
        "collection_date": "23 Aug 2026",
        "conclusion": "pathogen_detected",
        "organisms": [
            {
                "name": "Streptococcus intermedius",
                "abundance": 85.2,
                "identity": 99.8,
                "pathogenic": true
            }
        ],
        "qc_items": [
            {"name": "Read quality", "status": "pass"},
            {"name": "Read support", "status": "pass"}
        ],
        "lab_name": "Molecular Microbiology Laboratory",
        "reference_db": "Validated bacterial 16S reference database, version 2026.08",
        "authorized_by": "Dr. Smith"
    }

    conclusion values:
        "pathogen_detected" - known pathogen found
        "organism_detected" - organism found in normally sterile specimen
        "normal"            - normal flora only
        "not_detected"      - no bacterial DNA detected

    Each organism entry:
        "name"       - species name (required)
        "abundance"  - relative abundance percentage
        "identity"   - sequence identity percentage
        "pathogenic" - boolean

    Each qc_items entry has "name" and "status" ("pass" or "fail").
    """
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)

    # Validate required fields
    required = ["report_id", "patient_id", "specimen_id", "specimen_type", "collection_date"]
    missing = [k for k in required if k not in data]
    if missing:
        raise ValueError(f"Missing required fields in JSON: {', '.join(missing)}")

    # Normalize organisms: ensure it's a list
    if "organisms" not in data:
        data["organisms"] = []

    for org in data["organisms"]:
        if "name" not in org:
            raise ValueError("Each organism must have a 'name' field")

    return data


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

    # Build QC items list
    qc_items = data.get("qc_items")
    if qc_items is None:
        # Backward compat: build from legacy single qc_status field
        qc_items = [
            {"name": "Read quality", "status": data.get("qc_status", "pass")},
            {"name": "Read support", "status": data.get("qc_status", "pass")},
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


def run(args: argparse.Namespace) -> int:
    json_path = Path(args.input)

    if not json_path.is_file():
        print(f"Error: Input file not found: {json_path}", file=sys.stderr)
        return 1

    try:
        data = _load_report_data(json_path)
    except (json.JSONDecodeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    html = generate_medical_report(data)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"Medical report: {output_path}")

    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "report-medical",
        help="Generate a 16S rRNA medical detection report as HTML.",
        description=(
            "Generate a clinical-style HTML report for 16S rRNA bacterial detection. "
            "Input is a JSON file describing patient information, specimen details, "
            "detected organisms, and QC status."
        ),
    )
    parser.add_argument(
        "input",
        help="Path to JSON file with report data (patient info, organisms, etc.).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="report_medical.html",
        help="Path to write the generated HTML report (default: report_medical.html).",
    )

    parser.set_defaults(func=run)
    return parser
