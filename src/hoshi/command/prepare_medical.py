"""Generate a JSON input file for report-medical from EMU abundance output.

Reads the EMU *_rel-abundance.tsv, extracts top organisms, and combines
with clinical metadata from CLI arguments to produce the JSON structure
expected by `hoshi report-medical`.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

_META_TAX_IDS = {"unmapped", "mapped_filtered", "mapped_unclassified"}


def _read_emu_abundance(tsv_path: Path, top_n: int) -> list[dict]:
    """Read EMU abundance TSV and return top N organisms as dicts."""
    df = pd.read_csv(tsv_path, sep="\t")

    # Filter out meta rows
    df = df[~df["tax_id"].astype(str).isin(_META_TAX_IDS)].copy()

    # Sort by abundance descending, take top N
    df = df.sort_values("abundance", ascending=False).head(top_n)

    organisms = []
    for _, row in df.iterrows():
        name = row.get("species", "")
        if not name or pd.isna(name):
            # Fall back to genus if species is empty
            name = row.get("genus", "Unknown")
            if pd.isna(name):
                name = "Unknown"

        organisms.append({
            "name": str(name).strip(),
            "abundance": round(float(row["abundance"]) * 100, 2),
            "identity": None,
            "pathogenic": None,
        })

    return organisms


def _determine_conclusion(organisms: list[dict]) -> str:
    """Determine conclusion based on whether organisms were found.

    Returns a default; user can override via --conclusion flag.
    """
    if not organisms:
        return "not_detected"
    return "organism_detected"


def run(args: argparse.Namespace) -> int:
    tsv_path = Path(args.tsv)

    if not tsv_path.is_file():
        print(f"Error: TSV file not found: {tsv_path}", file=sys.stderr)
        return 1

    organisms = _read_emu_abundance(tsv_path, top_n=args.top)

    conclusion = args.conclusion or _determine_conclusion(organisms)

    # Build QC items
    qc_items = [
        {"name": "Read quality", "status": args.qc_read_quality},
        {"name": "Read support", "status": args.qc_read_support},
    ]

    # Build report data
    report_data = {
        "report_id": args.report_id,
        "patient_id": args.patient_id,
        "specimen_id": args.specimen_id,
        "specimen_type": args.specimen_type,
        "collection_date": args.collection_date,
        "conclusion": conclusion,
        "organisms": organisms,
        "qc_items": qc_items,
    }

    # Optional fields
    if args.lab_name:
        report_data["lab_name"] = args.lab_name
    if args.reference_db:
        report_data["reference_db"] = args.reference_db
    if args.authorized_by:
        report_data["authorized_by"] = args.authorized_by

    # Write JSON
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=2, ensure_ascii=False)
    print(f"Prepared: {output_path}")

    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "prepare-medical",
        help="Generate JSON input for report-medical from EMU abundance TSV.",
        description=(
            "Read an EMU *_rel-abundance.tsv file and clinical metadata to produce "
            "the JSON file expected by `hoshi report-medical`."
        ),
    )
    parser.add_argument(
        "tsv",
        help="Path to EMU *_rel-abundance.tsv file.",
    )

    # Clinical metadata
    meta = parser.add_argument_group("clinical metadata")
    meta.add_argument("--report-id", default="N/A", help="Report ID (e.g. 16S-2026-000184).")
    meta.add_argument("--patient-id", default="N/A", help="Patient ID (e.g. HN-123456).")
    meta.add_argument("--specimen-id", default="N/A", help="Specimen ID (e.g. SP-26-001842).")
    meta.add_argument("--specimen-type", default="N/A", help="Specimen type (e.g. Synovial fluid).")
    meta.add_argument(
        "--collection-date",
        default="N/A",
        help="Specimen collection date (e.g. '23 Aug 2026').",
    )

    # Optional overrides
    opts = parser.add_argument_group("options")
    opts.add_argument("--top", type=int, default=5, help="Number of top organisms to include (default: 5).")
    opts.add_argument(
        "--conclusion",
        choices=["pathogen_detected", "organism_detected", "normal", "not_detected"],
        help="Override auto-determined conclusion.",
    )
    opts.add_argument("--lab-name", help="Laboratory name.")
    opts.add_argument("--reference-db", help="Reference database description.")
    opts.add_argument("--authorized-by", help="Authorizing person name.")

    # QC
    qc = parser.add_argument_group("quality control")
    qc.add_argument("--qc-read-quality", default="pass", choices=["pass", "fail"], help="Read quality QC status (default: pass).")
    qc.add_argument("--qc-read-support", default="pass", choices=["pass", "fail"], help="Read support QC status (default: pass).")

    # Output
    parser.add_argument(
        "-o",
        "--output",
        default="medical_input.json",
        help="Path to write the generated JSON file (default: medical_input.json).",
    )

    parser.set_defaults(func=run)
    return parser
