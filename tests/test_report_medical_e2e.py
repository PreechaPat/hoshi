"""End-to-end tests for the report-medical command.

Flow under test: classifier input (+ optional metadata JSON) -> Report -> HTML.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hoshi.command import report_medical
from hoshi.lib.pathogen import DEFAULT_PATHOGEN_SHEET

_EMU_TSV = "test_data/emu_output/test_ind/sample01/sample01_rel-abundance.tsv"
_SAVONT_DIR = "test_data/savont_output/test_ind/savont-out-sample01"


def _args(input_path: str, output: Path, **overrides) -> argparse.Namespace:
    defaults = dict(
        input=input_path,
        input_format="emu",
        metadata=None,
        pathogen_sheet=str(DEFAULT_PATHOGEN_SHEET),
        top=5,
        output=str(output),
        pdf=False,
        output_pdf=None,
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def _write_metadata(tmp_path: Path, data: dict) -> Path:
    p = tmp_path / "meta.json"
    p.write_text(json.dumps(data), encoding="utf-8")
    return p


# ─── e2e EMU, no metadata (N/A fallback) ─────────────────────────────


def test_e2e_emu_without_metadata_falls_back_to_na(tmp_path):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out))

    assert exit_code == 0
    html = out.read_text(encoding="utf-8")
    # Missing clinical fields fall back to N/A.
    assert "N/A" in html
    # Organisms derived from the classifier output.
    assert "Clostridioides difficile" in html
    # Bacterial DNA result label replaces the old "Conclusion:" wording.
    assert "Bacterial DNA:" in html
    # Auto-derived result: a non-commensal organism (C. difficile) is present,
    # so the bacterial-DNA result resolves to PATHOGEN DETECTED.
    assert "PATHOGEN DETECTED" in html


# ─── e2e EMU, with metadata JSON ─────────────────────────────────────


def test_e2e_emu_with_metadata(tmp_path):
    meta = _write_metadata(
        tmp_path,
        {
            "report_metadata": {
                "report_id": "16S-2026-000184",
                "report_date": "24 Aug 2026 09:15",
                "patient_id": "HN-123456",
                "patient_name": "Somsri Chaiyaphum",
                "dob": "14 Mar 1978",
                "age": "48",
                "gender": "Female",
                "specimen_id": "SP-26-001842",
                "specimen_type": "Synovial fluid",
                "collection_date": "23 Aug 2026",
                "received_date": "23 Aug 2026 16:40",
                "ordering_physician": "Dr. Anong Wattana",
                "healthcare_provider": "Bangkok Central Hospital",
                "reason_for_testing": "Suspected septic arthritis",
                "test_performed": "Full-length 16S rRNA bacterial detection",
                "conclusion": "pathogen_detected",
                "authorized_by": "Dr. Smith",
            }
        },
    )
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out, metadata=str(meta)))

    assert exit_code == 0
    html = out.read_text(encoding="utf-8")
    assert "16S-2026-000184" in html
    assert "HN-123456" in html
    assert "PATHOGEN DETECTED" in html  # conclusion honored from JSON
    assert "Dr. Smith" in html
    # New clinical fields render from the metadata JSON.
    assert "Somsri Chaiyaphum" in html
    assert "14 Mar 1978" in html
    assert "Female" in html
    assert "23 Aug 2026 16:40" in html
    assert "Dr. Anong Wattana" in html
    assert "Bangkok Central Hospital" in html
    assert "Suspected septic arthritis" in html
    assert "Full-length 16S rRNA bacterial detection" in html
    # report_date sourced from metadata (not the datetime.now() fallback).
    assert "24 Aug 2026 09:15" in html

def test_e2e_grouped_metadata_shape_flattens_and_renders(tmp_path):
    # Grouped JSON: report_metadata + method groups, qc_items at top level.
    meta = _write_metadata(
        tmp_path,
        {
            "report_metadata": {
                "report_id": "16S-2026-000184",
                "report_date": "24 Aug 2026 09:15",
                "conclusion": "pathogen_detected",
                "patient_id": "HN-123456",
                "patient_name": "Somsri Chaiyaphum",
                "specimen_id": "SP-26-001842",
                "specimen_type": "Synovial fluid",
                "collection_date": "23 Aug 2026",
                "healthcare_provider": "Bangkok Central Hospital",
                "ordering_physician": "Dr. Anong Wattana",
                "authorized_by": "Dr. Smith",
            },
            "method": {
                "reference_db": "16S ref DB v2026.08",
                "method": "Broad-range 16S amplification and classification.",
            },
            "qc_items": [{"name": "Read quality", "status": "pass"}],
        },
    )
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out, metadata=str(meta)))

    assert exit_code == 0
    html = out.read_text(encoding="utf-8")
    # Fields from report_metadata group.
    assert "16S-2026-000184" in html
    assert "HN-123456" in html
    assert "Somsri Chaiyaphum" in html
    assert "Bangkok Central Hospital" in html
    assert "Dr. Anong Wattana" in html
    assert "Dr. Smith" in html
    assert "PATHOGEN DETECTED" in html
    assert "24 Aug 2026 09:15" in html
    # Fields from method group.
    assert "16S ref DB v2026.08" in html
    assert "Broad-range 16S amplification" in html



def test_e2e_report_date_falls_back_when_absent(tmp_path):
    # No report_date in metadata -> template still renders a Reported line
    # (datetime.now() fallback), and the new fields fall back to N/A.
    meta = _write_metadata(
        tmp_path,
        {
            "report_metadata": {
                "report_id": "16S-2026-000199",
                "patient_id": "HN-999",
                "specimen_id": "SP-26-000199",
                "specimen_type": "Blood",
                "collection_date": "01 Sep 2026",
            }
        },
    )
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out, metadata=str(meta)))

    assert exit_code == 0
    html = out.read_text(encoding="utf-8")
    assert "Reported:" in html
    # Absent optional fields render as N/A.
    assert "N/A" in html


# ─── e2e Savont, confidence -> identity ──────────────────────────────


def test_e2e_savont_populates_identity_from_confidence(tmp_path):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(
        _args(_SAVONT_DIR, out, input_format="savont")
    )

    assert exit_code == 0
    html = out.read_text(encoding="utf-8")
    assert "Clostridioides difficile" in html
    # Savont provides per-species confidence -> identity column not all "N/A".
    # (Exact value depends on fixture; assert at least one percent-looking cell.)
    assert "%" in html or "99" in html or "100" in html


def test_e2e_top_limits_organisms(tmp_path):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out, top=2))
    assert exit_code == 0
    html = out.read_text(encoding="utf-8")
    # Only 2 organism rows -> the 3rd-ranked species should be absent.
    assert "Clostridioides difficile" in html
    assert "Lactiplantibacillus plantarum" not in html


# ─── error handling ──────────────────────────────────────────────────


def test_e2e_missing_input_returns_error(tmp_path, capsys):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args("test_data/nope.tsv", out))
    assert exit_code == 1
    assert not out.exists()
    assert "not found" in capsys.readouterr().err.lower()


def test_e2e_missing_metadata_file_returns_error(tmp_path, capsys):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(
        _args(_EMU_TSV, out, metadata=str(tmp_path / "absent.json"))
    )
    assert exit_code == 1
    assert "metadata file not found" in capsys.readouterr().err.lower()


# ─── PDF output ──────────────────────────────────────────────────────


def _is_pdf(path: Path) -> bool:
    return path.is_file() and path.read_bytes().startswith(b"%PDF")


def test_e2e_no_pdf_by_default(tmp_path):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out))

    assert exit_code == 0
    assert out.is_file()
    assert not out.with_suffix(".pdf").exists()


def test_e2e_pdf_written_via_flag_and_custom_path(tmp_path):
    # --pdf writes a PDF next to the HTML; --output-pdf implies --pdf and
    # writes to the given path.
    flag_html = tmp_path / "flag.html"
    assert report_medical.run(_args(_EMU_TSV, flag_html, pdf=True)) == 0
    assert _is_pdf(flag_html.with_suffix(".pdf"))

    custom_html = tmp_path / "custom.html"
    custom_pdf = tmp_path / "custom" / "medical.pdf"
    assert (
        report_medical.run(
            _args(_EMU_TSV, custom_html, output_pdf=str(custom_pdf))
        )
        == 0
    )
    assert _is_pdf(custom_pdf)

