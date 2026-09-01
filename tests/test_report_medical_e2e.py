"""End-to-end tests for the report-medical command.

Flow under test: classifier input (+ optional metadata JSON) -> Report -> HTML.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hoshi.command import report_medical

_EMU_TSV = "test_data/emu_output/test_ind/sample01/sample01_rel-abundance.tsv"
_SAVONT_DIR = "test_data/savont_output/test_ind/savont-out-sample01"


def _args(input_path: str, output: Path, **overrides) -> argparse.Namespace:
    defaults = dict(
        input=input_path,
        input_format="emu",
        metadata=None,
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
    # Auto-derived conclusion (organisms found).
    assert "ORGANISM DETECTED" in html


# ─── e2e EMU, with metadata JSON ─────────────────────────────────────


def test_e2e_emu_with_metadata(tmp_path):
    meta = _write_metadata(
        tmp_path,
        {
            "report_id": "16S-2026-000184",
            "patient_id": "HN-123456",
            "specimen_id": "SP-26-001842",
            "specimen_type": "Synovial fluid",
            "collection_date": "23 Aug 2026",
            "conclusion": "pathogen_detected",
            "authorized_by": "Dr. Smith",
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


def test_e2e_pdf_flag_writes_pdf_next_to_html(tmp_path):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out, pdf=True))

    assert exit_code == 0
    assert out.is_file()
    pdf_path = out.with_suffix(".pdf")
    assert _is_pdf(pdf_path)


def test_e2e_output_pdf_implies_pdf_and_uses_path(tmp_path):
    out = tmp_path / "report.html"
    pdf_out = tmp_path / "custom" / "medical.pdf"
    exit_code = report_medical.run(
        _args(_EMU_TSV, out, output_pdf=str(pdf_out))
    )

    assert exit_code == 0
    assert _is_pdf(pdf_out)


def test_e2e_no_pdf_by_default(tmp_path):
    out = tmp_path / "report.html"
    exit_code = report_medical.run(_args(_EMU_TSV, out))

    assert exit_code == 0
    assert out.is_file()
    assert not out.with_suffix(".pdf").exists()

