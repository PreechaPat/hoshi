"""Tests for the medical Report intermediate and its Jinja2 rendering path.

Covers the flow: inputs -> Report (intermediate) -> flatten -> HTML.
"""

from __future__ import annotations

import pandas as pd
import pytest

from hoshi.command.report_medical import generate_medical_report_from_report
from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.medical import (
    build_medical_report,
    report_to_medical_data,
)
from hoshi.lib.report import Report


def _make_experiment() -> SummarizedExperiment:
    """A 3-feature (+1 meta row) single-sample experiment."""
    abundance = pd.DataFrame(
        {"s1": [0.6, 0.3, 0.1, 0.05]},
        index=["1496", "562", "1351", "unmapped"],
    )
    row_data = pd.DataFrame(
        {"species": ["Clostridioides difficile", "Escherichia coli", "Enterococcus faecalis", ""]},
        index=["1496", "562", "1351", "unmapped"],
    )
    return SummarizedExperiment(
        assays={"abundance": abundance},
        row_data=row_data,
        col_data=pd.DataFrame(index=["s1"]),
    )


_CLINICAL = {
    "report_id": "16S-2026-000184",
    "patient_id": "HN-123456",
    "specimen_id": "SP-26-001842",
    "specimen_type": "Synovial fluid",
    "collection_date": "23 Aug 2026",
    "conclusion": "pathogen_detected",
}


# ─── build_medical_report ────────────────────────────────────────────


def test_build_medical_report_stores_everything_in_metadata():
    report = build_medical_report(
        _make_experiment(),
        clinical=_CLINICAL,
        pathogens={"1496": True, "562": False},
        qc_items=[{"name": "Read quality", "status": "pass"}],
        confidence={"1496": 99.8},
    )
    assert isinstance(report, Report)
    # Clinical envelope lives in metadata (Option A).
    assert report.metadata["report_id"] == "16S-2026-000184"
    assert report.metadata["pathogens"] == {"1496": True, "562": False}
    assert report.metadata["qc_items"] == [{"name": "Read quality", "status": "pass"}]
    # Confidence uses the typed field.
    assert report.confidence == {"1496": 99.8}
    # Experiment is held unmodified.
    assert report.experiment.n_features == 4


# ─── report_to_medical_data (flatten) ────────────────────────────────


def test_flatten_derives_organisms_sorted_and_drops_meta_rows():
    report = build_medical_report(_make_experiment(), clinical=_CLINICAL)
    data = report_to_medical_data(report)

    names = [o["name"] for o in data["organisms"]]
    # 'unmapped' meta row dropped; sorted by abundance desc.
    assert names == ["Clostridioides difficile", "Escherichia coli", "Enterococcus faecalis"]
    assert data["organisms"][0]["abundance"] == pytest.approx(60.0)


def test_flatten_combines_pathogens_and_confidence_at_render_time():
    report = build_medical_report(
        _make_experiment(),
        clinical=_CLINICAL,
        pathogens={"1496": True, "562": False},
        confidence={"1496": 99.8},
    )
    data = report_to_medical_data(report)
    by_name = {o["name"]: o for o in data["organisms"]}

    assert by_name["Clostridioides difficile"]["pathogenic"] is True
    assert by_name["Clostridioides difficile"]["identity"] == pytest.approx(99.8)
    assert by_name["Escherichia coli"]["pathogenic"] is False
    # No confidence for this tax_id -> identity stays None.
    assert by_name["Escherichia coli"]["identity"] is None
    # No pathogen entry for this tax_id -> pathogenic stays None.
    assert by_name["Enterococcus faecalis"]["pathogenic"] is None


def test_flatten_top_limits_organisms():
    report = build_medical_report(_make_experiment(), clinical=_CLINICAL)
    data = report_to_medical_data(report, top=2)
    assert len(data["organisms"]) == 2


def test_flatten_requires_clinical_fields():
    report = Report(experiment=_make_experiment(), metadata={"report_id": "x"})
    with pytest.raises(ValueError, match="Missing required clinical fields"):
        report_to_medical_data(report)


# ─── end-to-end render ───────────────────────────────────────────────


def test_render_from_report_produces_html_with_report_values():
    report = build_medical_report(
        _make_experiment(),
        clinical=_CLINICAL,
        pathogens={"1496": True},
        confidence={"1496": 99.8},
        qc_items=[{"name": "Read quality", "status": "pass"}],
    )
    html = generate_medical_report_from_report(report, report_date="31 Aug 2026")

    assert "16S-2026-000184" in html
    assert "PATHOGEN DETECTED" in html
    assert "Clostridioides difficile" in html
    assert "99.8" in html
