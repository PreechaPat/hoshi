"""Tests for the medical Report intermediate and its Jinja2 rendering path.

Covers the flow: inputs -> Report (intermediate) -> flatten -> HTML.
"""

from __future__ import annotations

import pandas as pd
import pytest

from hoshi.command.report_medical import (
    _determine_conclusion,
    generate_medical_report_from_report,
)
from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.medical import (
    build_medical_report,
    report_to_medical_data,
)
from hoshi.lib.report import Report


def _make_experiment() -> SummarizedExperiment:
    """A 3-OTU (+1 meta row) single-sample per-OTU experiment.

    Feature index is per-OTU; tax_id (incl. the ``unmapped`` meta sentinel)
    lives in row_data.
    """
    otu_ids = ["s1:ASV0", "s1:ASV1", "s1:ASV2", "s1:ASV3"]
    abundance = pd.DataFrame({"s1": [0.6, 0.3, 0.1, 0.05]}, index=otu_ids)
    row_data = pd.DataFrame(
        {
            "tax_id": ["1496", "562", "1351", "unmapped"],
            "species": [
                "Clostridioides difficile",
                "Escherichia coli",
                "Enterococcus faecalis",
                "",
            ],
        },
        index=otu_ids,
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
        confidence={"s1:ASV0": 99.8},
    )
    assert isinstance(report, Report)
    # Clinical envelope lives in metadata (Option A).
    assert report.metadata["report_id"] == "16S-2026-000184"
    assert report.metadata["pathogens"] == {"1496": True, "562": False}
    assert report.metadata["qc_items"] == [{"name": "Read quality", "status": "pass"}]
    # Confidence uses the typed field (keyed by per-OTU feature id).
    assert report.confidence == {"s1:ASV0": 99.8}
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
        confidence={"s1:ASV0": 99.8},
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


def test_flatten_carries_new_optional_fields():
    clinical = dict(
        _CLINICAL,
        patient_name="Somsri Chaiyaphum",
        dob="14 Mar 1978",
        age="48",
        gender="Female",
        received_date="23 Aug 2026 16:40",
        ordering_physician="Dr. Anong Wattana",
        healthcare_provider="Bangkok Central Hospital",
        reason_for_testing="Suspected septic arthritis",
        test_performed="Full-length 16S rRNA detection",
        report_date="24 Aug 2026 09:15",
    )
    report = build_medical_report(_make_experiment(), clinical=clinical)
    data = report_to_medical_data(report)

    assert data["patient_name"] == "Somsri Chaiyaphum"
    assert data["dob"] == "14 Mar 1978"
    assert data["age"] == "48"
    assert data["gender"] == "Female"
    assert data["received_date"] == "23 Aug 2026 16:40"
    assert data["ordering_physician"] == "Dr. Anong Wattana"
    assert data["healthcare_provider"] == "Bangkok Central Hospital"
    assert data["reason_for_testing"] == "Suspected septic arthritis"
    assert data["test_performed"] == "Full-length 16S rRNA detection"
    assert data["report_date"] == "24 Aug 2026 09:15"


def test_flatten_new_optional_fields_default_to_none():
    report = build_medical_report(_make_experiment(), clinical=_CLINICAL)
    data = report_to_medical_data(report)

    for key in (
        "patient_name",
        "dob",
        "age",
        "gender",
        "received_date",
        "ordering_physician",
        "healthcare_provider",
        "reason_for_testing",
        "test_performed",
        "report_date",
    ):
        assert data[key] is None


# ─── end-to-end render ───────────────────────────────────────────────


def test_render_from_report_produces_html_with_report_values():
    report = build_medical_report(
        _make_experiment(),
        clinical=_CLINICAL,
        pathogens={"1496": True},
        confidence={"s1:ASV0": 99.8},
        qc_items=[{"name": "Read quality", "status": "pass"}],
    )
    html = generate_medical_report_from_report(report, report_date="31 Aug 2026")

    assert "16S-2026-000184" in html
    assert "PATHOGEN DETECTED" in html
    assert "Clostridioides difficile" in html
    assert "99.8" in html


# ─── _determine_conclusion (bacterial-DNA result logic) ──────────────


def test_determine_conclusion_no_organisms_is_not_detected():
    assert _determine_conclusion([]) == "not_detected"


def test_determine_conclusion_non_commensal_is_pathogen_detected():
    organisms = [
        {"name": "Escherichia coli", "pathogenic": "commensal (gut, stool)"},
        {"name": "Clostridioides difficile", "pathogenic": "primary"},
    ]
    assert _determine_conclusion(organisms) == "pathogen_detected"


def test_determine_conclusion_unknown_classification_is_pathogen_detected():
    # Unknown (None/empty) classification is treated as non-commensal.
    assert _determine_conclusion([{"name": "X", "pathogenic": None}]) == "pathogen_detected"
    assert _determine_conclusion([{"name": "Y", "pathogenic": ""}]) == "pathogen_detected"


def test_determine_conclusion_all_commensal_is_organism_detected():
    organisms = [
        {"name": "Escherichia coli", "pathogenic": "commensal (gut, stool)"},
        {"name": "Bacteroides fragilis", "pathogenic": "Commensal (gut)"},
    ]
    assert _determine_conclusion(organisms) == "organism_detected"
