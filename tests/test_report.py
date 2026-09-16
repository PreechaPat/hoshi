import pandas as pd
import pytest

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.report import Report


def _make_experiment() -> SummarizedExperiment:
    """A small 3-OTU, 1-sample per-OTU experiment for report tests."""
    otu_ids = ["s1:ASV0", "s1:ASV1", "s1:ASV2"]
    abundance = pd.DataFrame({"s1": [0.6, 0.3, 0.1]}, index=otu_ids)
    return SummarizedExperiment(
        assays={"abundance": abundance},
        row_data=pd.DataFrame({"tax_id": ["1496", "562", "1351"]}, index=otu_ids),
        col_data=pd.DataFrame(index=["s1"]),
        metadata={"source": "savont"},
    )


# ─── construction ────────────────────────────────────────────────────


def test_report_holds_experiment_confidence_and_metadata():
    se = _make_experiment()
    report = Report(
        experiment=se,
        confidence={"s1:ASV0": 100.0, "s1:ASV1": 98.0, "s1:ASV2": 90.0},
        metadata={"note": "test"},
    )
    assert report.experiment is se
    assert report.confidence == {"s1:ASV0": 100.0, "s1:ASV1": 98.0, "s1:ASV2": 90.0}
    assert report.metadata == {"note": "test"}


# ─── confidence_pct ──────────────────────────────────────────────────


def test_confidence_pct_none_when_no_confidence():
    report = Report(experiment=_make_experiment())
    assert report.confidence_pct is None


def test_confidence_pct_is_max_of_per_otu_confidence():
    # Reports render at species level; confidence collapses per-OTU with max().
    # TODO: revisit — placeholder aggregation (see Report.confidence_pct).
    se = _make_experiment()
    report = Report(
        experiment=se,
        confidence={"s1:ASV0": 100.0, "s1:ASV1": 90.0, "s1:ASV2": 50.0},
    )
    assert report.confidence_pct == pytest.approx(100.0)


def test_confidence_pct_max_without_abundance_assay():
    counts = pd.DataFrame({"s1": [10, 20]}, index=["s1:ASV0", "s1:ASV1"])
    se = SummarizedExperiment(
        assays={"counts": counts},
        row_data=pd.DataFrame({"tax_id": ["1496", "562"]}, index=["s1:ASV0", "s1:ASV1"]),
        col_data=pd.DataFrame(index=["s1"]),
    )
    report = Report(experiment=se, confidence={"s1:ASV0": 100.0, "s1:ASV1": 80.0})
    assert report.confidence_pct == pytest.approx(100.0)


# ─── immutable update helpers ────────────────────────────────────────


def test_with_confidence_returns_new_object_without_mutation():
    report = Report(experiment=_make_experiment())
    updated = report.with_confidence({"s1:ASV0": 99.0})
    assert report.confidence == {}  # original untouched
    assert updated.confidence == {"s1:ASV0": 99.0}
    assert updated is not report


def test_with_metadata_merges_without_mutation():
    report = Report(experiment=_make_experiment(), metadata={"a": 1})
    updated = report.with_metadata(b=2)
    assert report.metadata == {"a": 1}  # original untouched
    assert updated.metadata == {"a": 1, "b": 2}
