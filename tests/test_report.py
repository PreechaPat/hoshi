import pandas as pd
import pytest

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.report import Report


def _make_experiment() -> SummarizedExperiment:
    """A small 3-feature, 1-sample experiment for report tests."""
    abundance = pd.DataFrame(
        {"s1": [0.6, 0.3, 0.1]},
        index=["1496", "562", "1351"],
    )
    return SummarizedExperiment(
        assays={"abundance": abundance},
        row_data=pd.DataFrame(index=["1496", "562", "1351"]),
        col_data=pd.DataFrame(index=["s1"]),
        metadata={"source": "savont"},
    )


# ─── construction ────────────────────────────────────────────────────


def test_report_holds_experiment_confidence_and_metadata():
    se = _make_experiment()
    report = Report(
        experiment=se,
        confidence={"1496": 100.0, "562": 98.0, "1351": 90.0},
        metadata={"note": "test"},
    )
    assert report.experiment is se
    assert report.confidence == {"1496": 100.0, "562": 98.0, "1351": 90.0}
    assert report.metadata == {"note": "test"}


# ─── confidence_pct ──────────────────────────────────────────────────


def test_confidence_pct_none_when_no_confidence():
    report = Report(experiment=_make_experiment())
    assert report.confidence_pct is None


def test_confidence_pct_is_abundance_weighted():
    se = _make_experiment()  # abundances: 1496=0.6, 562=0.3, 1351=0.1
    report = Report(
        experiment=se,
        confidence={"1496": 100.0, "562": 90.0, "1351": 50.0},
    )
    # weighted mean = 0.6*100 + 0.3*90 + 0.1*50 = 60 + 27 + 5 = 92.0
    assert report.confidence_pct == pytest.approx(92.0)


def test_confidence_pct_falls_back_to_plain_mean_without_abundance():
    # Experiment with no "abundance" assay -> weights unavailable.
    counts = pd.DataFrame({"s1": [10, 20]}, index=["1496", "562"])
    se = SummarizedExperiment(
        assays={"counts": counts},
        row_data=pd.DataFrame(index=["1496", "562"]),
        col_data=pd.DataFrame(index=["s1"]),
    )
    report = Report(experiment=se, confidence={"1496": 100.0, "562": 80.0})
    # plain mean = (100 + 80) / 2 = 90.0
    assert report.confidence_pct == pytest.approx(90.0)


# ─── immutable update helpers ────────────────────────────────────────


def test_with_confidence_returns_new_object_without_mutation():
    report = Report(experiment=_make_experiment())
    updated = report.with_confidence({"1496": 99.0})
    assert report.confidence == {}  # original untouched
    assert updated.confidence == {"1496": 99.0}
    assert updated is not report


def test_with_metadata_merges_without_mutation():
    report = Report(experiment=_make_experiment(), metadata={"a": 1})
    updated = report.with_metadata(b=2)
    assert report.metadata == {"a": 1}  # original untouched
    assert updated.metadata == {"a": 1, "b": 2}
