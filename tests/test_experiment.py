from pathlib import Path

import pandas as pd
import pytest

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.ingress import emu_to_experiment, _derive_sample_name


# ─── SummarizedExperiment tests ──────────────────────────────────────


def test_empty_experiment():
    se = SummarizedExperiment()
    assert se.n_features == 0
    assert se.n_samples == 0
    assert se.assay_names == []
    assert se.shape == (0, 0)


def test_experiment_construction():
    abundance = pd.DataFrame(
        {"s1": [0.5, 0.3, 0.2], "s2": [0.1, 0.6, 0.3]},
        index=["tax_a", "tax_b", "tax_c"],
    )
    row_data = pd.DataFrame(
        {"genus": ["Genus_A", "Genus_B", "Genus_C"]},
        index=["tax_a", "tax_b", "tax_c"],
    )
    col_data = pd.DataFrame(
        {"condition": ["control", "treatment"]},
        index=["s1", "s2"],
    )

    se = SummarizedExperiment(
        assays={"abundance": abundance},
        row_data=row_data,
        col_data=col_data,
        metadata={"source": "test"},
    )

    assert se.n_features == 3
    assert se.n_samples == 2
    assert se.assay_names == ["abundance"]
    assert list(se.feature_ids) == ["tax_a", "tax_b", "tax_c"]
    assert list(se.sample_ids) == ["s1", "s2"]
    assert se.metadata == {"source": "test"}


def test_experiment_validation_mismatched_shapes():
    a1 = pd.DataFrame({"s1": [1, 2]}, index=["f1", "f2"])
    a2 = pd.DataFrame({"s1": [1, 2, 3]}, index=["f1", "f2", "f3"])

    with pytest.raises(ValueError, match="shape"):
        SummarizedExperiment(assays={"a1": a1, "a2": a2})


def test_experiment_validation_mismatched_row_data():
    assay = pd.DataFrame({"s1": [1, 2]}, index=["f1", "f2"])
    row_data = pd.DataFrame({"col": ["a", "b"]}, index=["f1", "f3"])

    with pytest.raises(ValueError, match="row_data index"):
        SummarizedExperiment(assays={"a": assay}, row_data=row_data)


def test_experiment_validation_mismatched_col_data():
    assay = pd.DataFrame({"s1": [1, 2]}, index=["f1", "f2"])
    col_data = pd.DataFrame({"col": ["x"]}, index=["s_wrong"])

    with pytest.raises(ValueError, match="col_data index"):
        SummarizedExperiment(assays={"a": assay}, col_data=col_data)


def test_experiment_subset_features():
    assay = pd.DataFrame(
        {"s1": [1, 2, 3], "s2": [4, 5, 6]},
        index=["f1", "f2", "f3"],
    )
    row_data = pd.DataFrame({"name": ["A", "B", "C"]}, index=["f1", "f2", "f3"])

    se = SummarizedExperiment(assays={"counts": assay}, row_data=row_data)
    sub = se.subset_features(["f1", "f3"])

    assert sub.n_features == 2
    assert sub.n_samples == 2
    assert list(sub.feature_ids) == ["f1", "f3"]
    assert list(sub.row_data["name"]) == ["A", "C"]


def test_experiment_subset_samples():
    assay = pd.DataFrame(
        {"s1": [1, 2], "s2": [3, 4], "s3": [5, 6]},
        index=["f1", "f2"],
    )
    col_data = pd.DataFrame(
        {"batch": ["A", "A", "B"]}, index=["s1", "s2", "s3"]
    )

    se = SummarizedExperiment(assays={"counts": assay}, col_data=col_data)
    sub = se.subset_samples(["s1", "s3"])

    assert sub.n_samples == 2
    assert list(sub.sample_ids) == ["s1", "s3"]
    assert list(sub.col_data["batch"]) == ["A", "B"]


# ─── _derive_sample_name tests ───────────────────────────────────────


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("sample01_rel-abundance.tsv", "sample01"),
        ("barcode11.fastq_rel-abundance.tsv", "barcode11"),
        ("emu-mock01.tsv", "emu-mock01"),
    ],
)
def test_derive_sample_name(filename, expected):
    assert _derive_sample_name(Path(filename)) == expected


# ─── emu_to_experiment tests ─────────────────────────────────────────


def test_emu_to_experiment_single_sample():
    se = emu_to_experiment(
        "test_data/emu_output/emu-mock01.tsv",
        sample_names=["mock01"],
    )

    assert se.n_samples == 1
    assert se.n_features > 0
    assert "abundance" in se.assay_names
    assert "counts" in se.assay_names
    assert list(se.sample_ids) == ["mock01"]
    assert se.metadata == {"source": "emu"}

    # Check that row_data has taxonomy
    assert "genus" in se.row_data.columns
    assert "species" in se.row_data.columns

    # Check col_data
    assert "source_file" in se.col_data.columns

    # Check abundance sums roughly to 1
    total = se.assays["abundance"]["mock01"].sum()
    assert total == pytest.approx(1.0, abs=0.01)


def test_emu_to_experiment_multi_sample():
    se = emu_to_experiment([
        "test_data/emu_output/test_ind/sample01/sample01_rel-abundance.tsv",
        "test_data/emu_output/test_ind/sample02/sample02_rel-abundance.tsv",
        "test_data/emu_output/test_ind/sample03/sample03_rel-abundance.tsv",
    ])

    assert se.n_samples == 3
    assert list(se.sample_ids) == ["sample01", "sample02", "sample03"]
    assert se.n_features >= 11  # union of all features across samples

    # Features missing from a sample should be filled with 0
    abundance = se.assays["abundance"]
    assert (abundance >= 0).all().all()


def test_emu_to_experiment_auto_derives_sample_names():
    se = emu_to_experiment(
        "test_data/emu_output/silva_count/barcode11.fastq_rel-abundance.tsv"
    )
    assert list(se.sample_ids) == ["barcode11"]


def test_emu_to_experiment_mismatched_names_raises():
    with pytest.raises(ValueError, match="Length mismatch"):
        emu_to_experiment(
            ["test_data/emu_output/emu-mock01.tsv"],
            sample_names=["a", "b"],
        )
