import pandas as pd
import pytest

from hoshi.lib.savont_reader import SavontReader

SAVONT_SAMPLE_DIR = "test_data/savont_output/test_ind/savont-out-sample01"


# ─── construction ────────────────────────────────────────────────────


def test_reader_requires_existing_directory():
    with pytest.raises(ValueError, match="Savont directory not found"):
        SavontReader("test_data/savont_output/does-not-exist")


# ─── per-feature abundance table ─────────────────────────────────────


def test_per_feature_abundance_is_one_row_per_asv():
    """The per-feature table keeps one row per ASV (no species collapsing)."""
    r = SavontReader(SAVONT_SAMPLE_DIR)
    df = r.per_feature_abundance()
    # feature_id (asv_header) is unique per row
    assert df["feature_id"].is_unique
    # abundance is normalised across features
    assert df["abundance"].sum() == pytest.approx(1.0)
    assert df["abundance"].is_monotonic_decreasing


# ─── species-level rollup ────────────────────────────────────────────


def test_species_abundance_does_not_leak_confidence():
    """Confidence (alignment_identity) must NOT leak into the species table."""
    r = SavontReader(SAVONT_SAMPLE_DIR)
    df = r.species_abundance()
    assert "alignment_identity" not in df.columns


def test_species_abundance_sorted_descending_and_sums_to_one():
    r = SavontReader(SAVONT_SAMPLE_DIR)
    df = r.species_abundance()
    assert df["abundance"].is_monotonic_decreasing
    assert df["abundance"].sum() == pytest.approx(1.0)


# ─── per-feature confidence ──────────────────────────────────────────


def test_feature_confidence_keyed_by_asv_in_percent_range():
    r = SavontReader(SAVONT_SAMPLE_DIR)
    conf = r.feature_confidence()
    assert conf  # non-empty for this sample
    assert all(isinstance(k, str) for k in conf)
    assert all(0.0 <= v <= 100.0 for v in conf.values())


def test_feature_confidence_is_per_asv_identity():
    """Confidence must equal each ASV's alignment_identity (first hit per ASV)."""
    r = SavontReader(SAVONT_SAMPLE_DIR)
    conf = r.feature_confidence()

    mapping = pd.read_csv(r.asv_mappings_path, sep="\t")
    first = mapping.drop_duplicates(subset=["asv_header"], keep="first").copy()
    expected = dict(zip(first["asv_header"].astype(str), first["alignment_identity"]))

    for asv, value in conf.items():
        assert value == pytest.approx(float(expected[asv]))
