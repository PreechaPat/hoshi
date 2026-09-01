
import pandas as pd
import pytest

from hoshi.lib.ingress import (
    read_emu_abundance,
    read_savont_abundance,
    read_savont_abundance_into_summarizedexperiment,
)


def test_read_emu_abundance_from_path():
    df = read_emu_abundance("test_data/emu_output/emu-mock01.tsv")

    expected_columns = {
        "tax_id",
        "abundance",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "estimated counts"
    }
    assert expected_columns.issubset(set(df.columns))

    abundance_value = df.loc[df["tax_id"] == "1290", "abundance"].iat[0]
    assert abundance_value == pytest.approx(0.2755103097)
    assert df["estimated counts"].dtype.kind == "f"


def test_read_emu_abundance_dataframe_input_adds_taxonomy_columns():
    raw = pd.DataFrame({"tax_id": ["1", "2"], "abundance": [0.1, 0.9]})
    df = read_emu_abundance(raw)

    for column in ("superkingdom", "phylum", "class", "order", "family", "genus", "species"):
        assert column in df.columns
        assert df[column].isna().all()

    assert df.loc[df["tax_id"] == "1", "abundance"].iat[0] == pytest.approx(0.1)


def test_read_emu_abundance_reorders_and_limits_columns():
    df = read_emu_abundance("test_data/emu_output/emu-mock01.tsv", reorder=True)

    expected_columns = [
        "tax_id",
        "abundance",
        "estimated counts",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]

    assert list(df.columns) == expected_columns
    assert df["estimated counts"].dtype.kind == "i"


# ─── Savont reader tests ─────────────────────────────────────────────────────

SAVONT_SAMPLE_DIR = "test_data/savont_output/test_ind/savont-out-sample01"
SAVONT_FEATURE_TABLE = SAVONT_SAMPLE_DIR + "/feature-table.tsv"
SAVONT_ASV_MAPPING = SAVONT_SAMPLE_DIR + "/asv_mappings.tsv"


def test_read_savont_abundance_returns_expected_columns():
    df = read_savont_abundance(SAVONT_FEATURE_TABLE, SAVONT_ASV_MAPPING)

    expected_columns = [
        "abundance",
        "tax_id",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom",
        "estimated counts",
    ]
    assert list(df.columns) == expected_columns


def test_read_savont_abundance_tax_id_is_regular_column():
    df = read_savont_abundance(SAVONT_FEATURE_TABLE, SAVONT_ASV_MAPPING)

    # tax_id should be a regular column, not the index
    assert "tax_id" in df.columns
    assert df.index.name != "tax_id"


def test_read_savont_abundance_values():
    df = read_savont_abundance(SAVONT_FEATURE_TABLE, SAVONT_ASV_MAPPING)

    # Total reads = 71 + 53 + 47 + 44 + 34 + 33 + 18 = 300
    assert df["estimated counts"].sum() == 300

    # Abundances should sum to 1.0
    assert df["abundance"].sum() == pytest.approx(1.0)

    # Clostridioides difficile (tax_id 1496) has ASVs 0 (71) + 1 (53) = 124 reads
    row = df[df["tax_id"] == "1496"]
    assert len(row) == 1
    assert row["estimated counts"].iat[0] == 124
    assert row["abundance"].iat[0] == pytest.approx(124 / 300)
    assert row["species"].iat[0] == "Clostridioides difficile"


def test_read_savont_abundance_sorted_by_abundance():
    df = read_savont_abundance(SAVONT_FEATURE_TABLE, SAVONT_ASV_MAPPING)

    abundances = df["abundance"].tolist()
    assert abundances == sorted(abundances, reverse=True)


def test_read_savont_abundance_missing_file():
    with pytest.raises(ValueError, match="Feature table not found"):
        read_savont_abundance("/nonexistent/feature-table.tsv", SAVONT_ASV_MAPPING)


def test_savont_into_summarizedexperiment_single_sample():
    se = read_savont_abundance_into_summarizedexperiment(SAVONT_SAMPLE_DIR, sample_names=["sample01"])

    assert se.metadata["source"] == "savont"
    assert se.n_samples == 1
    assert se.sample_ids[0] == "sample01"
    assert se.n_features > 0

    # row_data should have taxonomy columns (not tax_id since it's the index)
    assert "species" in se.row_data.columns
    assert "genus" in se.row_data.columns
    assert "phylum" in se.row_data.columns

    # Abundance matrix should sum to ~1.0 per sample
    assert se.assays["abundance"]["sample01"].sum() == pytest.approx(1.0)


def test_savont_into_summarizedexperiment_custom_sample_name():
    se = read_savont_abundance_into_summarizedexperiment(SAVONT_SAMPLE_DIR, sample_names=["my_sample"])

    assert se.sample_ids[0] == "my_sample"
    assert se.col_data.loc["my_sample", "sample_name"] == "my_sample"


def test_savont_into_summarizedexperiment_name_mismatch_raises():
    with pytest.raises(ValueError, match="Length mismatch"):
        read_savont_abundance_into_summarizedexperiment(SAVONT_SAMPLE_DIR, sample_names=["a", "b"])
