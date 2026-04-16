from pathlib import Path

import pandas as pd
import pytest

from hoshi.lib.ingress import read_emu_abundance


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
