"""Tests for the ``hoshi convert`` command.

Focuses on the CLI glue: format selection, file vs. stdout output, and the
Savont → species-count-table path added for issue #1.
"""

import pandas as pd

from hoshi import cli

SAVONT_SAMPLE_DIR = "test_data/savont_output/test_ind/savont-out-sample01"
EMU_SAMPLE = "test_data/emu_output/test_ind/sample02/sample02_rel-abundance.tsv"

_COUNT_TABLE_COLUMNS = [
    "relative_abundance",
    "estimated_count",
    "tax_id",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
]


def test_convert_savont_to_table_writes_file(tmp_path):
    out = tmp_path / "species_counts.tsv"
    code = cli.main(
        [
            "convert",
            "--input-format", "savont",
            "--output-format", "table",
            SAVONT_SAMPLE_DIR,
            "-o", str(out),
        ]
    )
    assert code == 0
    assert out.exists()

    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _COUNT_TABLE_COLUMNS
    # Savont sample has 300 reads total; tax_id + counts are restored.
    assert df["estimated_count"].sum() == 300
    assert all(str(int(t)).isdigit() for t in df["tax_id"])


def test_convert_savont_to_table_stdout(capsys):
    code = cli.main(
        [
            "convert",
            "--input-format", "savont",
            "--output-format", "table",
            SAVONT_SAMPLE_DIR,
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    header = out.splitlines()[0]
    assert header.split("\t") == _COUNT_TABLE_COLUMNS


def test_convert_savont_custom_name(tmp_path):
    """The -n/--name flag is accepted for savont input."""
    out = tmp_path / "named.tsv"
    code = cli.main(
        [
            "convert",
            "--input-format", "savont",
            "--output-format", "table",
            SAVONT_SAMPLE_DIR,
            "-n", "my_sample",
            "-o", str(out),
        ]
    )
    assert code == 0
    assert out.exists()


def test_convert_missing_input_returns_error(capsys):
    code = cli.main(
        [
            "convert",
            "--input-format", "savont",
            "--output-format", "table",
            "/nonexistent/savont-dir",
        ]
    )
    assert code == 1
    assert "not found" in capsys.readouterr().err


def test_convert_emu_to_table_still_works(tmp_path):
    """The new table output also works from the existing emu input path."""
    out = tmp_path / "emu_table.tsv"
    code = cli.main(
        [
            "convert",
            "--output-format", "table",
            EMU_SAMPLE,
            "-o", str(out),
        ]
    )
    assert code == 0
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _COUNT_TABLE_COLUMNS
