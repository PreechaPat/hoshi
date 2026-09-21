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
    "sequence_identity",
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
    # Savont supplies per-OTU alignment identity, so sequence_identity is emitted.
    assert list(df.columns) == _COUNT_TABLE_COLUMNS
    # Savont sample has 300 reads total; tax_id + counts are restored.
    assert df["estimated_count"].sum() == 300
    assert all(str(int(t)).isdigit() for t in df["tax_id"])
    # Sequence identity is a percentage in 0–100 for every species row.
    assert df["sequence_identity"].notna().all()
    assert ((df["sequence_identity"] >= 0) & (df["sequence_identity"] <= 100)).all()
