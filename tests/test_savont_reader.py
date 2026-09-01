import pandas as pd
import pytest

from hoshi.lib.savont_reader import SavontReader

SAVONT_SAMPLE_DIR = "test_data/savont_output/test_ind/savont-out-sample01"


# ─── construction ────────────────────────────────────────────────────


def test_reader_requires_existing_directory():
    with pytest.raises(ValueError, match="Savont directory not found"):
        SavontReader("test_data/savont_output/does-not-exist")


def test_reader_exposes_fixed_layout_paths():
    r = SavontReader(SAVONT_SAMPLE_DIR)
    assert r.feature_table_path.name == "feature-table.tsv"
    assert r.asv_mappings_path.name == "asv_mappings.tsv"
    assert r.final_asvs_path.name == "final_asvs.fasta"
    assert r.final_clusters_path.name == "final_clusters.tsv"
    # all under the given directory
    assert r.feature_table_path.parent == r.directory


# ─── abundance table ─────────────────────────────────────────────────


def test_abundance_has_expected_columns_and_no_identity_column():
    r = SavontReader(SAVONT_SAMPLE_DIR)
    df = r.abundance()
    expected = [
        "abundance", "tax_id", "species", "genus", "family",
        "order", "class", "phylum", "superkingdom", "estimated counts",
    ]
    assert list(df.columns) == expected
    # confidence must NOT leak into the core abundance table
    assert "alignment_identity" not in df.columns


def test_abundance_sorted_descending_and_sums_to_one():
    r = SavontReader(SAVONT_SAMPLE_DIR)
    df = r.abundance()
    assert df["abundance"].is_monotonic_decreasing
    assert df["abundance"].sum() == pytest.approx(1.0)


# ─── species confidence (max per species) ────────────────────────────


def test_species_confidence_keyed_by_tax_id_in_percent_range():
    r = SavontReader(SAVONT_SAMPLE_DIR)
    conf = r.species_confidence()
    assert conf  # non-empty for this sample
    assert all(isinstance(k, str) for k in conf)
    assert all(0.0 <= v <= 100.0 for v in conf.values())


def test_species_confidence_is_max_per_species():
    """Confidence must equal the MAX alignment_identity per tax_id."""
    r = SavontReader(SAVONT_SAMPLE_DIR)
    conf = r.species_confidence()

    # Recompute expected max directly from the mapping file (first hit per ASV).
    mapping = pd.read_csv(r.asv_mappings_path, sep="\t")
    first = mapping.drop_duplicates(subset=["asv_header"], keep="first").copy()
    first["tax_id"] = first["tax_id"].astype(int).astype(str)
    expected = first.groupby("tax_id")["alignment_identity"].max()

    for tax_id, value in conf.items():
        assert value == pytest.approx(float(expected[tax_id]))
