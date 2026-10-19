import pandas as pd
import pytest

from hoshi.lib.experiment import META_TAX_IDS
from hoshi.lib.savont_reader import UNMAPPED_TAX_ID, SavontReader

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


# ─── per-feature sequence identity ───────────────────────────────────


def test_feature_sequence_identity_keyed_by_asv_in_percent_range():
    r = SavontReader(SAVONT_SAMPLE_DIR)
    conf = r.feature_sequence_identity()
    assert conf  # non-empty for this sample
    assert all(isinstance(k, str) for k in conf)
    assert all(0.0 <= v <= 100.0 for v in conf.values())


def test_feature_sequence_identity_is_per_asv_identity():
    """Sequence identity must equal each ASV's alignment_identity (first hit per ASV)."""
    r = SavontReader(SAVONT_SAMPLE_DIR)
    conf = r.feature_sequence_identity()

    mapping = pd.read_csv(r.asv_mappings_path, sep="\t")
    first = mapping.drop_duplicates(subset=["asv_header"], keep="first").copy()
    expected = dict(zip(first["asv_header"].astype(str), first["alignment_identity"]))

    for asv, value in conf.items():
        assert value == pytest.approx(float(expected[asv]))



# ─── sub-threshold hits, unaligned ASVs, parity with `savont export` ──

SOIL_DIR = "test_data/savont_output/test_ind/savont-out-soil"
# Produced by `savont export --input-dirs savont-out-soil --output-dir savont-out-soil-export`
# (savont 0.8.0). Export resolves each ASV to the first asv_mappings.tsv row and
# counts reads per lineage string.
SOIL_EXPORT_TAXON_COUNTS = (
    "test_data/savont_output/test_ind/savont-out-soil-export/merged_taxon_counts.tsv"
)
_LINEAGE_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def _first_rows(reader: SavontReader) -> pd.DataFrame:
    mapping = pd.read_csv(reader.asv_mappings_path, sep="\t", dtype=str, keep_default_na=False)
    return mapping.drop_duplicates(subset=["asv_header"], keep="first").set_index("asv_header")


def test_unclassified_species_tax_id_is_tagged_with_unk_suffix():
    """Below-threshold hits keep the nearest tax_id but are tagged so they never
    pass as that species; the suffix counts how many ranks above species are
    also unresolved (000 = only species)."""
    r = SavontReader(SOIL_DIR)
    df = r.per_feature_abundance().set_index("feature_id")
    first = _first_rows(r)

    for asv, hit in first.iterrows():
        unresolved_above = sum(
            hit[rank].startswith("UNCLASSIFIED")
            for rank in ("genus", "family", "order", "class", "phylum", "superkingdom")
        )
        if hit["species"].startswith("UNCLASSIFIED"):
            expected = f"{hit['tax_id']}_UNK{unresolved_above:03d}"
        else:
            expected = hit["tax_id"]
        assert df.loc[asv, "tax_id"] == expected, asv

    tagged = df["tax_id"].str.contains("_UNK", na=False)
    assert (df.loc["final_consensus_0_depth_878", "tax_id"]) == "674703_UNK000"
    assert df.loc[tagged, "tax_id"].str.endswith("_UNK001").any()  # genus also unresolved
    assert not df.loc[~tagged & df["tax_id"].ne("unmapped"), "species"].str.startswith(
        "UNCLASSIFIED"
    ).any()


def test_unaligned_asvs_are_kept_as_unmapped():
    """ASVs with no alignment are absent from asv_mappings.tsv; they must still
    count toward the total (as savont does), flagged with the meta tax_id."""
    r = SavontReader(SOIL_DIR)
    df = r.per_feature_abundance().set_index("feature_id")
    feature_table = pd.read_csv(r.feature_table_path, sep="\t", index_col=0)

    assert set(df.index) == set(feature_table.index)
    assert UNMAPPED_TAX_ID in META_TAX_IDS  # so species_view never lists it
    assert df.loc["final_consensus_242_depth_14", "tax_id"] == "unmapped"
    assert df.loc["final_consensus_242_depth_14", "estimated counts"] == 14
    assert df["estimated counts"].sum() == feature_table.to_numpy().sum()


def test_count_table_matches_savont_export_per_lineage():
    """Read counts and percentages per lineage equal `savont export`."""
    from hoshi.lib.egress import build_count_table  # noqa: PLC0415

    r = SavontReader(SOIL_DIR)
    table = build_count_table(r.to_summarized_experiment())
    lineage = table[_LINEAGE_ORDER].fillna("").astype(str).agg(";".join, axis=1)
    ours = table.groupby(lineage)["estimated_count"].sum().astype(int)
    ours_pct = table.groupby(lineage)["relative_abundance"].sum()

    export = pd.read_csv(SOIL_EXPORT_TAXON_COUNTS, sep="\t", dtype={"taxon": str})
    export_counts = export.set_index("taxon").iloc[:, 0]
    total = export_counts.sum()
    # Unaligned reads are counted in the total but not listed as an organism.
    export_counts = export_counts.drop("Unclassified")

    pd.testing.assert_series_equal(
        ours.sort_index(), export_counts.sort_index(), check_names=False
    )
    pd.testing.assert_series_equal(
        ours_pct.sort_index(), (export_counts / total).sort_index(), check_names=False
    )
