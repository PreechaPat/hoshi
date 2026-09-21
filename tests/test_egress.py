"""Tests for hoshi.lib.egress — Kraken2 report conversion."""

import pandas as pd
import pytest

from hoshi.lib.egress import experiment_to_kraken2, write_kraken2_report
from hoshi.lib.egress import (
    experiment_to_count_table,
    build_count_table,
    write_count_table,
)
from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.ingress import read_emu_abundance_into_summarizedexperiment


# ─── Fixtures ────────────────────────────────────────────────────────


@pytest.fixture
def simple_experiment():
    """A minimal single-sample per-OTU experiment for testing.

    Feature index is per-OTU (one OTU per species here); tax_id lives in
    row_data as a nullable annotation.
    """
    otu_ids = ["sample1:ASV0", "sample1:ASV1", "sample1:ASV2"]
    counts = pd.DataFrame({"sample1": [100.0, 200.0, 50.0]}, index=otu_ids)
    abundance = pd.DataFrame({"sample1": [0.286, 0.571, 0.143]}, index=otu_ids)
    row_data = pd.DataFrame(
        {
            "tax_id": ["1351", "1496", "817"],
            "superkingdom": ["Bacteria", "Bacteria", "Bacteria"],
            "phylum": ["Bacillota", "Bacillota", "Bacteroidota"],
            "class": ["Bacilli", "Clostridia", "Bacteroidia"],
            "order": ["Lactobacillales", "Peptostreptococcales", "Bacteroidales"],
            "family": ["Enterococcaceae", "Peptostreptococcaceae", "Bacteroidaceae"],
            "genus": ["Enterococcus", "Clostridioides", "Bacteroides"],
            "species": [
                "Enterococcus faecalis",
                "Clostridioides difficile",
                "Bacteroides fragilis",
            ],
        },
        index=otu_ids,
    )
    col_data = pd.DataFrame(
        {"sample_name": ["sample1"]},
        index=["sample1"],
    )
    return SummarizedExperiment(
        assays={"abundance": abundance, "counts": counts},
        row_data=row_data,
        col_data=col_data,
        metadata={"source": "test"},
    )


@pytest.fixture
def real_experiment():
    """Load a real Emu sample from test data."""
    return read_emu_abundance_into_summarizedexperiment(
        "test_data/emu_output/test_ind/sample02/sample02_rel-abundance.tsv"
    )


# ─── Basic conversion tests ─────────────────────────────────────────


def test_experiment_to_kraken2_tab_delimited(simple_experiment):
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")
    for line in lines:
        fields = line.split("\t")
        assert len(fields) == 6, f"Expected 6 fields, got {len(fields)}: {line}"


def test_experiment_to_kraken2_has_unclassified_line(simple_experiment):
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")
    # First line should be unclassified
    fields = lines[0].split("\t")
    assert fields[3] == "U"
    assert "unclassified" in fields[5]


def test_experiment_to_kraken2_rank_codes(simple_experiment):
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")

    rank_codes_found = set()
    for line in lines:
        fields = line.split("\t")
        rank_codes_found.add(fields[3])

    # Should have at least: U (unclassified), D (domain/superkingdom),
    # P (phylum), and S (species)
    assert "U" in rank_codes_found
    assert "D" in rank_codes_found
    assert "P" in rank_codes_found
    assert "S" in rank_codes_found


def test_experiment_to_kraken2_species_direct_counts(simple_experiment):
    """Species-level entries should have non-zero direct counts."""
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")

    species_lines = [line for line in lines if line.split("\t")[3] == "S"]
    assert len(species_lines) == 3  # 3 species in our fixture

    for line in species_lines:
        fields = line.split("\t")
        direct_count = int(fields[2])
        assert direct_count > 0


def test_experiment_to_kraken2_higher_ranks_zero_direct(simple_experiment):
    """Non-species ranks should have 0 direct counts."""
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")

    for line in lines:
        fields = line.split("\t")
        rank_code = fields[3]
        direct_count = int(fields[2])
        if rank_code not in ("S", "U"):
            assert direct_count == 0, (
                f"Rank {rank_code} should have 0 direct counts, got {direct_count}"
            )


def test_experiment_to_kraken2_clade_aggregation(simple_experiment):
    """Clade counts should aggregate child counts properly."""
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")

    # Find the Bacteria (superkingdom/domain) line — should have all 350 counts
    domain_lines = [line for line in lines if line.split("\t")[3] == "D"]
    assert len(domain_lines) == 1  # Only one domain: Bacteria

    fields = domain_lines[0].split("\t")
    clade_count = int(fields[1])
    assert clade_count == 350  # 100 + 200 + 50

    # Bacillota phylum should have 300 (100 + 200)
    phylum_lines = [line for line in lines if line.split("\t")[3] == "P"]
    bacillota_line = [line for line in phylum_lines if "Bacillota" in line]
    assert len(bacillota_line) == 1
    fields = bacillota_line[0].split("\t")
    assert int(fields[1]) == 300


def test_experiment_to_kraken2_percentage_sums(simple_experiment):
    """Percentages at the top level should sum to ~100 (excluding unclassified)."""
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")

    # Domain-level lines should sum to 100%
    domain_lines = [line for line in lines if line.split("\t")[3] == "D"]
    total_pct = sum(float(line.split("\t")[0]) for line in domain_lines)
    assert total_pct == pytest.approx(100.0, abs=0.1)


# ─── Error handling tests ────────────────────────────────────────────


def test_experiment_to_kraken2_no_counts_assay_raises():
    """Should raise if no 'counts' assay exists."""
    abundance = pd.DataFrame({"s1": [0.5, 0.5]}, index=["f1", "f2"])
    se = SummarizedExperiment(assays={"abundance": abundance})

    with pytest.raises(ValueError, match="counts"):
        experiment_to_kraken2(se)


def test_experiment_to_kraken2_multi_sample_no_selection_raises():
    """Should raise if multiple samples and no sample specified."""
    counts = pd.DataFrame(
        {"s1": [10, 20], "s2": [30, 40]}, index=["f1", "f2"]
    )
    row_data = pd.DataFrame(
        {"species": ["Sp A", "Sp B"]}, index=["f1", "f2"]
    )
    se = SummarizedExperiment(
        assays={"counts": counts}, row_data=row_data
    )

    with pytest.raises(ValueError, match="samples"):
        experiment_to_kraken2(se)


def test_experiment_to_kraken2_invalid_sample_raises(simple_experiment):
    """Should raise if the specified sample doesn't exist."""
    with pytest.raises(ValueError, match="not found"):
        experiment_to_kraken2(simple_experiment, sample="nonexistent")


def test_experiment_to_kraken2_no_taxonomy_raises():
    """Should raise if row_data has no taxonomy columns."""
    counts = pd.DataFrame({"s1": [10, 20]}, index=["f1", "f2"])
    row_data = pd.DataFrame({"other_col": ["x", "y"]}, index=["f1", "f2"])
    se = SummarizedExperiment(assays={"counts": counts}, row_data=row_data)

    with pytest.raises(ValueError, match="taxonomy column"):
        experiment_to_kraken2(se)


def test_experiment_to_kraken2_empty_counts_returns_empty():
    """All-zero counts should return empty string."""
    counts = pd.DataFrame({"s1": [0.0, 0.0]}, index=["f1", "f2"])
    row_data = pd.DataFrame(
        {"species": ["Sp A", "Sp B"]}, index=["f1", "f2"]
    )
    se = SummarizedExperiment(assays={"counts": counts}, row_data=row_data)

    result = experiment_to_kraken2(se)
    assert result == ""


# ─── Real data integration tests ────────────────────────────────────


def test_real_emu_to_kraken2(real_experiment):
    """Convert a real Emu output file and verify structure."""
    result = experiment_to_kraken2(real_experiment)

    lines = result.strip().split("\n")
    assert len(lines) > 5  # Should have multiple taxonomy lines

    # All lines should be tab-delimited with 6 fields
    for line in lines:
        fields = line.split("\t")
        assert len(fields) == 6

    # First line is unclassified
    assert lines[0].split("\t")[3] == "U"

    # Should contain known species from sample02
    assert any("Streptococcus agalactiae" in line for line in lines)
    assert any("Lactobacillus gasseri" in line for line in lines)


def test_real_emu_to_kraken2_tax_ids(real_experiment):
    """Tax IDs at species level should be numeric NCBI IDs."""
    result = experiment_to_kraken2(real_experiment)
    lines = result.strip().split("\n")

    species_lines = [line for line in lines if line.split("\t")[3] == "S"]
    for line in species_lines:
        fields = line.split("\t")
        tax_id = fields[4]
        assert tax_id.isdigit(), f"Species tax_id should be numeric, got: {tax_id}"


def test_real_emu_to_kraken2_counts_consistent(real_experiment):
    """Clade counts at top-level rank should equal total estimated counts."""
    result = experiment_to_kraken2(real_experiment)
    lines = result.strip().split("\n")

    # Find the top-level rank (first non-U rank that appears)
    # In sample02, superkingdom is blank so root is phylum
    non_u_lines = [line for line in lines if line.split("\t")[3] != "U"]
    assert len(non_u_lines) > 0, "Should have taxonomy output"

    # The top-level rank is the first rank code that appears
    top_rank = non_u_lines[0].split("\t")[3]
    top_lines = [line for line in lines if line.split("\t")[3] == top_rank]
    total_from_report = sum(int(line.split("\t")[1]) for line in top_lines)

    # Compare with actual total counts from the experiment
    total_counts = int(round(
        real_experiment.assays["counts"][real_experiment.sample_ids[0]].sum()
    ))
    assert total_from_report == total_counts


# ─── write_kraken2_report tests ──────────────────────────────────────


def test_write_kraken2_report(simple_experiment, tmp_path):
    """Should write report to file."""
    output_file = tmp_path / "report.txt"
    write_kraken2_report(simple_experiment, output_file)

    assert output_file.exists()
    content = output_file.read_text()
    assert len(content) > 0
    assert "\t" in content

    # Should match what experiment_to_kraken2 returns
    expected = experiment_to_kraken2(simple_experiment)
    assert content == expected


# ─── Species count table (per-OTU egress unit tests) ─────────────────
#
# The experiment feature axis is per-OTU; the count table aggregates OTUs up to
# species by tax_id. These build hand-made per-OTU SEs so they stay
# classifier-agnostic. The Savont end-to-end path is covered in test_ingress.py.

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


def _per_otu_se(counts_map, tax_ids, taxonomy, *, samples=("s1",)):
    """Build a per-OTU SE: index = otu ids, tax_id + taxonomy in row_data."""
    otu_ids = list(counts_map.keys())
    counts = pd.DataFrame({samples[0]: list(counts_map.values())}, index=otu_ids)
    row = {"tax_id": [tax_ids[o] for o in otu_ids]}
    for rank, values in taxonomy.items():
        row[rank] = [values[o] for o in otu_ids]
    row_data = pd.DataFrame(row, index=otu_ids)
    total = counts[samples[0]].sum()
    abundance = pd.DataFrame(
        {samples[0]: counts[samples[0]] / total if total else 0.0}, index=otu_ids
    )
    return SummarizedExperiment(
        assays={"abundance": abundance, "counts": counts}, row_data=row_data
    )


@pytest.fixture
def per_otu_experiment():
    """Three OTUs; two share tax_id 1496 (should collapse to one species)."""
    return _per_otu_se(
        counts_map={"s1:ASV0": 100.0, "s1:ASV1": 24.0, "s1:ASV2": 50.0},
        tax_ids={"s1:ASV0": "1496", "s1:ASV1": "1496", "s1:ASV2": "817"},
        taxonomy={
            "species": {
                "s1:ASV0": "Clostridioides difficile",
                "s1:ASV1": "Clostridioides difficile",
                "s1:ASV2": "Bacteroides fragilis",
            },
            "genus": {
                "s1:ASV0": "Clostridioides",
                "s1:ASV1": "Clostridioides",
                "s1:ASV2": "Bacteroides",
            },
            "superkingdom": {
                "s1:ASV0": "Bacteria",
                "s1:ASV1": "Bacteria",
                "s1:ASV2": "Bacteria",
            },
        },
    )


def test_count_table_aggregates_otus_to_species_by_tax_id(per_otu_experiment):
    """OTUs sharing a tax_id collapse; counts sum; no OTU column exposed."""
    df = build_count_table(per_otu_experiment).set_index("tax_id")

    assert "otu_id" not in df.columns
    assert set(df.index) == {"1496", "817"}  # two ASVs of 1496 merged
    assert df.loc["1496", "estimated_count"] == pytest.approx(124.0)  # 100 + 24
    assert df.loc["817", "estimated_count"] == pytest.approx(50.0)


def test_count_table_keeps_taxonomy_uncollapsed(per_otu_experiment):
    df = build_count_table(per_otu_experiment).set_index("tax_id")
    assert df.loc["1496", "species"] == "Clostridioides difficile"
    assert df.loc["1496", "genus"] == "Clostridioides"
    assert df.loc["817", "superkingdom"] == "Bacteria"


def test_count_table_sorted_by_abundance_descending(per_otu_experiment):
    df = build_count_table(per_otu_experiment)
    assert df["relative_abundance"].is_monotonic_decreasing
    assert df.iloc[0]["tax_id"] == "1496"  # 124 reads > 50


def test_count_table_blank_tax_id_rows_stay_separate():
    """OTUs with no tax_id are not collapsed together (option iii)."""
    se = _per_otu_se(
        counts_map={"s1:ASV0": 10.0, "s1:ASV1": 20.0, "s1:ASV2": 30.0},
        tax_ids={"s1:ASV0": "1496", "s1:ASV1": pd.NA, "s1:ASV2": pd.NA},
        taxonomy={
            "species": {
                "s1:ASV0": "Clostridioides difficile",
                "s1:ASV1": "",
                "s1:ASV2": "",
            },
        },
    )
    df = build_count_table(se)
    # one collapsed row for 1496 + two separate blank-tax_id rows
    assert len(df) == 3
    blanks = df[df["tax_id"].astype(str).isin(["", "<NA>"])]
    assert len(blanks) == 2
    # counts still total correctly
    assert df["estimated_count"].astype(float).sum() == pytest.approx(60.0)


def test_count_table_derives_abundance_when_assay_absent():
    counts = pd.DataFrame({"s1": [30.0, 10.0]}, index=["s1:ASV0", "s1:ASV1"])
    row_data = pd.DataFrame(
        {"tax_id": ["10", "20"], "species": ["Sp A", "Sp B"], "genus": ["A", "B"]},
        index=["s1:ASV0", "s1:ASV1"],
    )
    se = SummarizedExperiment(assays={"counts": counts}, row_data=row_data)

    df = build_count_table(se).set_index("tax_id")
    assert df.loc["10", "relative_abundance"] == pytest.approx(0.75)
    assert df.loc["20", "relative_abundance"] == pytest.approx(0.25)


def test_count_table_no_counts_assay_raises():
    abundance = pd.DataFrame({"s1": [0.5, 0.5]}, index=["f1", "f2"])
    se = SummarizedExperiment(assays={"abundance": abundance})
    with pytest.raises(ValueError, match="counts"):
        build_count_table(se)


def test_count_table_no_taxonomy_raises():
    counts = pd.DataFrame({"s1": [10, 20]}, index=["f1", "f2"])
    row_data = pd.DataFrame({"other_col": ["x", "y"]}, index=["f1", "f2"])
    se = SummarizedExperiment(assays={"counts": counts}, row_data=row_data)
    with pytest.raises(ValueError, match="taxonomy column"):
        build_count_table(se)


def test_count_table_multi_sample_no_selection_raises():
    counts = pd.DataFrame({"s1": [10, 20], "s2": [30, 40]}, index=["f1", "f2"])
    row_data = pd.DataFrame(
        {"tax_id": ["1", "2"], "species": ["Sp A", "Sp B"]}, index=["f1", "f2"]
    )
    se = SummarizedExperiment(assays={"counts": counts}, row_data=row_data)
    with pytest.raises(ValueError, match="samples"):
        build_count_table(se)


def test_count_table_sequence_identity_na_when_absent(per_otu_experiment):
    """No sequence-identity metadata (e.g. EMU) → column present but N/A."""
    df = build_count_table(per_otu_experiment)
    assert list(df.columns) == _COUNT_TABLE_COLUMNS
    assert df["sequence_identity"].isna().all()


def test_count_table_emits_sequence_identity_aggregated_by_max():
    """When the experiment carries per-OTU sequence identity, emit it per species,
    aggregated with max() over the OTUs that collapse into that species."""
    se = _per_otu_se(
        counts_map={"s1:ASV0": 100.0, "s1:ASV1": 24.0, "s1:ASV2": 50.0},
        tax_ids={"s1:ASV0": "1496", "s1:ASV1": "1496", "s1:ASV2": "817"},
        taxonomy={
            "species": {
                "s1:ASV0": "Clostridioides difficile",
                "s1:ASV1": "Clostridioides difficile",
                "s1:ASV2": "Bacteroides fragilis",
            },
        },
    )
    se.metadata["sequence_identity"] = {
        "s1": {"s1:ASV0": 88.0, "s1:ASV1": 97.0, "s1:ASV2": 99.5}
    }
    df = build_count_table(se).set_index("tax_id")

    assert "sequence_identity" in df.columns
    # 1496 collapses ASV0 (88) + ASV1 (97) → max = 97.0
    assert df.loc["1496", "sequence_identity"] == pytest.approx(97.0)
    assert df.loc["817", "sequence_identity"] == pytest.approx(99.5)


def test_experiment_to_count_table_is_tab_delimited(per_otu_experiment):
    tsv = experiment_to_count_table(per_otu_experiment)
    header = tsv.splitlines()[0]
    assert header.split("\t") == _COUNT_TABLE_COLUMNS


def test_write_count_table(per_otu_experiment, tmp_path):
    output_file = tmp_path / "species_counts.tsv"
    write_count_table(per_otu_experiment, output_file)

    assert output_file.exists()
    assert output_file.read_text() == experiment_to_count_table(per_otu_experiment)
