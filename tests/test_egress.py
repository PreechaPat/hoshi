"""Tests for hoshi.lib.egress — Kraken2 report conversion."""

import pandas as pd
import pytest

from hoshi.lib.egress import experiment_to_kraken2, write_kraken2_report
from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.ingress import read_emu_abundance_into_summarizedexperiment


# ─── Fixtures ────────────────────────────────────────────────────────


@pytest.fixture
def simple_experiment():
    """A minimal single-sample experiment for testing."""
    counts = pd.DataFrame(
        {"sample1": [100.0, 200.0, 50.0]},
        index=["1351", "1496", "817"],
    )
    abundance = pd.DataFrame(
        {"sample1": [0.286, 0.571, 0.143]},
        index=["1351", "1496", "817"],
    )
    row_data = pd.DataFrame(
        {
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
        index=["1351", "1496", "817"],
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


def test_experiment_to_kraken2_returns_string(simple_experiment):
    result = experiment_to_kraken2(simple_experiment)
    assert isinstance(result, str)
    assert len(result) > 0


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


def test_experiment_to_kraken2_indentation(simple_experiment):
    """Deeper ranks should be more indented than shallower ones."""
    result = experiment_to_kraken2(simple_experiment)
    lines = result.strip().split("\n")

    # Domain lines have 2 leading spaces (depth=1)
    # Phylum lines have 4 leading spaces (depth=2)
    # Species lines have 14 leading spaces (depth=7)
    for line in lines:
        fields = line.split("\t")
        rank_code = fields[3]
        name = fields[5]
        leading_spaces = len(name) - len(name.lstrip())

        if rank_code == "U":
            assert leading_spaces == 0
        elif rank_code == "D":
            assert leading_spaces == 2
        elif rank_code == "P":
            assert leading_spaces == 4


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
