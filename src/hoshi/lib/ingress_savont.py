"""Savont ingress — read Savont output directories into a ``SummarizedExperiment``.

Savont folder / species quirk
-----------------------------
Savont output has the same fixed directory structure (no sample-specific information),
so the sample name must be supplied explicitly via ``sample_names``.

Each Savont sample directory contains, after the initial (denoise/cluster)
run::

    feature-table.tsv    # ASV read counts            (used here)
    final_asvs.fasta     # representative ASV sequences
    final_clusters.tsv   # cluster membership

and, after the classify run, additional taxonomy outputs including::

    asv_mappings.tsv        # per-ASV taxonomy assignments  (used here)
    genus_abundance.tsv     # relative abundance at genus level
    species_abundance.tsv   # relative abundance at species level

Species resolution quirk
------------------------
1. The ASV mapping file can list multiple taxonomy hits per ASV. We resolve each
ASV to a **single species by taking the first hit** (``keep="first"``), then
aggregate read counts by ``tax_id`` to compute relative abundance.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.savont_reader import SavontReader


def read_savont_abundance(
    feature_table: str | Path,
    asv_mapping: str | Path,
) -> pd.DataFrame:
    """
    Load Savont feature table and ASV mapping to build a species-level abundance table.

    Reads the feature table (ASV read counts) and ASV mapping (taxonomy
    assignments per ASV), resolves each ASV to a single species (first hit),
    and aggregates read counts by species to compute relative abundance.

    This is a thin wrapper over :class:`SavontReader`. It accepts explicit file
    paths (rather than a directory) for backward compatibility; the two files
    must live in the same directory (the standard Savont layout).

    Parameters
    ----------
    feature_table : str or Path
        Path to ``feature-table.tsv`` containing ASV read counts.

    asv_mapping : str or Path
        Path to ``asv_mappings.tsv`` containing taxonomy assignments per ASV.

    Returns
    -------
    pd.DataFrame
        Species-level abundance table with columns:
        abundance, tax_id, species, genus, family, order, class, phylum,
        superkingdom, estimated counts.
        Indexed by a positional integer index (tax_id is a regular column).

    Raises
    ------
    ValueError
        If required files are missing or cannot be parsed.
    """
    feature_path = Path(feature_table)
    mapping_path = Path(asv_mapping)

    if not feature_path.is_file():
        raise ValueError(f"Feature table not found: {feature_path}")
    if not mapping_path.is_file():
        raise ValueError(f"ASV mapping not found: {mapping_path}")

    return SavontReader(feature_path.parent).abundance()


def read_savont_abundance_into_summarizedexperiment(
    input_dirs: str | Path | list[str | Path],
    *,
    sample_names: list[str],
) -> SummarizedExperiment:
    """
    Load one or more Savont output directories into a SummarizedExperiment.

    Input layout (Savont)
    ---------------------
    Savont output is 100% the same fixed directory structure for every sample so
    the sample name must be supplied explicitly via ``sample_names``.

    Each Savont sample directory contains, after the initial (denoise/cluster)
    run:
      - ``feature-table.tsv``   (ASV read counts — used here)
      - ``final_asvs.fasta``    (representative ASV sequences)
      - ``final_clusters.tsv``  (cluster membership)
    and, after the classify run, additional taxonomy outputs including:
      - ``asv_mappings.tsv``    (per-ASV taxonomy assignments — used here)
      - ``genus_abundance.tsv``  (relative abundance at the genus level)
      - ``species_abundance.tsv`` (relative abundance at the species level)

    Only ``feature-table.tsv`` and ``asv_mappings.tsv`` are consumed today. The
    other outputs (``final_asvs.fasta``, ``final_clusters.tsv``, and any further
    classify-stage files) are intentionally noted here because we will need them
    later for more detailed reports.

    Parameters
    ----------
    input_dirs : str, Path, or list thereof
        Path(s) to Savont output directories. Each directory represents one
        sample and must contain ``feature-table.tsv`` and ``asv_mappings.tsv``.

    sample_names : list of str
        Sample names corresponding to each directory.

    Returns
    -------
    SummarizedExperiment
        Container with:
        - assays["abundance"]: relative abundance matrix (features × samples)
        - assays["counts"]: estimated counts matrix (features × samples)
        - row_data: taxonomy annotations per feature (indexed by tax_id)
        - col_data: sample metadata (indexed by sample name)
        - metadata: {"source": "savont"}
    """
    SAVONT_TAXONOMY_COLUMNS = (
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom",
    )

    # Normalise to list
    if isinstance(input_dirs, (str, Path)):
        input_dirs = [input_dirs]
    input_dirs = [Path(d) for d in input_dirs]

    if len(sample_names) != len(input_dirs):
        raise ValueError(
            f"Length mismatch: {len(input_dirs)} directories but "
            f"{len(sample_names)} sample names."
        )

    # Collect per-sample data
    abundance_series: dict[str, pd.Series] = {}
    counts_series: dict[str, pd.Series] = {}
    taxonomy_frames: list[pd.DataFrame] = []
    # Per-sample species-calling confidence (keyed by sample -> {tax_id: pct}),
    # stored as non-core key-value side data in metadata. The Report composite
    # consumes this; the abundance table itself stays confidence-free.
    species_confidence: dict[str, dict[str, float]] = {}

    for dir_path, name in zip(input_dirs, sample_names):
        reader = SavontReader(dir_path)
        df = reader.abundance().set_index("tax_id")

        abundance_series[name] = df["abundance"]
        counts_series[name] = df["estimated counts"].astype(float)
        species_confidence[name] = reader.species_confidence()

        tax_cols = [c for c in SAVONT_TAXONOMY_COLUMNS if c in df.columns]
        taxonomy_frames.append(df[tax_cols])

    # Build assay matrices (features × samples), filling missing features with 0
    abundance_matrix = pd.DataFrame(abundance_series).fillna(0.0)
    counts_matrix = pd.DataFrame(counts_series).fillna(0.0)

    # Build row_data from merged taxonomy (take first non-NA per tax_id)
    row_data = pd.concat(taxonomy_frames).groupby(level=0).first()
    row_data = row_data.reindex(abundance_matrix.index)

    # Build col_data
    col_data = pd.DataFrame(
        {
            "sample_name": sample_names,
            "source_file": [str(d) for d in input_dirs],
        },
        index=sample_names,
    )

    return SummarizedExperiment(
        assays={"abundance": abundance_matrix, "counts": counts_matrix},
        row_data=row_data,
        col_data=col_data,
        metadata={"source": "savont", "species_confidence": species_confidence},
    )
