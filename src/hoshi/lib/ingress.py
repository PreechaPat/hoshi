from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from hoshi.lib.experiment import SummarizedExperiment


def read_input_table(input_data, required_columns=None, sep=None):
    """
    Read a TSV/CSV file or a pandas DataFrame and validate required columns.

    Parameters
    ----------
    input_data : str or pd.DataFrame
        Input file path or already-loaded DataFrame.

    required_columns : list of str, optional
        List of column names that must exist in the table.
        If any are missing, raises ValueError.

    sep : str, optional
        Delimiter for file reading (e.g., '\t' for TSV). If None, will auto-detect.

    Returns
    -------
    pd.DataFrame
        The loaded and validated DataFrame.

    Raises
    ------
    ValueError
        If required columns are missing or file cannot be parsed.
    """
    if isinstance(input_data, pd.DataFrame):
        df = input_data.copy()
    elif isinstance(input_data, str):
        if not os.path.isfile(input_data):
            raise ValueError(f"File not found: {input_data}")
        try:
            df = pd.read_csv(input_data, sep=sep, engine="python")  # auto-detect by default
        except Exception as e:
            raise ValueError(f"Failed to read file '{input_data}': {e}")
    else:
        raise ValueError("input_data must be a file path or a pandas DataFrame.")

    if required_columns:
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

    return df


def read_emu_abundance(input_data, *, sep: str | None = "\t", reorder: bool = False) -> pd.DataFrame:
    """
    Load an EMU rel-abundance table and normalise its column layout.
    """

    EMU_REQUIRED_COLUMNS = ("tax_id", "abundance")
    EMU_TAXONOMY_COLUMNS = (
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    )
    EMU_OUTPUT_COLUMNS = (
        "tax_id",
        "abundance",
        "estimated counts",
        *EMU_TAXONOMY_COLUMNS,
    )
    EMU_NUMERIC_COLUMNS = ("abundance", "estimated counts")
    df = read_input_table(input_data, required_columns=list(EMU_REQUIRED_COLUMNS), sep=sep)

    # Only pick up those superkingdom -> species taxonomy
    for column in EMU_TAXONOMY_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA

    for column in EMU_NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    if reorder:
        if "estimated counts" in df.columns:
            df = df.rename(columns={"estimated counts": "estimated counts"})
            df["estimated counts"] = pd.to_numeric(df["estimated counts"], errors="coerce").astype(int)
        if "estimated counts" not in df.columns:
            df["estimated counts"] = pd.NA

        for column in EMU_OUTPUT_COLUMNS:
            if column not in df.columns:
                df[column] = pd.NA

        df = df.loc[:, list(EMU_OUTPUT_COLUMNS)]

    return df


def emu_to_experiment(
    input_paths: str | Path | list[str | Path],
    *,
    sample_names: list[str] | None = None,
    superkingdom: str | None = "Bacteria",
) -> SummarizedExperiment:
    """
    Load one or more Emu rel-abundance TSV files into a SummarizedExperiment.

    Parameters
    ----------
    input_paths : str, Path, or list thereof
        Path(s) to Emu rel-abundance TSV files. Each file represents one sample.

    sample_names : list of str, optional
        Sample names corresponding to each file. If None, names are derived
        from the filename (stem without '_rel-abundance' suffix).

    superkingdom : str or None, default "Bacteria"
        Value to fill in the superkingdom column when Emu leaves it empty.
        Emu 16S output typically does not populate superkingdom, so this
        override ensures downstream formats (e.g., Kraken2) get a proper
        domain line. Set to None to skip filling.

    Returns
    -------
    SummarizedExperiment
        Container with:
        - assays["abundance"]: relative abundance matrix (features × samples)
        - assays["counts"]: estimated counts matrix (features × samples)
        - row_data: taxonomy annotations per feature (indexed by tax_id)
        - col_data: sample metadata (indexed by sample name)
        - metadata: {"source": "emu"}
    """
    EMU_TAXONOMY_COLUMNS = (
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    )

    # Normalise to list
    if isinstance(input_paths, (str, Path)):
        input_paths = [input_paths]
    input_paths = [Path(p) for p in input_paths]

    if sample_names is None:
        sample_names = [_derive_sample_name(p) for p in input_paths]

    if len(sample_names) != len(input_paths):
        raise ValueError(
            f"Length mismatch: {len(input_paths)} files but "
            f"{len(sample_names)} sample names."
        )

    # Collect per-sample data
    abundance_series: dict[str, pd.Series] = {}
    counts_series: dict[str, pd.Series] = {}
    taxonomy_frames: list[pd.DataFrame] = []

    for path, name in zip(input_paths, sample_names):
        df = read_emu_abundance(str(path), reorder=True)

        # Filter out unmapped/unclassified control rows
        df = df[pd.to_numeric(df["tax_id"], errors="coerce").notna()].copy()
        df["tax_id"] = df["tax_id"].astype(int).astype(str)
        df = df.set_index("tax_id")

        abundance_series[name] = df["abundance"]
        counts_series[name] = df["estimated counts"].astype(float)

        # Collect taxonomy (will merge/deduplicate later)
        tax_cols = [c for c in EMU_TAXONOMY_COLUMNS if c in df.columns]
        taxonomy_frames.append(df[tax_cols])

    # Build assay matrices (features × samples), filling missing features with 0
    abundance_matrix = pd.DataFrame(abundance_series).fillna(0.0)
    counts_matrix = pd.DataFrame(counts_series).fillna(0.0)

    # Build row_data from merged taxonomy (take first non-NA per tax_id)
    row_data = pd.concat(taxonomy_frames).groupby(level=0).first()
    row_data = row_data.reindex(abundance_matrix.index)

    # Fill empty superkingdom values with the override
    if superkingdom and "superkingdom" in row_data.columns:
        row_data["superkingdom"] = row_data["superkingdom"].fillna(superkingdom)
        row_data["superkingdom"] = row_data["superkingdom"].replace("", superkingdom)

    # Build col_data
    col_data = pd.DataFrame(
        {
            "sample_name": sample_names,
            "source_file": [str(p) for p in input_paths],
        },
        index=sample_names,
    )

    return SummarizedExperiment(
        assays={"abundance": abundance_matrix, "counts": counts_matrix},
        row_data=row_data,
        col_data=col_data,
        metadata={"source": "emu"},
    )


def read_savont_abundance(
    feature_table: str | Path,
    asv_mapping: str | Path,
) -> pd.DataFrame:
    """
    Load Savont feature table and ASV mapping to build a species-level abundance table.

    Reads the feature table (ASV read counts) and ASV mapping (taxonomy
    assignments per ASV), resolves each ASV to a single species (first hit),
    and aggregates read counts by species to compute relative abundance.

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
    SAVONT_TAXONOMY_COLUMNS = (
        "tax_id",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom",
    )

    feature_path = Path(feature_table)
    mapping_path = Path(asv_mapping)

    if not feature_path.is_file():
        raise ValueError(f"Feature table not found: {feature_path}")
    if not mapping_path.is_file():
        raise ValueError(f"ASV mapping not found: {mapping_path}")

    # Read feature table: columns are "#OTU ID" and one or more sample columns
    feature_df = pd.read_csv(feature_path, sep="\t")
    if "#OTU ID" not in feature_df.columns:
        raise ValueError(
            f"feature-table.tsv missing '#OTU ID' column. "
            f"Found: {list(feature_df.columns)}"
        )

    # Extract depth per ASV (sum across sample columns for multi-column case)
    sample_cols = [c for c in feature_df.columns if c != "#OTU ID"]
    if not sample_cols:
        raise ValueError("feature-table.tsv has no sample columns.")

    feature_df = feature_df.set_index("#OTU ID")
    asv_counts = feature_df[sample_cols].sum(axis=1).astype(int)
    asv_counts.name = "depth"

    # Read ASV mappings
    mapping_df = pd.read_csv(mapping_path, sep="\t")
    required_mapping_cols = ["asv_header", "tax_id", "species"]
    missing = [c for c in required_mapping_cols if c not in mapping_df.columns]
    if missing:
        raise ValueError(f"asv_mappings.tsv missing columns: {missing}")

    # Resolve each ASV to a single tax_id (first hit per ASV)
    first_hits = mapping_df.drop_duplicates(subset=["asv_header"], keep="first")

    resolved_rows = []
    for _, hit in first_hits.iterrows():
        asv_id = hit["asv_header"]
        if asv_id not in asv_counts.index:
            continue

        row = {"asv_header": asv_id, "depth": asv_counts[asv_id]}
        for col in SAVONT_TAXONOMY_COLUMNS:
            row[col] = hit.get(col, pd.NA)
        resolved_rows.append(row)

    if not resolved_rows:
        raise ValueError("No ASVs could be resolved from the mapping file.")

    resolved_df = pd.DataFrame(resolved_rows)
    resolved_df["tax_id"] = resolved_df["tax_id"].astype(int).astype(str)

    # Aggregate by tax_id: sum depths, keep first taxonomy
    agg_df = (
        resolved_df.groupby("tax_id", as_index=False)
        .agg(
            {
                "depth": "sum",
                "species": "first",
                "genus": "first",
                "family": "first",
                "order": "first",
                "class": "first",
                "phylum": "first",
                "superkingdom": "first",
            }
        )
    )

    # Compute relative abundance
    total_depth = agg_df["depth"].sum()
    agg_df["abundance"] = agg_df["depth"] / total_depth if total_depth > 0 else 0.0
    agg_df = agg_df.rename(columns={"depth": "estimated counts"})

    # Reorder columns
    output_columns = [
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
    for col in output_columns:
        if col not in agg_df.columns:
            agg_df[col] = pd.NA

    agg_df = agg_df[output_columns].sort_values("abundance", ascending=False).reset_index(drop=True)

    return agg_df


def savont_to_experiment(
    input_dirs: str | Path | list[str | Path],
    *,
    sample_names: list[str],
) -> SummarizedExperiment:
    """
    Load one or more Savont output directories into a SummarizedExperiment.

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

    for dir_path, name in zip(input_dirs, sample_names):
        df = read_savont_abundance(
            dir_path / "feature-table.tsv",
            dir_path / "asv_mappings.tsv",
        )
        df = df.set_index("tax_id")

        abundance_series[name] = df["abundance"]
        counts_series[name] = df["estimated counts"].astype(float)

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
        metadata={"source": "savont"},
    )


def _derive_sample_name(path: Path) -> str:
    """Derive a sample name from an Emu output filename."""
    stem = path.stem
    # Strip common Emu suffixes
    for suffix in ("_rel-abundance", ".fastq_rel-abundance", "_rel_abundance"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    # Also handle patterns like "barcode11.fastq_rel-abundance"
    if stem.endswith(".fastq"):
        stem = stem[: -len(".fastq")]
    return stem
