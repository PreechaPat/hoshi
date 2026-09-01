"""EMU ingress — read EMU rel-abundance output into a ``SummarizedExperiment``.

EMU folder / filename quirk
---------------------------
EMU emits **one TSV file per sample**, and the *filename itself* encodes the
sample name as a prefix, e.g.::

    sample01_rel-abundance.tsv
    barcode11.fastq_rel-abundance.tsv

There is **no fixed parent-directory structure** — the files can live anywhere
and are addressed individually. Consequently, when the caller does not pass
``sample_names`` explicitly we derive each sample name from the filename stem
(see :func:`_derive_sample_name`, which strips the ``_rel-abundance`` /
``.fastq`` suffixes).

This is the key structural difference from Savont (see ``ingress_savont``),
which uses a fixed per-sample *directory* layout where the filenames carry no
sample identity, so Savont requires ``sample_names`` to be supplied. Keeping the
two sources in separate modules lets each document its own quirk without the
files growing unwieldy.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.ingress import read_input_table


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


def read_emu_abundance_into_summarizedexperiment(
    input_paths: str | Path | list[str | Path],
    *,
    sample_names: list[str] | None = None,
    superkingdom: str | None = "Bacteria",
) -> SummarizedExperiment:
    """
    Load one or more Emu rel-abundance TSV files into a SummarizedExperiment.

    Input layout (EMU)
    ------------------
    EMU emits one TSV file per sample, and the *filename itself* encodes the
    sample name as a prefix (e.g. ``sample01_rel-abundance.tsv``,
    ``barcode11.fastq_rel-abundance.tsv``). There is no fixed parent-directory
    structure — files may live anywhere. Because of this, when ``sample_names``
    is not supplied we derive each sample name from the filename stem (see
    ``_derive_sample_name``, which strips the ``_rel-abundance`` / ``.fastq``
    suffixes). This is the key structural difference from Savont, which instead
    uses a fixed per-sample directory layout (see
    ``read_savont_abundance_into_summarizedexperiment``).

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
