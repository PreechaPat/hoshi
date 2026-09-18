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
            df["estimated counts"] = pd.to_numeric(
                df["estimated counts"], errors="coerce"
            ).astype(int)
        else:
            df["estimated counts"] = pd.NA

        for column in EMU_OUTPUT_COLUMNS:
            if column not in df.columns:
                df[column] = pd.NA

        # Drop every non-output column and fix their order in one step.
        df = read_input_table(df, required_columns=list(EMU_OUTPUT_COLUMNS), keep_only_required=True)

    return df


def read_emu_abundance_into_summarizedexperiment(
    input_path: str | Path,
    *,
    sample_name: str | None = None,
    superkingdom: str | None = "Bacteria",
) -> SummarizedExperiment:
    """
    Load a single Emu rel-abundance TSV file into a SummarizedExperiment.

    Input layout (EMU)
    ------------------
    EMU emits one TSV file per sample, and the *filename itself* encodes the
    sample name as a prefix (e.g. ``sample01_rel-abundance.tsv``,
    ``barcode11.fastq_rel-abundance.tsv``). There is no fixed parent-directory
    structure — files may live anywhere. Because of this, when ``sample_name``
    is not supplied we derive it from the filename stem (see
    ``_derive_sample_name``, which strips the ``_rel-abundance`` / ``.fastq``
    suffixes). This is the key structural difference from Savont, which instead
    uses a fixed per-sample directory layout (see
    ``read_savont_abundance_into_summarizedexperiment``).

    This reads exactly one sample and returns a single-sample experiment.
    Combining several samples into one experiment (and the sample-scoping that
    entails) is a separate, explicit step — see
    :func:`hoshi.lib.experiment.combine_experiments`.

    Parameters
    ----------
    input_path : str or Path
        Path to an Emu rel-abundance TSV file (one sample).

    sample_name : str, optional
        Sample name for the resulting experiment. If None, it is derived from
        the filename (stem without '_rel-abundance' suffix).

    superkingdom : str or None, default "Bacteria"
        Value to fill in the superkingdom column when Emu leaves it empty.
        Emu 16S output typically does not populate superkingdom, so this
        override ensures downstream formats (e.g., Kraken2) get a proper
        domain line. Set to None to skip filling.

    Returns
    -------
    SummarizedExperiment
        Single-sample container with:
        - assays["abundance"]: relative abundance matrix (OTUs × 1)
        - assays["counts"]: estimated counts matrix (OTUs × 1)
        - row_data: tax_id + taxonomy annotations per OTU (raw feature id index)
        - col_data: sample metadata (indexed by sample name)
        - metadata: {"source": "emu"}

    Notes
    -----
    EMU has no ASV/OTU concept — each row is a ``tax_id``. To align EMU with the
    per-OTU model used for Savont, **one tax_id/species is treated as one OTU**
    (no sequence identity). The raw feature id therefore *is* the ``tax_id``;
    ``tax_id`` is also carried as a nullable row_data column.
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

    path = Path(input_path)
    name = sample_name if sample_name is not None else _derive_sample_name(path)

    df = read_emu_abundance(str(path), reorder=True)

    # Filter out unmapped/unclassified control rows.
    df = df[pd.to_numeric(df["tax_id"], errors="coerce").notna()].copy()
    df["tax_id"] = df["tax_id"].astype(int).astype(str)

    # Raw feature id = tax_id (one tax_id == one feature for EMU).
    df.index = pd.Index(df["tax_id"].astype(str), name="feature_id")

    abundance_matrix = pd.DataFrame({name: df["abundance"]}).fillna(0.0)
    counts_matrix = pd.DataFrame(
        {name: df["estimated counts"].astype(float)}
    ).fillna(0.0)

    row_cols = ["tax_id", *[c for c in EMU_TAXONOMY_COLUMNS if c in df.columns]]
    row_data = df[row_cols].reindex(abundance_matrix.index)

    # Fill empty superkingdom values with the override (EMU 16S usually omits it).
    if superkingdom and "superkingdom" in row_data.columns:
        row_data["superkingdom"] = row_data["superkingdom"].fillna(superkingdom)
        row_data["superkingdom"] = row_data["superkingdom"].replace("", superkingdom)

    col_data = pd.DataFrame(
        {"sample_name": [name], "source_file": [str(path)]},
        index=[name],
    )

    return SummarizedExperiment(
        assays={"abundance": abundance_matrix, "counts": counts_matrix},
        row_data=row_data,
        col_data=col_data,
        metadata={"source": "emu"},
    )


def _derive_sample_name(path: Path) -> str:
    """Derive a sample name from an Emu output filename.

    Canonical implementation lives in :func:`hoshi.lib.emu_reader._derive_sample_name`
    (the EMU-specific naming rule belongs with the EMU reader). Re-exported here
    for backward compatibility with ``from hoshi.lib.ingress import _derive_sample_name``.
    """
    from hoshi.lib.emu_reader import _derive_sample_name as _derive  # noqa: PLC0415

    return _derive(path)
