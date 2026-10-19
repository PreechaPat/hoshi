"""SavontReader — an interface to a single Savont output directory.

Holds a reference to one Savont folder, knows its fixed file layout, reads each
file once (cached), and exposes derived products (per-feature abundance table,
per-feature estimated sequence identity) so higher layers don't touch files
directly.
Species-level rollups are an explicit downstream step
(:func:`hoshi.lib.experiment.species_view`), not a reader concern.

Savont folder layout (fixed for every sample)
---------------------------------------------
Initial (denoise/cluster) run::

    feature-table.tsv    # ASV read counts            (used)
    final_asvs.fasta     # representative ASV sequences (not yet used)
    final_clusters.tsv   # cluster membership           (not yet used)

Classify run adds::

    asv_mappings.tsv        # per-ASV taxonomy assignments   (used)
    genus_abundance.tsv     # relative abundance, genus level (not yet used)
    species_abundance.tsv   # relative abundance, species lvl (not yet used)

Quirks handled here
-------------------
Resolution mirrors ``savont export`` so counts per lineage match it:

- Multiple reference hits per ASV: take the **first row** per ``asv_header``
  (the file is sorted by Savont's EM abundance; on an exact EM tie Savont's own
  ``species_abundance.tsv`` picks the last tied row instead — see
  ``docs/research/savont-hit-selection.md``).
- Hits below the species threshold (``species`` = ``UNCLASSIFIED``) keep the
  nearest reference's tax_id with an ``_UNK<nnn>`` tag, so they never merge
  with, or are reported as, that species.
- ASVs with no alignment are missing from ``asv_mappings.tsv``; they are kept
  as ``unmapped`` so relative abundance uses every read in the feature table.
- Multiple ASVs collapsing onto one species: depths are **summed** per
  ``tax_id``; sequence identity takes the **max** ``alignment_identity``
  (best-matching ASV wins per species).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import pandas as pd

from hoshi.lib.experiment import (
    SummarizedExperiment,
)

# Fixed filenames within a Savont sample directory.
_FEATURE_TABLE = "feature-table.tsv"
_ASV_MAPPINGS = "asv_mappings.tsv"
_FINAL_ASVS = "final_asvs.fasta"
_FINAL_CLUSTERS = "final_clusters.tsv"
_GENUS_ABUNDANCE = "genus_abundance.tsv"
_SPECIES_ABUNDANCE = "species_abundance.tsv"

_TAXONOMY_COLUMNS = (
    "tax_id",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
)

# Savont writes this (or ``UNCLASSIFIED-(<asv>)`` with --detailed-unclassified)
# for every rank the hit's identity is too low to call.
_UNCLASSIFIED_MARKER = "UNCLASSIFIED"
_RANKS_ABOVE_SPECIES = ("genus", "family", "order", "class", "phylum", "superkingdom")

# tax_id for ASVs with no alignment: counted in the total, never shown as an
# organism (it is one of META_TAX_IDS, which species_view drops).
UNMAPPED_TAX_ID = "unmapped"


def _is_unclassified(value: object) -> bool:
    return isinstance(value, str) and value.startswith(_UNCLASSIFIED_MARKER)


def _clean_tax_id(value: object) -> str | None:
    """Normalise a raw tax_id cell to a digit string; blank/non-numeric -> None."""
    number = pd.to_numeric(value, errors="coerce")
    return None if pd.isna(number) else str(int(number))


def _tag_unresolved_tax_id(tax_id: str | None, hit: pd.Series) -> str | None:
    """Tag the nearest-reference tax_id of a hit Savont could not call to species.

    Savont keeps the nearest reference's tax_id even when the species (and maybe
    higher ranks) are ``UNCLASSIFIED``. Using it bare would count those reads as
    that species and match it in the pathogen sheet. The tag
    ``<tax_id>_UNK<nnn>`` keeps the nearest reference traceable while making it
    a different key; ``nnn`` is how many ranks above species are also
    unresolved (``000`` = only species), so ASVs resolved to different depths
    never merge under one lineage.
    """
    if tax_id is None or not _is_unclassified(hit.get("species")):
        return tax_id
    unresolved_above = sum(_is_unclassified(hit.get(rank)) for rank in _RANKS_ABOVE_SPECIES)
    return f"{tax_id}_UNK{unresolved_above:03d}"


@dataclass(frozen=True)
class SavontReader:
    """Interface to one Savont output directory.

    Parameters
    ----------
    directory : str or Path
        Path to a single Savont output directory.
    sample_name : str, optional
        Sample name used when producing a :class:`SummarizedExperiment`.
        Defaults to the directory name.

    Raises
    ------
    ValueError
        If ``directory`` does not exist or is not a directory.
    """

    directory: Path
    sample_name: str

    def __init__(self, directory: str | Path, *, sample_name: str | None = None) -> None:
        path = Path(directory)
        if not path.is_dir():
            raise ValueError(f"Savont directory not found: {path}")
        # frozen dataclass: set via object.__setattr__
        object.__setattr__(self, "directory", path)
        object.__setattr__(self, "sample_name", sample_name or path.name)

    # ─── Fixed-layout path properties ────────────────────────────────

    @property
    def feature_table_path(self) -> Path:
        return self.directory / _FEATURE_TABLE

    @property
    def asv_mappings_path(self) -> Path:
        return self.directory / _ASV_MAPPINGS

    @property
    def final_asvs_path(self) -> Path:
        return self.directory / _FINAL_ASVS

    @property
    def final_clusters_path(self) -> Path:
        return self.directory / _FINAL_CLUSTERS

    @property
    def genus_abundance_path(self) -> Path:
        return self.directory / _GENUS_ABUNDANCE

    @property
    def species_abundance_path(self) -> Path:
        return self.directory / _SPECIES_ABUNDANCE

    # ─── Raw file readers (cached: read each file once) ──────────────

    @cached_property
    def feature_table(self) -> pd.DataFrame:
        """Raw ``feature-table.tsv`` indexed by ``#OTU ID`` (ASV id)."""
        from hoshi.lib.ingress import read_input_table  # noqa: PLC0415

        path = self.feature_table_path
        if not path.is_file():
            raise ValueError(f"Feature table not found: {path}")
        # Sample count columns are dynamic (one per sample), so we validate the
        # index column only and keep everything else.
        df = read_input_table(str(path), required_columns=["#OTU ID"], sep="\t")
        return df.set_index("#OTU ID")

    @cached_property
    def asv_mappings(self) -> pd.DataFrame:
        """Raw ``asv_mappings.tsv`` (per-ASV taxonomy assignments)."""
        from hoshi.lib.ingress import read_input_table  # noqa: PLC0415

        path = self.asv_mappings_path
        if not path.is_file():
            raise ValueError(f"ASV mapping not found: {path}")
        return read_input_table(
            str(path),
            required_columns=["asv_header", "tax_id", "species"],
            sep="\t",
        )

    # ─── Derived: per-ASV depth ──────────────────────────────────────

    @cached_property
    def asv_depths(self) -> pd.Series:
        """Read depth per ASV, summed across all sample columns."""
        feature_df = self.feature_table
        sample_cols = list(feature_df.columns)
        if not sample_cols:
            raise ValueError("feature-table.tsv has no sample columns.")
        depths = feature_df[sample_cols].sum(axis=1).astype(int)
        depths.name = "depth"
        return depths

    # ─── Derived: resolved per-ASV table (first hit per ASV) ─────────

    @cached_property
    def _resolved_asvs(self) -> pd.DataFrame:
        """One row per ASV in the feature table, resolved like ``savont export``.

        Mapped ASVs take their first ``asv_mappings.tsv`` row. Below-threshold
        hits get a tagged ``tax_id`` (see :func:`_tag_unresolved_tax_id`). ASVs
        missing from the mapping file (no alignment) are kept with the
        ``unmapped`` meta tax_id and no taxonomy, so their reads stay in the
        total. Carries ``depth`` and ``alignment_identity`` (NaN when absent).
        """
        mapping_df = self.asv_mappings
        depths = self.asv_depths
        has_identity = "alignment_identity" in mapping_df.columns

        # First hit per ASV (mapping files can list several references per ASV).
        first_hits = mapping_df.drop_duplicates(subset=["asv_header"], keep="first")
        first_hits = first_hits.set_index("asv_header")

        rows = []
        for asv_id, depth in depths.items():
            row = {"asv_header": asv_id, "depth": depth, "alignment_identity": pd.NA}
            if asv_id not in first_hits.index:
                row.update({col: pd.NA for col in _TAXONOMY_COLUMNS})
                row["tax_id"] = UNMAPPED_TAX_ID
                rows.append(row)
                continue

            hit = first_hits.loc[asv_id]
            if has_identity:
                row["alignment_identity"] = pd.to_numeric(
                    hit.get("alignment_identity"), errors="coerce"
                )
            for col in _TAXONOMY_COLUMNS:
                row[col] = hit.get(col, pd.NA)
            row["tax_id"] = _tag_unresolved_tax_id(_clean_tax_id(hit.get("tax_id")), hit)
            rows.append(row)

        if not rows:
            raise ValueError(f"No ASVs found in {self.feature_table_path}.")

        return pd.DataFrame(rows)

    # ─── Derived: per-OTU (per-ASV) abundance table ──────────────────

    def per_feature_abundance(self) -> pd.DataFrame:
        """Per-feature (per-OTU/ASV) abundance table (one row per ``asv_header``).

        This is the feature-level table: each ASV keeps its own read depth,
        relative abundance (depth / total depth), taxonomy assignment, tax_id
        (nullable), and ``alignment_identity``. No aggregation to species is done
        here — that is an explicit downstream step (see
        :func:`hoshi.lib.experiment.aggregate_to_species`).

        Columns: feature_id, abundance, estimated counts, tax_id, species,
        genus, family, order, class, phylum, superkingdom, alignment_identity.
        Sorted by abundance descending. ``feature_id`` here is Savont's own
        ``asv_header`` (the reader scopes it per sample when building an
        experiment).
        """
        resolved = self._resolved_asvs.copy()
        total_depth = resolved["depth"].sum()
        resolved["abundance"] = (
            resolved["depth"] / total_depth if total_depth > 0 else 0.0
        )
        resolved = resolved.rename(
            columns={"depth": "estimated counts", "asv_header": "feature_id"}
        )

        # Every column below is guaranteed present: _resolved_asvs always builds
        # depth (→ estimated counts), asv_header (→ feature_id), tax_id, the
        # taxonomy ranks, and alignment_identity; abundance is added just above.
        columns = [
            "feature_id",
            "abundance",
            "estimated counts",
            *self._SE_TAXONOMY_COLUMNS,
            "tax_id",
            "alignment_identity",
        ]

        return (
            resolved[columns]
            .sort_values("abundance", ascending=False)
            .reset_index(drop=True)
        )

    # ─── Derived: per-feature estimated sequence identity ────────────

    def feature_sequence_identity(self) -> dict[str, float]:
        """Per-feature sequence identity, keyed by ``asv_header`` (percent 0–100).

        This is each ASV's ``alignment_identity`` — a genuine per-feature value.
        Features with no identity data are omitted.
        """
        resolved = self._resolved_asvs
        if "alignment_identity" not in resolved.columns:
            return {}
        identity = pd.to_numeric(resolved["alignment_identity"], errors="coerce")
        return {
            str(asv): float(val)
            for asv, val in zip(resolved["asv_header"], identity)
            if pd.notna(val)
        }

    # ─── Unified reader interface ────────────────────────────────────

    _SE_TAXONOMY_COLUMNS = (
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom",
    )

    def to_summarized_experiment(self) -> SummarizedExperiment:
        """Build a single-sample per-OTU :class:`SummarizedExperiment`.

        Features are OTUs/ASVs — one row per ``asv_header`` — keyed by Savont's
        **raw** feature id (unscoped). Taxonomy and ``tax_id`` are nullable
        ``row_data`` annotations; species aggregation is an explicit downstream
        step. Per-OTU estimated sequence identity (:meth:`feature_sequence_identity`,
        keyed by the raw feature id) is carried in
        ``metadata["sequence_identity"][sample_name]``.

        Feature ids are scoped per sample only when several samples are merged
        (see :func:`hoshi.lib.experiment.combine_experiments`), so a single
        sample keeps its classifier-native ids here.

        Returns
        -------
        SummarizedExperiment
            Container with:
            - assays["abundance"]: relative abundance matrix (OTUs × 1)
            - assays["counts"]: estimated counts matrix (OTUs × 1)
            - row_data: tax_id + taxonomy annotations per OTU (indexed by the
              raw OTU id)
            - col_data: sample metadata (indexed by sample name)
            - metadata: {"source": "savont", "sequence_identity": {...}}
        """
        name = self.sample_name
        df = self.per_feature_abundance().copy()

        # Raw (unscoped) feature id — Savont's own ASV id.
        df.index = pd.Index(df["feature_id"].astype(str), name="feature_id")

        abundance_matrix = pd.DataFrame({name: df["abundance"]}).fillna(0.0)
        counts_matrix = pd.DataFrame(
            {name: df["estimated counts"].astype(float)}
        ).fillna(0.0)

        row_cols = [c for c in ("tax_id", *self._SE_TAXONOMY_COLUMNS) if c in df.columns]
        row_data = df[row_cols].reindex(abundance_matrix.index)

        # Per-feature sequence identity, keyed by the raw feature id.
        sequence_identity = dict(self.feature_sequence_identity())

        col_data = pd.DataFrame(
            {
                "sample_name": [name],
                "source_file": [str(self.directory)],
            },
            index=[name],
        )

        return SummarizedExperiment(
            assays={"abundance": abundance_matrix, "counts": counts_matrix},
            row_data=row_data,
            col_data=col_data,
            metadata={
                "source": "savont",
                "sequence_identity": {name: sequence_identity},
            },
        )