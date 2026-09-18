"""SavontReader — an interface to a single Savont output directory.

Holds a reference to one Savont folder, knows its fixed file layout, reads each
file once (cached), and exposes derived products (species-level abundance table,
per-species calling confidence) so higher layers don't touch files directly.

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
- Multiple reference hits per ASV: take the **first hit** per ``asv_header``.
- Multiple ASVs collapsing onto one species: depths are **summed** per
  ``tax_id``; confidence takes the **max** ``alignment_identity`` (best-matching
  ASV wins per species).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import pandas as pd

from hoshi.lib.experiment import (
    SummarizedExperiment,
    aggregate_to_species,
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

_ABUNDANCE_OUTPUT_COLUMNS = (
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
)


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
        """One row per ASV present in the feature table, first taxonomy hit.

        Carries ``depth`` and ``alignment_identity`` (NaN when the mapping file
        predates that column) alongside the taxonomy columns.
        """
        mapping_df = self.asv_mappings
        depths = self.asv_depths
        has_identity = "alignment_identity" in mapping_df.columns

        # First hit per ASV (mapping files can list several references per ASV).
        first_hits = mapping_df.drop_duplicates(subset=["asv_header"], keep="first")

        rows = []
        for _, hit in first_hits.iterrows():
            asv_id = hit["asv_header"]
            if asv_id not in depths.index:
                continue
            row = {
                "asv_header": asv_id,
                "depth": depths[asv_id],
                "alignment_identity": (
                    pd.to_numeric(hit.get("alignment_identity"), errors="coerce")
                    if has_identity
                    else pd.NA
                ),
            }
            for col in _TAXONOMY_COLUMNS:
                row[col] = hit.get(col, pd.NA)
            rows.append(row)

        if not rows:
            raise ValueError("No ASVs could be resolved from the mapping file.")

        resolved = pd.DataFrame(rows)
        # tax_id is a nullable annotation now, not a join key: parse tolerantly
        # so ASVs with a missing/blank/non-numeric tax_id survive as <NA>.
        tax_num = pd.to_numeric(resolved["tax_id"], errors="coerce").astype("Int64")
        resolved["tax_id"] = tax_num.astype(str).where(tax_num.notna(), pd.NA)
        return resolved

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

        columns = [
            "feature_id",
            "abundance",
            "estimated counts",
            *self._SE_TAXONOMY_COLUMNS,
            "tax_id",
            "alignment_identity",
        ]
        for col in columns:
            if col not in resolved.columns:
                resolved[col] = pd.NA

        return (
            resolved[columns]
            .sort_values("abundance", ascending=False)
            .reset_index(drop=True)
        )

    # ─── Derived: species-level abundance table (rollup) ─────────────

    def species_abundance(self) -> pd.DataFrame:
        """Species-level abundance table (one row per ``tax_id``).

        A convenience rollup of :meth:`per_feature_abundance` that aggregates
        features to species by ``tax_id`` (features with no tax_id stay as their
        own rows, so counts still total correctly). Kept for callers that want a
        species view without going through a full experiment.

        Columns: abundance, tax_id, species, genus, family, order, class,
        phylum, superkingdom, estimated counts. Sorted by abundance descending.
        """
        per_feature = self.per_feature_abundance()
        rolled = aggregate_to_species(per_feature)

        for col in _ABUNDANCE_OUTPUT_COLUMNS:
            if col not in rolled.columns:
                rolled[col] = pd.NA
        return rolled[list(_ABUNDANCE_OUTPUT_COLUMNS)].reset_index(drop=True)

    # ─── Derived: per-feature calling confidence ─────────────────────

    def feature_confidence(self) -> dict[str, float]:
        """Per-feature calling confidence, keyed by ``asv_header`` (percent 0–100).

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
        step. Per-OTU calling confidence (:meth:`feature_confidence`, keyed by
        the raw feature id) is carried in ``metadata["confidence"][sample_name]``.

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
            - metadata: {"source": "savont", "confidence": {...}}
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

        # Per-feature confidence, keyed by the raw feature id.
        confidence = dict(self.feature_confidence())

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
                "confidence": {name: confidence},
            },
        )
