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

    Raises
    ------
    ValueError
        If ``directory`` does not exist or is not a directory.
    """

    directory: Path

    def __init__(self, directory: str | Path) -> None:
        path = Path(directory)
        if not path.is_dir():
            raise ValueError(f"Savont directory not found: {path}")
        # frozen dataclass: set via object.__setattr__
        object.__setattr__(self, "directory", path)

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
        path = self.feature_table_path
        if not path.is_file():
            raise ValueError(f"Feature table not found: {path}")
        df = pd.read_csv(path, sep="\t")
        if "#OTU ID" not in df.columns:
            raise ValueError(
                f"feature-table.tsv missing '#OTU ID' column. "
                f"Found: {list(df.columns)}"
            )
        return df.set_index("#OTU ID")

    @cached_property
    def asv_mappings(self) -> pd.DataFrame:
        """Raw ``asv_mappings.tsv`` (per-ASV taxonomy assignments)."""
        path = self.asv_mappings_path
        if not path.is_file():
            raise ValueError(f"ASV mapping not found: {path}")
        df = pd.read_csv(path, sep="\t")
        required = ["asv_header", "tax_id", "species"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"asv_mappings.tsv missing columns: {missing}")
        return df

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
        resolved["tax_id"] = resolved["tax_id"].astype(int).astype(str)
        return resolved

    # ─── Derived: species-level abundance table ──────────────────────

    def abundance(self) -> pd.DataFrame:
        """Species-level abundance table (one row per ``tax_id``).

        Columns: abundance, tax_id, species, genus, family, order, class,
        phylum, superkingdom, estimated counts. Sorted by abundance desc, with a
        positional integer index (``tax_id`` is a regular column).
        """
        resolved = self._resolved_asvs

        agg = (
            resolved.groupby("tax_id", as_index=False).agg(
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

        total_depth = agg["depth"].sum()
        agg["abundance"] = agg["depth"] / total_depth if total_depth > 0 else 0.0
        agg = agg.rename(columns={"depth": "estimated counts"})

        for col in _ABUNDANCE_OUTPUT_COLUMNS:
            if col not in agg.columns:
                agg[col] = pd.NA

        return (
            agg[list(_ABUNDANCE_OUTPUT_COLUMNS)]
            .sort_values("abundance", ascending=False)
            .reset_index(drop=True)
        )

    # ─── Derived: species-calling confidence ─────────────────────────

    def species_confidence(self) -> dict[str, float]:
        """Per-species calling confidence, keyed by ``tax_id`` (percent 0–100).

        For each species, this is the **max** ``alignment_identity`` across the
        ASVs that resolve to it — the identity of the best-matching ASV. Species
        with no identity data (all-NaN) are omitted.
        """
        resolved = self._resolved_asvs
        if "alignment_identity" not in resolved.columns:
            return {}

        identity = pd.to_numeric(resolved["alignment_identity"], errors="coerce")
        per_species = (
            resolved.assign(_identity=identity)
            .groupby("tax_id")["_identity"]
            .max()
            .dropna()
        )
        return {str(tax_id): float(value) for tax_id, value in per_species.items()}
