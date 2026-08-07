"""
SummarizedExperiment-like container for abundance data.

Inspired by Bioconductor's SummarizedExperiment in R, this provides a
structured container that keeps assay matrices, feature (row) metadata,
and sample (column) metadata together.

Typical usage with Emu output:
    - assays: {"abundance": DataFrame, "counts": DataFrame}
      Each assay is features (rows) × samples (columns).
    - row_data: taxonomy annotations per feature (tax_id, species, genus, ...)
    - col_data: sample-level metadata (sample_name, source_file, ...)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd


@dataclass
class SummarizedExperiment:
    """
    A lightweight SummarizedExperiment container.

    Attributes
    ----------
    assays : dict[str, pd.DataFrame]
        Named assay matrices. Each DataFrame has features as rows and
        samples as columns. All assays must share the same row and
        column indices.

    row_data : pd.DataFrame
        Feature-level metadata. Index must match the row index of the
        assay matrices (e.g., tax_id).

    col_data : pd.DataFrame
        Sample-level metadata. Index must match the column index of
        the assay matrices (e.g., sample names).

    metadata : dict[str, Any]
        Experiment-level metadata (e.g., tool version, database used).
    """

    assays: dict[str, pd.DataFrame] = field(default_factory=dict)
    row_data: pd.DataFrame = field(default_factory=lambda: pd.DataFrame())
    col_data: pd.DataFrame = field(default_factory=lambda: pd.DataFrame())
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._validate()

    def _validate(self) -> None:
        """Check that dimensions are consistent across assays and metadata."""
        if not self.assays:
            return

        # All assays must have the same shape
        reference_name = next(iter(self.assays))
        reference = self.assays[reference_name]
        n_features, n_samples = reference.shape

        for name, assay in self.assays.items():
            if assay.shape != (n_features, n_samples):
                raise ValueError(
                    f"Assay '{name}' has shape {assay.shape}, expected "
                    f"{(n_features, n_samples)} (from '{reference_name}')."
                )
            if not assay.index.equals(reference.index):
                raise ValueError(
                    f"Assay '{name}' row index does not match '{reference_name}'."
                )
            if not assay.columns.equals(reference.columns):
                raise ValueError(
                    f"Assay '{name}' column index does not match '{reference_name}'."
                )

        # row_data index must match assay row index (if row_data is non-empty)
        if not self.row_data.empty:
            if not self.row_data.index.equals(reference.index):
                raise ValueError(
                    "row_data index does not match assay row index."
                )

        # col_data index must match assay column index (if col_data is non-empty)
        if not self.col_data.empty:
            if not self.col_data.index.equals(reference.columns):
                raise ValueError(
                    "col_data index does not match assay column index."
                )

    # ─── Convenience properties ──────────────────────────────────────

    @property
    def n_features(self) -> int:
        """Number of features (rows)."""
        if not self.assays:
            return 0
        return next(iter(self.assays.values())).shape[0]

    @property
    def n_samples(self) -> int:
        """Number of samples (columns)."""
        if not self.assays:
            return 0
        return next(iter(self.assays.values())).shape[1]

    @property
    def feature_ids(self) -> pd.Index:
        """Row index (feature identifiers)."""
        if not self.assays:
            return pd.Index([])
        return next(iter(self.assays.values())).index

    @property
    def sample_ids(self) -> pd.Index:
        """Column index (sample identifiers)."""
        if not self.assays:
            return pd.Index([])
        return next(iter(self.assays.values())).columns

    @property
    def assay_names(self) -> list[str]:
        """Names of stored assays."""
        return list(self.assays.keys())

    @property
    def shape(self) -> tuple[int, int]:
        """(n_features, n_samples)."""
        return (self.n_features, self.n_samples)

    # ─── Subset operations ───────────────────────────────────────────

    def subset_features(self, feature_ids: list | pd.Index) -> SummarizedExperiment:
        """Return a new SummarizedExperiment with only the specified features."""
        new_assays = {name: assay.loc[feature_ids] for name, assay in self.assays.items()}
        new_row_data = self.row_data.loc[feature_ids] if not self.row_data.empty else self.row_data
        return SummarizedExperiment(
            assays=new_assays,
            row_data=new_row_data,
            col_data=self.col_data.copy(),
            metadata=self.metadata.copy(),
        )

    def subset_samples(self, sample_ids: list | pd.Index) -> SummarizedExperiment:
        """Return a new SummarizedExperiment with only the specified samples."""
        new_assays = {name: assay[sample_ids] for name, assay in self.assays.items()}
        new_col_data = self.col_data.loc[sample_ids] if not self.col_data.empty else self.col_data
        return SummarizedExperiment(
            assays=new_assays,
            row_data=self.row_data.copy(),
            col_data=new_col_data,
            metadata=self.metadata.copy(),
        )

    # ─── Reconstruction ─────────────────────────────────────────────

    def to_dataframe(self, sample: str | None = None) -> pd.DataFrame:
        """
        Reconstruct a flat DataFrame combining assays and row_data.

        For a single-sample SE or when `sample` is specified, produces a
        DataFrame with columns: tax_id, abundance, estimated counts, plus
        all row_data columns (taxonomy). This matches the format expected
        by downstream functions like compute_diversity() and get_sankey_data().

        Parameters
        ----------
        sample : str, optional
            Which sample column to extract. Required if n_samples > 1.
            If the SE has exactly 1 sample, it is used automatically.

        Returns
        -------
        pd.DataFrame
            Flat DataFrame with tax_id as a column (not index).
        """
        if sample is None:
            if self.n_samples == 1:
                sample = self.sample_ids[0]
            else:
                raise ValueError(
                    f"SE has {self.n_samples} samples; specify which one with `sample=`."
                )

        if sample not in self.sample_ids:
            raise ValueError(f"Sample '{sample}' not found. Available: {list(self.sample_ids)}")

        # Start with assay values for this sample
        data = pd.DataFrame(index=self.feature_ids)
        if "abundance" in self.assays:
            data["abundance"] = self.assays["abundance"][sample]
        if "counts" in self.assays:
            data["estimated counts"] = self.assays["counts"][sample]

        # Add row_data (taxonomy)
        if not self.row_data.empty:
            data = data.join(self.row_data)

        # Move index (tax_id) to column
        data = data.reset_index().rename(columns={"index": "tax_id"})

        return data

    # ─── Representation ──────────────────────────────────────────────

    def __repr__(self) -> str:
        assay_info = ", ".join(f"'{k}'" for k in self.assay_names)
        return (
            f"SummarizedExperiment("
            f"n_features={self.n_features}, "
            f"n_samples={self.n_samples}, "
            f"assays=[{assay_info}])"
        )
