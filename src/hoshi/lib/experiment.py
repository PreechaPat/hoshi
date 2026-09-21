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

# Canonical taxonomy ranks, specific → broad. Taxonomy (and tax_id) are nullable
# per-OTU annotations, not the feature identity: an OTU/ASV is one sequence that
# *may* be assigned a lineage (possibly only down to genus, or not at all).
TAXONOMY_RANKS: tuple[str, ...] = (
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
)

# Column that carries the per-feature id when a flat table is produced from an
# experiment (see :meth:`SummarizedExperiment.to_dataframe`). Features are
# OTUs/ASVs — one sequence each — identified by the classifier's own id, scoped
# per sample to stay globally unique across a multi-sample experiment.
FEATURE_ID_COLUMN = "feature_id"

# tax_id values classifiers use for control/meta rows that are not organisms
# (EMU emits these). Dropped before any species-level view is rendered.
META_TAX_IDS: frozenset[str] = frozenset(
    {"unmapped", "mapped_filtered", "mapped_unclassified"}
)


def make_feature_ids(n: int, *, prefix: str = "F", start: int = 1) -> list[str]:
    """Return ``n`` running feature identifiers (``F00001``, ``F00002`` …).

    Only used as a fallback when a classifier provides no id of its own. Ids are
    zero-padded to at least five digits and widen automatically beyond 99999, so
    string ordering matches numeric ordering. ``prefix`` lets callers scope ids
    per sample.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    width = max(5, len(str(start + n - 1))) if n else 5
    return [f"{prefix}{i:0{width}d}" for i in range(start, start + n)]


def aggregate_to_species(
    df: pd.DataFrame,
    *,
    value_columns: tuple[str, ...] = ("abundance", "estimated counts"),
    sequence_identity_column: str = "sequence_identity",
    keep_blank_tax_id_separate: bool = True,
) -> pd.DataFrame:
    """Aggregate a per-OTU flat table up to one row per species (``tax_id``).

    OTUs that share a ``tax_id`` are summed across ``value_columns`` and keep
    the first taxonomy lineage seen. When ``keep_blank_tax_id_separate`` is True
    (option iii), OTUs with a missing/blank ``tax_id`` are **not** collapsed
    together — each stays as its own row with an empty ``tax_id`` — so read
    counts still total correctly and unassigned sequences are never merged into
    one bogus taxon.

    Reports render at species/genus level, so when a ``sequence_identity_column``
    is present it is aggregated with ``max()`` (the best-matching OTU wins for the
    species).

    TODO: max() is a placeholder for sequence-identity aggregation; revisit later.

    The input is expected to have a ``tax_id`` column plus any of the
    :data:`TAXONOMY_RANKS`. Non-value, non-taxonomy columns are dropped.
    """
    if "tax_id" not in df.columns:
        raise ValueError("aggregate_to_species requires a 'tax_id' column.")

    work = df.copy()
    tax_id = work["tax_id"].astype("string")
    blank = tax_id.isna() | (tax_id.str.strip() == "")

    present_values = [c for c in value_columns if c in work.columns]
    for col in present_values:
        work[col] = pd.to_numeric(work[col], errors="coerce").fillna(0.0)

    has_identity = sequence_identity_column in work.columns
    tax_cols = [c for c in TAXONOMY_RANKS if c in work.columns]

    agg_map: dict[str, str] = {c: "sum" for c in present_values}
    agg_map.update({c: "first" for c in tax_cols})
    if has_identity:
        agg_map[sequence_identity_column] = "max"

    out_cols = ["tax_id", *present_values, *tax_cols]
    if has_identity:
        out_cols.append(sequence_identity_column)

    # Rows with a real tax_id collapse by tax_id.
    assigned = work[~blank].copy()
    assigned["tax_id"] = tax_id[~blank]
    collapsed = (
        assigned.groupby("tax_id", as_index=False, dropna=False).agg(agg_map)
        if not assigned.empty
        else assigned.reindex(columns=out_cols)
    )

    # Unassigned rows are kept individually (option iii) with a blank tax_id.
    unassigned = work[blank].copy()
    if not unassigned.empty:
        unassigned["tax_id"] = ""
        unassigned = unassigned[out_cols]

    out = pd.concat(
        [f for f in (collapsed, unassigned) if not f.empty],
        ignore_index=True,
    ) if (not collapsed.empty or not unassigned.empty) else collapsed

    sort_col = present_values[0] if present_values else None
    if sort_col:
        out = out.sort_values(sort_col, ascending=False)
    return out.reset_index(drop=True)


def species_view(
    df: pd.DataFrame,
    *,
    top: int | None = None,
    drop_meta: bool = True,
) -> pd.DataFrame:
    """Collapse a per-OTU flat table into the species view reports render.

    This is the shared "OTU rows → species rows" rollup used by every report and
    the count-table export. It is the single place that:

    1. drops classifier control/meta rows (:data:`META_TAX_IDS`) so they never
       surface as organisms,
    2. aggregates OTUs up to one row per ``tax_id`` (:func:`aggregate_to_species`;
       blank tax_ids stay individual so counts still total correctly), and
    3. returns rows sorted by abundance descending, optionally trimmed to ``top``.

    ``df`` is a per-OTU flat frame as produced by
    :meth:`SummarizedExperiment.to_dataframe` (``tax_id`` + taxonomy + value
    columns, optionally a ``sequence_identity`` column which is aggregated with
    ``max()``). Callers keep the un-aggregated ``df`` for OTU-level products
    (diversity, Sankey); this only owns the species collapse.
    """
    work = df
    if drop_meta and "tax_id" in work.columns:
        work = work[~work["tax_id"].astype(str).isin(META_TAX_IDS)].copy()

    rolled = aggregate_to_species(work)
    if top is not None:
        rolled = rolled.head(top).reset_index(drop=True)
    return rolled


def _scope_feature_id(sample: str, feature_id: Any) -> str:
    """Prefix a raw feature id with its sample name to keep it globally unique.

    Feature ids from a classifier (e.g. Savont's ``final_consensus_0`` or EMU's
    ``tax_id``) are only unique *within* one sample. When several samples are
    combined into one experiment the same raw id can appear in more than one
    sample, so we scope it as ``<sample>:<feature_id>`` at merge time. Producers
    keep the raw id; scoping is applied only here (see :func:`combine_experiments`).
    """
    return f"{sample}:{feature_id}"


def combine_experiments(
    experiments: list[SummarizedExperiment],
) -> SummarizedExperiment:
    """Combine per-sample experiments into one experiment.

    Each input is expected to be a single-sample :class:`SummarizedExperiment`
    whose feature index carries the classifier's **raw** feature ids (unscoped).

    Feature ids are only unique *within* a sample, so once the combined result
    holds more than one sample the same raw id could collide across samples.
    This is the one place that resolves that: when the result spans multiple
    samples, feature ids are scoped per sample (``<sample>:<feature_id>``); a
    single-sample result keeps its raw ids untouched. Producers therefore never
    scope, and a lone sample never carries a pointless sample prefix.

    The result is the block-sparse union used across the codebase:

    - ``assays`` — every named assay is aligned into a (features × samples)
      matrix; a feature absent from a sample is filled with ``0``. All inputs
      must expose the same assay names.
    - ``row_data`` — per-feature annotations, unioned across samples (first
      value wins for a given id; ids are unique per sample, and scoped when
      multi-sample, so there is no real conflict).
    - ``col_data`` — sample metadata rows concatenated in input order.
    - ``metadata`` — ``source`` is carried through when all inputs agree;
      ``sequence_identity`` is merged into ``{sample: {feature_id: pct}}`` with keys
      matching the (raw or scoped) feature index.
    """
    if not experiments:
        return SummarizedExperiment()

    assay_names = list(experiments[0].assay_names)
    for exp in experiments:
        if list(exp.assay_names) != assay_names:
            raise ValueError(
                "combine_experiments requires all inputs to share the same "
                f"assay names; got {list(exp.assay_names)} vs {assay_names}."
            )

    # Scope ids only when the combined result actually spans multiple samples;
    # a single sample keeps its classifier-native ids.
    total_samples = sum(exp.n_samples for exp in experiments)
    scope = total_samples > 1

    def key(sample: str, feature_id: Any) -> str:
        return _scope_feature_id(sample, feature_id) if scope else str(feature_id)

    per_assay_series: dict[str, dict[str, pd.Series]] = {n: {} for n in assay_names}
    row_data_frames: list[pd.DataFrame] = []
    col_data_frames: list[pd.DataFrame] = []
    sequence_identity: dict[str, dict[str, float]] = {}
    sources: set[str] = set()

    for exp in experiments:
        for sample in exp.sample_ids:
            sample = str(sample)
            index = exp.feature_ids
            keyed_index = pd.Index(
                [key(sample, fid) for fid in index],
                name=FEATURE_ID_COLUMN,
            )

            for name in assay_names:
                col = exp.assays[name][sample]
                per_assay_series[name][sample] = pd.Series(
                    col.to_numpy(), index=keyed_index
                )

            if not exp.row_data.empty:
                rd = exp.row_data.copy()
                rd.index = keyed_index
                row_data_frames.append(rd)

            # Carry per-feature sequence identity, keyed to match the feature index
            # (raw feature id, or sample-scoped when multi-sample). Tolerate an
            # already-scoped key so re-combining is idempotent.
            per_sample_identity = exp.metadata.get("sequence_identity", {})
            raw_identity = per_sample_identity.get(sample, {})
            if raw_identity:
                sequence_identity[sample] = {
                    key(sample, _strip_sample_prefix(sample, fid)): pct
                    for fid, pct in raw_identity.items()
                }

        if not exp.col_data.empty:
            col_data_frames.append(exp.col_data)
        src = exp.metadata.get("source")
        if src is not None:
            sources.add(src)

    assays: dict[str, pd.DataFrame] = {}
    reference_index: pd.Index | None = None
    for name in assay_names:
        matrix = pd.DataFrame(per_assay_series[name]).fillna(0.0)
        matrix.index.name = FEATURE_ID_COLUMN
        assays[name] = matrix
        if reference_index is None:
            reference_index = matrix.index

    if row_data_frames:
        row_data = pd.concat(row_data_frames)
        row_data = row_data[~row_data.index.duplicated(keep="first")]
        row_data = row_data.reindex(reference_index)
        row_data.index.name = FEATURE_ID_COLUMN
    else:
        row_data = pd.DataFrame()

    col_data = (
        pd.concat(col_data_frames) if col_data_frames else pd.DataFrame()
    )

    metadata: dict[str, Any] = {}
    if len(sources) == 1:
        metadata["source"] = next(iter(sources))
    if sequence_identity:
        metadata["sequence_identity"] = sequence_identity

    return SummarizedExperiment(
        assays=assays,
        row_data=row_data,
        col_data=col_data,
        metadata=metadata,
    )


def _strip_sample_prefix(sample: str, feature_id: Any) -> str:
    """Return ``feature_id`` without a leading ``<sample>:`` prefix, if present.

    Producers emit raw (unscoped) sequence-identity keys, but tolerate an already
    ``<sample>:``-scoped key so re-combining an experiment is idempotent.
    """
    fid = str(feature_id)
    prefix = f"{sample}:"
    return fid[len(prefix):] if fid.startswith(prefix) else fid


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

    def to_dataframe(
        self,
        sample: str | None = None,
        *,
        sequence_identity: dict[str, float] | None = None,
    ) -> pd.DataFrame:
        """
        Reconstruct a flat DataFrame combining assays and row_data.

        For a single-sample SE or when `sample` is specified, produces a
        DataFrame with columns: feature_id, abundance, estimated counts, plus
        all row_data columns (``tax_id`` and taxonomy). The feature index is the
        per-feature (OTU/ASV) id; ``tax_id`` is a regular (nullable) column
        carried in row_data. This matches the format expected by downstream
        functions like compute_diversity() and get_sankey_data(), which key on
        the ``tax_id`` column rather than the index.

        This is an un-pivot, not an aggregation: one OTU row in, one row out.
        Species collapsing is a separate, explicit step (:func:`species_view`).

        Parameters
        ----------
        sample : str, optional
            Which sample column to extract. Required if n_samples > 1.
            If the SE has exactly 1 sample, it is used automatically.
        sequence_identity : dict[str, float], optional
            Per-feature estimated sequence identity keyed by ``feature_id``. When
            given, it is attached as a ``sequence_identity`` column (features with
            no entry get ``NaN``), so callers no longer re-join it by hand.

        Returns
        -------
        pd.DataFrame
            Flat DataFrame with ``feature_id`` as a column (the former index)
            and ``tax_id`` as a regular column when present in row_data.
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

        # Add row_data (tax_id + taxonomy)
        if not self.row_data.empty:
            data = data.join(self.row_data)

        # Surface the feature index (OTU/ASV id) as a regular ``feature_id``
        # column. The feature index is always named ``feature_id`` (single-sample
        # producers emit raw ids; combine_experiments scopes them but keeps the
        # name), so reset_index() yields that column directly. Fall back to
        # renaming a bare ``index`` column defensively.
        had_named_index = data.index.name is not None
        data = data.reset_index()
        if not had_named_index and "index" in data.columns:
            data = data.rename(columns={"index": "feature_id"})

        # Attach per-feature sequence identity (keyed by feature_id) when supplied,
        # callers don't re-join it by hand before aggregating to species.
        if sequence_identity and "feature_id" in data.columns:
            data["sequence_identity"] = data["feature_id"].astype(str).map(sequence_identity)

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