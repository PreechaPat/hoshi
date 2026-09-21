"""
Egress module — export SummarizedExperiment to various output formats.

Currently supported:
  - Kraken2 report format
  - Species count table (flat TSV: abundance, estimated counts, optional
    sequence identity, tax_id, taxonomy)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hoshi.lib.experiment import SummarizedExperiment, species_view

# Mapping from Emu taxonomy column names to Kraken2 rank codes
_RANK_MAP: dict[str, str] = {
    "superkingdom": "D",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
}

# Ordered list of taxonomic ranks (broad → specific)
_RANK_ORDER: list[str] = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def experiment_to_kraken2(
    experiment: SummarizedExperiment,
    sample: str | None = None,
) -> str:
    """
    Convert a SummarizedExperiment to Kraken2 report format.

    Produces a tab-delimited string matching Kraken2's sample report format:
      1. Percentage of fragments covered by the clade rooted at this taxon
      2. Number of fragments covered by the clade rooted at this taxon
      3. Number of fragments assigned directly to this taxon
      4. Rank code (D, P, C, O, F, G, S)
      5. NCBI taxonomic ID
      6. Indented scientific name

    Since Emu provides species-level classifications (not a full read-level
    tree walk), this function:
      - Uses estimated counts as fragment counts per species.
      - Aggregates counts up the taxonomy tree to compute clade totals.
      - Reports each taxonomic rank with proper indentation.

    Parameters
    ----------
    experiment : SummarizedExperiment
        Must have an 'counts' assay and row_data with taxonomy columns.

    sample : str, optional
        Which sample to export. Required if n_samples > 1.

    Returns
    -------
    str
        Kraken2 report as a tab-delimited string.

    Raises
    ------
    ValueError
        If required assays or taxonomy columns are missing.
    """
    if "counts" not in experiment.assay_names:
        raise ValueError(
            "Experiment must have a 'counts' assay for Kraken2 conversion."
        )

    if sample is None:
        if experiment.n_samples == 1:
            sample = experiment.sample_ids[0]
        else:
            raise ValueError(
                f"Experiment has {experiment.n_samples} samples; "
                f"specify which one with `sample=`."
            )

    if sample not in experiment.sample_ids:
        raise ValueError(
            f"Sample '{sample}' not found. Available: {list(experiment.sample_ids)}"
        )

    # Get counts for this sample
    counts = experiment.assays["counts"][sample]
    row_data = experiment.row_data

    # Check that we have at least some taxonomy columns
    available_ranks = [r for r in _RANK_ORDER if r in row_data.columns]
    if not available_ranks:
        raise ValueError(
            "row_data must have at least one taxonomy column "
            f"({', '.join(_RANK_ORDER)})."
        )

    total_counts = counts.sum()
    if total_counts == 0:
        return ""

    # Build the taxonomy tree by aggregating counts at each rank level
    # Each node: (rank, name) → direct_counts, clade_counts
    tree = _build_taxonomy_tree(counts, row_data, available_ranks)

    # Format as Kraken2 report
    lines = _format_kraken2_report(tree, total_counts, available_ranks)

    return "\n".join(lines) + "\n"


def _build_taxonomy_tree(
    counts: pd.Series,
    row_data: pd.DataFrame,
    available_ranks: list[str],
) -> list[dict]:
    """
    Collect nonzero-count features as entries with their taxonomy lineage.

    Each entry is a dict with keys ``tax_id``, ``count``, and ``lineage``
    (a rank → name mapping, skipping blank ranks).
    """
    # Collect per-feature entries with their lineage. The feature index is a
    # per-OTU id, so the NCBI tax_id is read from row_data (nullable) rather than
    # the index; features without a tax_id fall back to "0".
    entries: list[dict] = []

    for feature_id in counts.index:
        count = counts[feature_id]
        if count <= 0:
            continue

        row = row_data.loc[feature_id]
        lineage = {}
        for rank in available_ranks:
            val = row.get(rank)
            if pd.notna(val) and str(val).strip():
                lineage[rank] = str(val).strip()

        tax_id_val = row.get("tax_id") if "tax_id" in row_data.columns else None
        tax_id = (
            str(tax_id_val).strip()
            if pd.notna(tax_id_val) and str(tax_id_val).strip()
            else "0"
        )

        entries.append({
            "tax_id": tax_id,
            "count": float(count),
            "lineage": lineage,
        })

    return entries


def _format_kraken2_report(
    entries: list[dict],
    total_counts: float,
    available_ranks: list[str],
) -> list[str]:
    """
    Format entries into Kraken2 report lines.

    Strategy:
    - Walk through each rank from broad to specific.
    - At each rank, group by the taxon name.
    - Clade counts = sum of all species that pass through this taxon.
    - Direct counts = 0 for non-species ranks; actual counts for species.
    """
    # Build clade counts at each rank level
    # clade_nodes[(rank, name)] = total counts of all species under this node
    clade_counts: dict[tuple[str, str], float] = {}
    # Track tax_ids at species level
    species_info: dict[str, dict] = {}  # name → {tax_id, count}

    for entry in entries:
        lineage = entry["lineage"]
        count = entry["count"]

        for rank in available_ranks:
            name = lineage.get(rank)
            if name:
                key = (rank, name)
                clade_counts[key] = clade_counts.get(key, 0.0) + count

        # Species-level direct assignment
        species_name = lineage.get("species")
        if species_name:
            if species_name not in species_info:
                species_info[species_name] = {
                    "tax_id": entry["tax_id"],
                    "count": 0.0,
                }
            species_info[species_name]["count"] += count

    # Now build the hierarchical report
    # We output an unclassified line first (with 0), then walk the tree
    lines: list[str] = []

    # Add unclassified line (0 for Emu since all reads are classified)
    lines.append(_format_line(0.0, 0, 0, "U", "0", "unclassified"))

    # Add root line — all classified reads fall under root (taxid 1)
    classified_count = int(round(total_counts))
    pct = 100.0 if total_counts > 0 else 0.0
    lines.append(_format_line(pct, classified_count, 0, "R", "1", "root"))

    # Walk through ranks in order, outputting taxa sorted by clade counts
    # Group entries by their lineage path for proper ordering
    _output_tree(
        lines, entries, available_ranks, total_counts, species_info
    )

    return lines


def _output_tree(
    lines: list[str],
    entries: list[dict],
    available_ranks: list[str],
    total_counts: float,
    species_info: dict[str, dict],
) -> None:
    """
    Output the taxonomy tree in Kraken2 hierarchical order.

    Uses a depth-first traversal: for each rank, group by ancestors,
    then output children sorted by clade count (descending).
    """
    # Build tree structure from entries
    # Strategy: for each rank level, collect unique taxa and their clade counts
    rank_nodes: dict[str, dict[str, dict]] = {}  # rank → name → node

    for rank in available_ranks:
        rank_nodes[rank] = {}

    # First pass: collect clade counts and lineage info
    for entry in entries:
        lineage = entry["lineage"]
        count = entry["count"]

        for rank in available_ranks:
            name = lineage.get(rank)
            if not name:
                continue
            if name not in rank_nodes[rank]:
                rank_nodes[rank][name] = {
                    "name": name,
                    "rank": rank,
                    "tax_id": entry["tax_id"] if rank == "species" else "0",
                    "direct_count": 0.0,
                    "clade_count": 0.0,
                }
            rank_nodes[rank][name]["clade_count"] += count

    # Set direct counts for species
    for name, info in species_info.items():
        if name in rank_nodes.get("species", {}):
            rank_nodes["species"][name]["direct_count"] = info["count"]
            rank_nodes["species"][name]["tax_id"] = info["tax_id"]

    # Build parent-child relationships based on lineage paths
    # For each entry, trace the lineage and establish links
    children_map: dict[tuple[str, str], set[tuple[str, str]]] = {}
    # (parent_rank, parent_name) → set of (child_rank, child_name)

    for entry in entries:
        lineage = entry["lineage"]
        # Find consecutive rank pairs in the lineage
        lineage_pairs = []
        for rank in available_ranks:
            name = lineage.get(rank)
            if name:
                lineage_pairs.append((rank, name))

        for i in range(len(lineage_pairs) - 1):
            parent = lineage_pairs[i]
            child = lineage_pairs[i + 1]
            if parent not in children_map:
                children_map[parent] = set()
            children_map[parent].add(child)

    # Find root nodes: the first rank that actually has populated entries.
    # (superkingdom is often blank in Emu output)
    root_rank = None
    for rank in available_ranks:
        if rank_nodes[rank]:
            root_rank = rank
            break

    if root_rank is None:
        return  # No taxonomy data at all

    # If superkingdom is not populated, inject an "unspecified" domain line
    has_domain = "superkingdom" in available_ranks and rank_nodes.get("superkingdom")
    starting_depth = 1  # depth under root

    if not has_domain:
        classified_count = int(round(total_counts))
        pct = 100.0 if total_counts > 0 else 0.0
        indented_name = "  " * starting_depth + "unspecified"
        lines.append(
            _format_line(pct, classified_count, 0, "D", "0", indented_name)
        )
        starting_depth += 1

    root_taxa = sorted(
        rank_nodes[root_rank].values(),
        key=lambda x: x["clade_count"],
        reverse=True,
    )

    # Depth-first output
    for node in root_taxa:
        _output_node(
            lines, node, rank_nodes, children_map,
            available_ranks, total_counts, depth=starting_depth,
        )


def _output_node(
    lines: list[str],
    node: dict,
    rank_nodes: dict[str, dict[str, dict]],
    children_map: dict[tuple[str, str], set[tuple[str, str]]],
    available_ranks: list[str],
    total_counts: float,
    depth: int,
) -> None:
    """Recursively output a node and its children."""
    rank = node["rank"]
    name = node["name"]
    rank_code = _RANK_MAP.get(rank, "U")

    pct = (node["clade_count"] / total_counts) * 100.0 if total_counts > 0 else 0.0
    clade = int(round(node["clade_count"]))
    direct = int(round(node["direct_count"]))
    tax_id = node["tax_id"]

    # Indentation: 2 spaces per depth level
    indented_name = "  " * depth + name

    lines.append(_format_line(pct, clade, direct, rank_code, tax_id, indented_name))

    # Output children
    key = (rank, name)
    if key in children_map:
        child_keys = children_map[key]
        # Sort children by clade count descending
        child_nodes = []
        for child_rank, child_name in child_keys:
            if child_name in rank_nodes.get(child_rank, {}):
                child_nodes.append(rank_nodes[child_rank][child_name])

        child_nodes.sort(key=lambda x: x["clade_count"], reverse=True)

        for child in child_nodes:
            _output_node(
                lines, child, rank_nodes, children_map,
                available_ranks, total_counts, depth=depth + 1,
            )


def _format_line(
    pct: float,
    clade_count: int,
    direct_count: int,
    rank_code: str,
    tax_id: str,
    name: str,
) -> str:
    """Format a single Kraken2 report line."""
    return f"{pct:6.2f}\t{clade_count}\t{direct_count}\t{rank_code}\t{tax_id}\t{name}"


def write_kraken2_report(
    experiment: SummarizedExperiment,
    output: str | Path,
    sample: str | None = None,
) -> None:
    """
    Write a Kraken2-format report to a file.

    Parameters
    ----------
    experiment : SummarizedExperiment
        The experiment to convert.
    output : str or Path
        Output file path.
    sample : str, optional
        Which sample to export. Required if n_samples > 1.
    """
    report = experiment_to_kraken2(experiment, sample=sample)
    Path(output).write_text(report)


# ── Species count table ──────────────────────────────────────────────
#
# Savont's native species output (``species_abundance.tsv``) lists relative
# abundance next to the taxonomy lineage but drops both the NCBI ``tax_id`` and
# any read count. This egress mirrors that species table while adding the two
# columns Savont omits: ``tax_id`` and ``estimated counts``. Taxonomy stays in
# separate rank columns (one row per species) rather than being collapsed into a
# single lineage string.

# Column order of the emitted species count table. Leads with the values of
# interest (relative abundance, estimated count, estimated sequence identity)
# and the identifier (tax_id), followed by the taxonomy lineage specific → broad.
# ``sequence_identity`` is only emitted when the classifier supplies it (Savont
# does; EMU does not) — see :func:`build_count_table`.
_COUNT_TABLE_COLUMNS: list[str] = [
    "relative_abundance",
    "estimated_count",
    "sequence_identity",
    "tax_id",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
]


def build_count_table(
    experiment: SummarizedExperiment,
    sample: str | None = None,
) -> pd.DataFrame:
    """Build a species-level count table from a per-OTU experiment.

    OTUs are rolled up to one row per species via the shared
    :func:`hoshi.lib.experiment.species_view` (drops classifier meta rows,
    aggregates by ``tax_id``, sorts by abundance; OTUs with a blank ``tax_id``
    stay individual so counts still total correctly).

    Mirrors Savont's ``species_abundance.tsv`` (relative abundance + taxonomy
    lineage) but adds the columns Savont drops: the NCBI ``tax_id``, the
    ``estimated_count``, and the per-species estimated ``sequence_identity``.
    ``sequence_identity`` is always a column for a stable schema; classifiers that
    do not report it (EMU) leave it N/A. Taxonomy is kept as separate rank columns;
    there is no OTU column.

    Parameters
    ----------
    experiment : SummarizedExperiment
        Per-OTU experiment. Must have a ``counts`` assay; ``abundance`` is used
        when present and otherwise derived from counts. ``row_data`` supplies
        ``tax_id`` + taxonomy columns; ``metadata["sequence_identity"]`` supplies the
        optional per-OTU sequence identity (Savont).
    sample : str, optional
        Which sample to export. Required if the experiment has > 1 sample.

    Returns
    -------
    pandas.DataFrame
        Columns: relative_abundance, estimated_count, sequence_identity, tax_id,
        species, genus, family, order, class, phylum, superkingdom.
        ``sequence_identity`` is N/A for classifiers that do not report it (EMU).
        Sorted by relative_abundance (descending) with a positional integer
        index.

    Raises
    ------
    ValueError
        If the ``counts`` assay is missing, no taxonomy column is present, or
        the requested sample is absent / ambiguous.
    """
    if "counts" not in experiment.assay_names:
        raise ValueError(
            "Experiment must have a 'counts' assay for count-table conversion."
        )
    if not any(r in experiment.row_data.columns for r in _RANK_ORDER):
        raise ValueError(
            "row_data must have at least one taxonomy column "
            f"({', '.join(_RANK_ORDER)})."
        )
    if "tax_id" not in experiment.row_data.columns:
        raise ValueError("row_data must have a 'tax_id' column.")

    # Per-OTU flat frame (sample resolution/validation lives in to_dataframe),
    # then roll up to species via the shared view so meta/control rows are
    # dropped consistently with the reports. Per-OTU sequence identity (Savont) is
    # attached here and aggregated to the species with max() by species_view.
    per_sample_identity = experiment.metadata.get("sequence_identity") or {}
    resolved = sample if sample is not None else (
        experiment.sample_ids[0] if experiment.n_samples == 1 else None
    )
    sequence_identity = per_sample_identity.get(str(resolved)) if resolved is not None else None

    per_otu = experiment.to_dataframe(sample=sample, sequence_identity=sequence_identity or None)
    if "abundance" not in per_otu.columns:
        counts = per_otu["estimated counts"].astype(float)
        total = counts.sum()
        per_otu["abundance"] = counts / total if total > 0 else counts * 0.0

    species = species_view(per_otu).rename(
        columns={"abundance": "relative_abundance", "estimated counts": "estimated_count"}
    )

    # Whole-read counts (Savont) render as ints; keep fractional EM counts (Emu).
    counts_series = species["estimated_count"]
    if counts_series.dropna().mod(1).eq(0).all():
        species["estimated_count"] = counts_series.round().astype("Int64")

    # ``sequence_identity`` is always a column for a stable schema. Savont
    # supplies it (aggregated to the species with max() by species_view); EMU
    # does not, so those rows stay N/A.
    if "sequence_identity" in species.columns:
        species["sequence_identity"] = pd.to_numeric(
            species["sequence_identity"], errors="coerce"
        ).round(1)

    # Guarantee a stable column set/order even when some ranks/identity absent.
    for col in _COUNT_TABLE_COLUMNS:
        if col not in species.columns:
            species[col] = pd.NA

    return (
        species[_COUNT_TABLE_COLUMNS]
        .sort_values("relative_abundance", ascending=False)
        .reset_index(drop=True)
    )


def experiment_to_count_table(
    experiment: SummarizedExperiment,
    sample: str | None = None,
) -> str:
    """Render :func:`build_count_table` as a tab-delimited string."""
    table = build_count_table(experiment, sample=sample)
    return table.to_csv(sep="\t", index=False)


def write_count_table(
    experiment: SummarizedExperiment,
    output: str | Path,
    sample: str | None = None,
) -> None:
    """Write the species count table to ``output`` as a TSV file."""
    Path(output).write_text(experiment_to_count_table(experiment, sample=sample))