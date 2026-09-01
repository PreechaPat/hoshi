"""
Egress module — export SummarizedExperiment to various output formats.

Currently supported:
  - Kraken2 report format
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hoshi.lib.experiment import SummarizedExperiment

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
) -> dict[tuple[str, str, str], float]:
    """
    Build a mapping of (rank, name, tax_id) → direct counts.

    For species-level entries, direct counts are the actual estimated counts.
    For higher-level taxa, direct counts are 0 (they only have clade counts).

    Returns a dict keyed by (rank, name, tax_id) with direct counts as values.
    Also returns a separate dict for clade counts.
    """
    # Collect species-level entries with their lineage
    entries: list[dict] = []

    for tax_id in counts.index:
        count = counts[tax_id]
        if count <= 0:
            continue

        row = row_data.loc[tax_id]
        lineage = {}
        for rank in available_ranks:
            val = row.get(rank)
            if pd.notna(val) and str(val).strip():
                lineage[rank] = str(val).strip()

        entries.append({
            "tax_id": str(tax_id),
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
