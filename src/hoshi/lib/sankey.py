from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.graph_objects as go

# Constants
TAXONOMY_LEVELS = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")
VALUE_COLUMNS = ("estimated counts", "abundance")


def hex_color_from_label(label: str) -> str:
    """Generate a consistent hex color for a label using MD5."""
    return f"#{hashlib.md5(label.encode()).hexdigest()[:6]}"


def get_sankey_data(tsv_path: str | Path) -> dict[str, Any]:
    """
    Converts TSV data into a dictionary structure (JSON-compatible)
    representing nodes and links for a Sankey diagram.
    """
    df = pd.read_csv(tsv_path, sep="\t")

    levels = [l for l in TAXONOMY_LEVELS if l in df.columns]
    val_col = next((c for c in VALUE_COLUMNS if c in df.columns), None)

    if len(levels) < 2 or not val_col:
        raise ValueError("Insufficient taxonomy levels or abundance data found in TSV.")

    # FIX 1: Strictly force the value column to numeric, matching your original code
    df[val_col] = pd.to_numeric(df[val_col], errors="coerce").fillna(0.0)

    # FIX 2: Convert missing values to actual empty strings for safe filtering
    for level in levels:
        df[level] = df[level].fillna("").astype(str).str.strip()

    links_list = []
    node_labels = set()
    node_level_map = {}

    for i in range(len(levels) - 1):
        src, tgt = levels[i], levels[i + 1]

        # Filter out rows where source or target is an empty string
        subset = df[(df[src] != "") & (df[tgt] != "")]
        grouped = subset.groupby([src, tgt], as_index=False)[val_col].sum()

        for _, row in grouped.iterrows():
            if row[val_col] <= 0:
                continue

            links_list.append(
                {
                    "source": row[src],
                    "target": row[tgt],
                    "value": row[val_col],  # Now strictly numeric!
                }
            )
            node_labels.update([row[src], row[tgt]])
            node_level_map[row[src]] = i
            node_level_map[row[tgt]] = i + 1

    sorted_labels = sorted(list(node_labels))
    label_to_idx = {label: i for i, label in enumerate(sorted_labels)}

    nodes = []
    for label in sorted_labels:
        nodes.append(
            {
                "label": label,
                "color": hex_color_from_label(label),
                "x_pos": node_level_map.get(label, 0) / (len(levels) - 1),
            }
        )

    formatted_links = [
        {
            "source": label_to_idx[link["source"]],
            "target": label_to_idx[link["target"]],
            "value": link["value"],
            "color": hex_color_from_label(link["source"]),
        }
        for link in links_list
    ]

    return {"nodes": nodes, "links": formatted_links, "levels": levels}


def render_sankey_figure(data: dict[str, Any], title: str = "Sankey Diagram") -> go.Figure:
    """Creates a Plotly Figure object from processed Sankey data."""

    # Calculate how many nodes are in the busiest column to set a safe height
    x_positions = [n["x_pos"] for n in data["nodes"]]
    max_nodes_in_col = max([x_positions.count(x) for x in set(x_positions)]) if x_positions else 1

    # 40 pixels per node guarantees enough room for the node thickness + padding
    dynamic_height = max(600, max_nodes_in_col * 40)

    fig = go.Figure(
        go.Sankey(
            arrangement="snap",
            node=dict(
                pad=15,  # <-- FIX: Reduced from 50 to prevent squashing nodes
                thickness=20,
                label=[n["label"] for n in data["nodes"]],
                color=[n["color"] for n in data["nodes"]],
                x=[n["x_pos"] for n in data["nodes"]],
                line=dict(color="black", width=0.5),
            ),
            link=dict(
                source=[l["source"] for l in data["links"]],
                target=[l["target"] for l in data["links"]],
                value=[float(l["value"]) for l in data["links"]],  # Ensuring float values
                color=[l["color"] for l in data["links"]],
            ),
        )
    )

    # Add Level Annotations
    for i, level in enumerate(data["levels"]):
        fig.add_annotation(
            x=i / (len(data["levels"]) - 1),
            y=1.05,
            text=f"<b>{level.capitalize()}</b>",
            showarrow=False,
            xanchor="center",
            font=dict(size=14),
        )

    fig.update_layout(
        title_text=title,
        font_size=12,
        height=dynamic_height,  # <-- FIX: Apply the dynamic height
        margin=dict(t=80, l=50, r=50, b=50),
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description="Convert EMU TSV to Sankey HTML")
    parser.add_argument("input", help="Path to input TSV file")
    parser.add_argument("-o", "--output", help="Path to output HTML (default: input.html)")
    parser.add_argument("-t", "--title", help="Plot title")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path.with_suffix(".html")

    print(f"Processing {input_path}...")

    try:
        data = get_sankey_data(input_path)
        fig = render_sankey_figure(data, title=args.title or f"Taxonomy Flow: {input_path.name}")
        fig.write_html(output_path)
        print(f"Success! Saved to {output_path}")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
