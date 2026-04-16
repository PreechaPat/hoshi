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

    # Force abundance column to float
    df[val_col] = pd.to_numeric(df[val_col], errors="coerce").fillna(0.0)
    
    link_aggregations = {}
    node_dict = {}

    # Iterate row by row to build unbroken chains, skipping empty taxonomies
    for _, row in df.iterrows():
        val = row[val_col]
        if val <= 0:
            continue
            
        # Extract the sequence of valid, non-empty labels for this specific bacteria
        valid_nodes = []
        for i, level in enumerate(levels):
            label = str(row[level]).strip()
            if label and label.lower() != "nan":
                valid_nodes.append((i, label))
                
        # Create links that bridge directly between the valid levels
        for j in range(len(valid_nodes) - 1):
            src_i, src_label = valid_nodes[j]
            tgt_i, tgt_label = valid_nodes[j+1]
            
            src_id = f"{src_i}::{src_label}"
            tgt_id = f"{tgt_i}::{tgt_label}"
            
            # Store the clean label and column index for positioning
            node_dict[src_id] = {"label": src_label, "level": src_i}
            node_dict[tgt_id] = {"label": tgt_label, "level": tgt_i}
            
            # Accumulate the abundance values for identical links
            pair = (src_id, tgt_id)
            link_aggregations[pair] = link_aggregations.get(pair, 0.0) + val

    # Assign integer indices based on the unique IDs
    sorted_node_ids = sorted(list(node_dict.keys()))
    id_to_idx = {node_id: idx for idx, node_id in enumerate(sorted_node_ids)}
    
    nodes = []
    for node_id in sorted_node_ids:
        label = node_dict[node_id]["label"]
        level_idx = node_dict[node_id]["level"]
        nodes.append({
            "label": label,
            "color": hex_color_from_label(label), 
            "x_pos": level_idx / (len(levels) - 1)
        })

    formatted_links = [
        {
            "source": id_to_idx[src_id],
            "target": id_to_idx[tgt_id],
            "value": value,
            "color": hex_color_from_label(node_dict[src_id]["label"])
        }
        for (src_id, tgt_id), value in link_aggregations.items()
    ]

    return {
        "nodes": nodes,
        "links": formatted_links,
        "levels": levels
    }

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
