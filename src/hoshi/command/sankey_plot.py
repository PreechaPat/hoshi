import pandas as pd
import plotly.graph_objects as go
import hashlib
import argparse


def hex_color_from_label(label):
    h = hashlib.md5(label.encode()).hexdigest()[:6]
    return f"#{h}"


def sankey_from_dataframe(df, level_cols, value_col="value", title=None):
    links = []
    node_levels = {}
    for i in range(len(level_cols) - 1):
        src_col = level_cols[i]
        tgt_col = level_cols[i + 1]
        step_df = df[[src_col, tgt_col, value_col]].copy()
        step_df.columns = ["source", "target", "value"]
        step_df = step_df[step_df["target"] != ""]
        links.append(step_df)

        for label in df[src_col].unique():
            node_levels[label] = i
        for label in df[tgt_col].unique():
            node_levels[label] = i + 1

    all_links = pd.concat(links, ignore_index=True)
    grouped = all_links.groupby(["source", "target"], as_index=False)["value"].sum()

    labels = pd.unique(grouped[["source", "target"]].values.ravel()).tolist()
    label_to_index = {label: i for i, label in enumerate(labels)}
    node_colors = [hex_color_from_label(label) for label in labels]

    sources = grouped["source"].map(label_to_index).tolist()
    targets = grouped["target"].map(label_to_index).tolist()
    values = grouped["value"].tolist()
    link_colors = [node_colors[src] for src in sources]

    # Assign only x positions
    node_x = [node_levels.get(label, 0) / (len(level_cols) - 1) for label in labels]

    fig = go.Figure(
        go.Sankey(
            arrangement="snap",
            node=dict(
                pad=50,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=labels,
                color=node_colors,
                x=node_x,  # Only x set, y is auto-handled
            ),
            link=dict(source=sources, target=targets, value=values, color=link_colors),
        )
    )

    # Axis-level annotations
    for i, level in enumerate(level_cols):
        fig.add_annotation(
            x=i / (len(level_cols) - 1),
            y=1.05,
            text=f"<b>{level}</b>",
            showarrow=False,
            font=dict(size=16),
            xanchor="center",
        )

    fig.update_layout(
        width=1200,
        height=800,
        font=dict(size=16),
        margin=dict(t=100, l=50, r=50, b=50),
        title=dict(
            text=f"<b>{title or 'Sankey Diagram'}</b>",
            font=dict(size=22),
            x=0.5,
            y=0.95,
            xanchor="center",
            yanchor="top",
        ),
    )

    return fig


# ------------------------
# Script Entry Point
# ------------------------


def main():
    parser = argparse.ArgumentParser(description="Plot a Sankey diagram from hierarchical TSV data.")
    parser.add_argument(
        "input_file",
        help="Input TSV file with hierarchical columns and abundance column.",
    )
    parser.add_argument(
        "--value_col",
        default="abundance",
        help="Column representing flow values (default: 'abundance').",
    )
    parser.add_argument("-o", "--output_file", help="Path to save the Sankey plot HTML.")
    parser.add_argument("--title", help="Plot title.")
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["species", "genus", "family", "order", "class", "phylum", "kingdom"],
        help="Ordered taxonomic ranks to use as levels.",
    )

    args = parser.parse_args()

    # Load TSV data
    df = pd.read_csv(args.input_file, sep="\t").fillna("")

    # Verify required columns
    required_cols = args.levels + [args.value_col]
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in input: {missing}")

    # Create plot
    fig = sankey_from_dataframe(df, level_cols=args.levels, value_col=args.value_col, title=args.title)

    if args.output_file:
        fig.write_html(args.output_file)
        print(f"Sankey plot saved to {args.output_file}")
    else:
        fig.show()


if __name__ == "__main__":
    main()
