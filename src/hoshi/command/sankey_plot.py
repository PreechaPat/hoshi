from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Sequence

import pandas as pd
import plotly.graph_objects as go

from hoshi.lib.sankey import SankeyInput, load_sankey_input


def hex_color_from_label(label: str) -> str:
    h = hashlib.md5(label.encode()).hexdigest()[:6]
    return f"#{h}"


def sankey_from_dataframe(df: pd.DataFrame, level_cols: Sequence[str], value_col: str, title: str | None = None) -> go.Figure:
    links: list[pd.DataFrame] = []
    node_levels: dict[str, int] = {}
    for index in range(len(level_cols) - 1):
        src_col = level_cols[index]
        tgt_col = level_cols[index + 1]
        step_df = df[[src_col, tgt_col, value_col]].copy()
        step_df.columns = ["source", "target", "value"]
        step_df = step_df[step_df["target"] != ""]
        links.append(step_df)

        for label in df[src_col].unique():
            if label != "":
                node_levels[label] = index
        for label in df[tgt_col].unique():
            if label != "":
                node_levels[label] = index + 1

    if not links:
        raise ValueError("At least two taxonomy levels with data are required to build a Sankey plot.")

    all_links = pd.concat(links, ignore_index=True)
    grouped = all_links.groupby(["source", "target"], as_index=False)["value"].sum()

    labels = pd.unique(grouped[["source", "target"]].values.ravel()).tolist()
    label_to_index = {label: i for i, label in enumerate(labels)}
    node_colors = [hex_color_from_label(label) for label in labels]

    sources = grouped["source"].map(label_to_index).tolist()
    targets = grouped["target"].map(label_to_index).tolist()
    values = grouped["value"].tolist()
    link_colors = [node_colors[src] for src in sources]

    node_x = [node_levels.get(label, 0) / max(len(level_cols) - 1, 1) for label in labels]

    fig = go.Figure(
        go.Sankey(
            arrangement="snap",
            node=dict(
                pad=50,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=labels,
                color=node_colors,
                x=node_x,
            ),
            link=dict(source=sources, target=targets, value=values, color=link_colors),
        )
    )

    for index, level in enumerate(level_cols):
        fig.add_annotation(
            x=index / max(len(level_cols) - 1, 1),
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


def render_sankey_html(tsv_file: str | Path, *, title: str | None = None) -> str:
    sankey_input: SankeyInput = load_sankey_input(tsv_file)
    figure = sankey_from_dataframe(
        sankey_input.frame, level_cols=sankey_input.levels, value_col=sankey_input.value_column, title=title
    )
    return figure.to_html(full_html=True, include_plotlyjs="cdn")


def main(tsv_file: str, output_html: str | None = None, *, title: str | None = None) -> Path:
    output_path = Path(output_html) if output_html else Path(tsv_file).with_suffix(".html")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    html = render_sankey_html(tsv_file, title=title or f"Sankey: {Path(tsv_file).name}")
    output_path.write_text(html, encoding="utf-8")
    return output_path


def build_parser(parser: argparse.ArgumentParser | None = None) -> argparse.ArgumentParser:
    parser = parser or argparse.ArgumentParser(description="Render a Sankey diagram from an EMU abundance table.")
    parser.add_argument("tsv_file", help="EMU abundance TSV file.")
    parser.add_argument("-o", "--output", help="Optional HTML output path. Defaults to replacing the TSV suffix with .html.")
    parser.add_argument("--title", help="Plot title to display in the rendered HTML.")
    return parser


def run(args: argparse.Namespace) -> int:
    destination = main(args.tsv_file, args.output, title=args.title)
    print(f"Sankey plot saved to {destination}")
    return 0


def cli(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(cli())
