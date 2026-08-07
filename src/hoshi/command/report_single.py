from __future__ import annotations

import argparse
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.diversity import compute_diversity
from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.ingress import emu_to_experiment
from hoshi.lib.sankey import get_sankey_data, render_sankey_figure

_REPORT_TEMPLATE_NAME = "single_report.html.j2"
_TEMPLATE_DIR = Path(__file__).parent / "templates"


def _find_emu_files(sample_dir: Path) -> dict[str, Path | None]:
    """Discover EMU output files in a sample directory.

    Returns a dict with keys: 'abundance', 'read_assignments', 'unclassified', 'unmapped'.
    Values are Path objects or None if not found.
    """
    files: dict[str, Path | None] = {
        "abundance": None,
        "read_assignments": None,
        "unclassified": None,
        "unmapped": None,
    }

    for f in sample_dir.iterdir():
        if not f.is_file():
            continue
        name = f.name
        if name.endswith("_rel-abundance.tsv"):
            files["abundance"] = f
        elif name.endswith("_read-assignment-distributions.tsv"):
            files["read_assignments"] = f
        elif name.endswith("_unclassified_mapped.fastq.gz"):
            files["unclassified"] = f
        elif name.endswith("_unmapped.fastq.gz"):
            files["unmapped"] = f

    return files


def generate_single_html_report(
    se: SummarizedExperiment,
    *,
    page_title: str | None = None,
    sample_name: str | None = None,
) -> str:
    """Generate an HTML report from a SummarizedExperiment (single sample).

    Parameters
    ----------
    se : SummarizedExperiment
        A single-sample experiment container with assays and row_data.
    page_title : str, optional
        Override the full page title.
    sample_name : str, optional
        Sample name for the report header. Falls back to SE's sample_id.
    """
    # Determine the sample to report on
    name = sample_name or str(se.sample_ids[0])
    title = page_title or name

    # Reconstruct flat DataFrame from SE (for downstream functions)
    df = se.to_dataframe()

    # Compute diversity statistics from the full DataFrame
    stats = compute_diversity(df)

    # Filter out meta rows for the display table, sort by abundance descending
    meta_ids = {"unmapped", "mapped_filtered", "mapped_unclassified"}
    df_display = df[~df["tax_id"].astype(str).isin(meta_ids)].copy()
    df_display = df_display.sort_values("abundance", ascending=False).reset_index(drop=True)
    df_subset = df_display[["species", "tax_id", "abundance", "estimated counts"]]

    # Top-5 and Top-10 tables for the report
    top5_html = df_subset.head(5).to_html(index=False, border=0, classes="data-table")
    top10_html = df_subset.head(10).to_html(index=False, border=0, classes="data-table")

    # Generate Sankey plot HTML from the flat DataFrame
    try:
        sankey_data = get_sankey_data(df)
        fig = render_sankey_figure(sankey_data, title=f"Taxonomy Flow: {name}")
        sankey_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    except Exception as e:
        sankey_html = f"<p>Sankey plot could not be generated: {str(e)}</p>"

    # Source file info from col_data if available
    source_file = ""
    sample_dir = ""
    if not se.col_data.empty and "source_file" in se.col_data.columns:
        source_file = se.col_data.iloc[0]["source_file"]
        sample_dir = str(Path(source_file).parent.resolve())

    report_data = {
        "name": name,
        "source_path": source_file,
        "sample_dir": sample_dir,
        "top5_html": top5_html,
        "top10_html": top10_html,
        "sankey_html": sankey_html,
    }

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2"), default_for_string=True),
    )
    template = env.get_template(_REPORT_TEMPLATE_NAME)
    return template.render(page_title=title, table=report_data, stats=stats)


def run(args: argparse.Namespace) -> int:
    sample_dir = Path(args.sample_dir)

    if not sample_dir.is_dir():
        raise NotADirectoryError(f"Expected a sample directory, got: {sample_dir}")

    # Discover Emu files
    emu_files = _find_emu_files(sample_dir)
    abundance_file = emu_files["abundance"]

    if abundance_file is None:
        raise FileNotFoundError(
            f"No *_rel-abundance.tsv file found in {sample_dir}. "
            "Is this an EMU output directory?"
        )

    # Pipeline: Emu TSV → SummarizedExperiment → generate report
    name = args.name or sample_dir.name
    se = emu_to_experiment(abundance_file, sample_names=[name])

    html = generate_single_html_report(se, page_title=args.title, sample_name=name)

    if args.output:
        Path(args.output).write_text(html, encoding="utf-8")
    else:
        print(html)
    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "report-single",
        help="Generate a single-sample HTML report from an EMU output folder.",
        description="Generate an HTML microbiome report from an EMU output directory. "
        "Automatically discovers *_rel-abundance.tsv and other EMU files in the folder.",
    )
    parser.add_argument(
        "sample_dir",
        help="Path to the EMU sample output directory (contains *_rel-abundance.tsv).",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to write the generated HTML. Prints to stdout if omitted.",
    )
    parser.add_argument(
        "-n",
        "--name",
        help="Sample name displayed in the report header (e.g. 'sample01'). "
        "Defaults to the folder name.",
    )
    parser.add_argument(
        "--title",
        help="Override the full page title (takes precedence over --name).",
    )

    parser.set_defaults(func=run)
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    build_parser(subparsers)
    raise SystemExit(run(parser.parse_args()))
