from __future__ import annotations

import argparse
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.ingress import read_emu_abundance
from hoshi.lib.sankey import get_sankey_data, render_sankey_figure

_REPORT_TEMPLATE_NAME = "single_report.html.j2"
_TEMPLATE_DIR = Path(__file__).parent / "templates"


def _load_table_data(path: Path) -> dict[str, str]:
    df = read_emu_abundance(str(path), reorder=True)
    table_html = df.to_html(index=False, border=0, classes="data-table")

    try:
        sankey_data = get_sankey_data(path)
        fig = render_sankey_figure(sankey_data, title=f"Sankey: {path.name}")
        sankey_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    except Exception as e:
        sankey_html = f"<p>Sankey plot could not be generated: {str(e)}</p>"

    return {
        "name": path.name,
        "html": table_html,
        "sankey_html": sankey_html,
    }


def generate_single_html_report(tsv_file: Path, *, page_title: str | None = None) -> str:
    if not tsv_file.exists():
        raise FileNotFoundError(f"TSV file not found: {tsv_file}")

    table_data = _load_table_data(tsv_file)
    title = page_title or f"Report: {tsv_file.name}"

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2"), default_for_string=True),
    )
    template = env.get_template(_REPORT_TEMPLATE_NAME)
    return template.render(page_title=title, table=table_data)


def run(args: argparse.Namespace) -> int:
    tsv_file = Path(args.tsv_file)
    html = generate_single_html_report(tsv_file, page_title=args.title)

    if args.output:
        Path(args.output).write_text(html, encoding="utf-8")
    else:
        print(html)
    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser("report-single", help="Render a single HTML report from one TSV file.")
    parser.add_argument("tsv_file", help="The TSV file to report on.")
    parser.add_argument(
        "-o",
        "--output",
        help="Path to write the generated HTML. Prints to stdout if omitted.",
    )
    parser.add_argument(
        "--title",
        help="Page title for the generated report.",
    )

    parser.set_defaults(func=run)
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    build_parser(subparsers)
    raise SystemExit(run(parser.parse_args()))
