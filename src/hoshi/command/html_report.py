from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Sequence
from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.ingress import read_emu_abundance

_REPORT_TEMPLATE_NAME = "basic_report.html.j2"
_TEMPLATE_DIR = Path(__file__).with_name("templates")


def _load_tables(tsv_files: Iterable[Path]) -> list[dict[str, str]]:
    tables: list[dict[str, str]] = []
    for index, path in enumerate(tsv_files):
        df = read_emu_abundance(str(path), reorder = True)
        # Use pandas HTML rendering for quick table output with consistent styling.
        table_html = df.to_html(index=False, border=0, classes="data-table")
        tables.append(
            {
                "name": path.name,
                "html": table_html,
                "dom_id": f"table-{index}",
            }
        )
    return tables


def generate_basic_html_report(
    tsv_files: Sequence[Path],
    *,
    page_title: str = "TSV Report",
) -> str:
    """
    Render a simple HTML report with one table section per TSV file.

    Parameters
    ----------
    tsv_files : Sequence[pathlib.Path]
        Collection of TSV file paths to render.
    page_title : str, optional
        Title to display at the top of the report.

    Returns
    -------
    str
        Rendered HTML document as a string.

    Raises
    ------
    ValueError
        If ``tsv_files`` is empty.
    """
    if not tsv_files:
        raise ValueError("At least one TSV file must be provided.")

    tables = _load_tables(tsv_files)

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2"), default_for_string=True),
    )
    template = env.get_template(_REPORT_TEMPLATE_NAME)
    return template.render(page_title=page_title, tables=tables)


def build_parser(parser: argparse.ArgumentParser | None = None) -> argparse.ArgumentParser:
    parser = parser or argparse.ArgumentParser(description="Render a simple HTML report from TSV files.")
    parser.add_argument(
        "tsv_files",
        nargs="+",
        help="One or more TSV files to include in the report.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to write the generated HTML. Prints to stdout if omitted.",
    )
    parser.add_argument(
        "--title",
        default="TSV Report",
        help="Page title for the generated report (default: 'TSV Report').",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    tsv_files = [Path(path) for path in args.tsv_files]
    html = generate_basic_html_report(tsv_files, page_title=args.title)

    if args.output:
        Path(args.output).write_text(html, encoding="utf-8")
    else:
        print(html)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
