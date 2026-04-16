from __future__ import annotations

import argparse
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape

from hoshi.lib.ingress import read_emu_abundance

_REPORT_TEMPLATE_NAME = "multi_report.html.j2"
_TEMPLATE_DIR = Path(__file__).parent / "templates"


def run(args: argparse.Namespace) -> int:
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "j2")),
    )

    reports = []
    for path_str in args.tsv_files:
        path = Path(path_str)
        df = read_emu_abundance(str(path))
        reports.append({"name": path.name, "html": df.to_html(index=False, classes="data-table")})

    template = env.get_template(_REPORT_TEMPLATE_NAME)
    html = template.render(page_title=args.title or "Multi-Sample Report", reports=reports)

    if args.output:
        Path(args.output).write_text(html, encoding="utf-8")
    else:
        print(html)
    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser("report-multi", help="Generate a combined HTML report for multiple samples.")
    parser.add_argument("tsv_files", nargs="+", help="TSV files to include in the report.")
    parser.add_argument("-o", "--output", help="Output HTML file path.")
    parser.add_argument("--title", help="Title for the summary report.")

    parser.set_defaults(func=run)
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    build_parser(subparsers)
    raise SystemExit(run(parser.parse_args()))
