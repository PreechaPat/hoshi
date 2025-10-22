from __future__ import annotations

import argparse
from typing import Sequence

from hoshi.lib.ingress import read_emu_abundance


def build_parser(parser: argparse.ArgumentParser | None = None) -> argparse.ArgumentParser:
    parser = parser or argparse.ArgumentParser(description="Compare EMU abundance tables against a negative control.")
    parser.add_argument(
        "-negative",
        required=True,
        help="Path to the negative control EMU abundance table.",
    )
    parser.add_argument(
        "emu_abundance",
        nargs="+",
        help="One or more EMU abundance tables to compare.",
    )
    return parser


def run_emu_compare(args: argparse.Namespace) -> int:
    """CLI entry point for EMU abundance comparison."""
    negative_df = read_emu_abundance(args.negative)
    positive_df = []
    for path in args.emu_abundance:
        pdf = read_emu_abundance(path)
        positive_df.append(pdf)

    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_emu_compare(args)


if __name__ == "__main__":
    raise SystemExit(main())
