from __future__ import annotations

import argparse

from hoshi.lib.ingress import read_emu_abundance


def run(args: argparse.Namespace) -> int:
    negative_df = read_emu_abundance(args.negative)
    positive_df = [read_emu_abundance(path) for path in args.emu_abundance]

    # Logic goes here
    return 0


def build_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "contrast",
        help="Contrast EMU abundance tables against a negative control.",
    )

    parser.add_argument(
        "-negative",
        required=True,
        help="Path to the negative control EMU abundance table.",
    )
    parser.add_argument(
        "emu_abundance",
        nargs="+",
        help="One or more EMU abundance tables to contrast.",
    )

    parser.set_defaults(func=run)
    return parser


if __name__ == "__main__":
    # Minimal setup for direct file execution/testing
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    build_parser(subparsers)
    raise SystemExit(run(parser.parse_args()))
