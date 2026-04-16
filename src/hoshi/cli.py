from __future__ import annotations

import argparse
import importlib.metadata
import sys
from typing import Sequence

from hoshi.command import contrast
from hoshi.command import enrich
from hoshi.command import report_multi
from hoshi.command import report_single


def _resolve_version() -> str:
    try:
        return importlib.metadata.version("hoshi")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="hoshi",
        description="Hoshi command line interface.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        dest="show_version",
        help="Show the installed hoshi version and exit.",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    subparsers.required = True

    # Register sub-commands from modules
    contrast.build_parser(subparsers)
    enrich.build_parser(subparsers)
    report_multi.build_parser(subparsers)
    report_single.build_parser(subparsers)

    args = parser.parse_args(argv)

    if args.show_version:
        print(f"{parser.prog} {_resolve_version()}")
        return 0

    if hasattr(args, "func"):
        try:
            result = args.func(args)
        except SystemExit as exc:
            code = exc.code
            if isinstance(code, int):
                return code
            if code is None:
                return 0
            print(code, file=sys.stderr)
            return 1

        if isinstance(result, int):
            return result
    else:
        # This path should not be reachable if subparsers are required and set up correctly
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
