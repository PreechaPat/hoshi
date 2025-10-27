from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import sys
from pathlib import Path
from typing import Callable, Sequence

import tomllib

COMMAND_PACKAGE = "hoshi.command"
MainEntry = Callable[..., int | None]

def _discover_commands() -> dict[str, MainEntry]:
    package = importlib.import_module(COMMAND_PACKAGE)
    command_dir = Path(package.__file__).parent
    commands: dict[str, MainEntry] = {}

    for path in sorted(command_dir.glob("*.py")):
        name = path.stem
        if name.startswith("_"):
            continue
        try:
            module = importlib.import_module(f"{COMMAND_PACKAGE}.{name}")
        except ModuleNotFoundError:
            continue
        # Add as module if build_parser and main_entry is avail.
        build_parser = getattr(module, "build_parser", None)
        main_entry = getattr(module, "main", None)
        if not callable(build_parser) or not callable(main_entry):
            continue
        commands[name] = main_entry

    return commands


def _resolve_version() -> str:
    try:
        return importlib.metadata.version("hoshi")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def main(argv: Sequence[str] | None = None) -> int:
    commands = _discover_commands()
    lines = ["Hoshi command line interface."]
    if commands:
        lines.extend(["", "Commands:"])
        lines.extend(f"  {name}" for name in sorted(commands))
    description = "\n".join(lines)
    parser = argparse.ArgumentParser(
        prog="hoshi",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        dest="show_version",
        help="Show the installed hoshi version and exit.",
    )
    parser.add_argument("command", nargs="?", help="Command to execute.")
    parser.add_argument(
        "command_args",
        nargs=argparse.REMAINDER,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    if args.show_version:
        print(f"{parser.prog} {_resolve_version()}")
        return 0

    if not args.command:
        parser.print_help()
        return 0

    entry = commands.get(args.command)
    if entry is None:
        parser.error(f"Unknown command: {args.command}")
    command_name = args.command

    remaining = list(args.command_args)
    if remaining and remaining[0] == "--":
        remaining = remaining[1:]
    original_argv = sys.argv[:]
    sys.argv = [f"{parser.prog} {command_name}", *remaining]
    try:
        code = getattr(entry, "__code__", None)
        defaults = getattr(entry, "__defaults__", ()) or ()
        if code is None:
            accepts_no_args = True
        else:
            positional_args = code.co_argcount
            accepts_no_args = positional_args == 0 or positional_args == len(defaults)
        if accepts_no_args:
            result = entry()
        else:
            result = entry(remaining)
    except SystemExit as exc:
        code = exc.code
        if isinstance(code, int):
            return code
        if code is None:
            return 0
        print(code, file=sys.stderr)
        return 1
    finally:
        sys.argv = original_argv

    if isinstance(result, int):
        return result
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
