from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import pkgutil
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Callable, Sequence

import tomllib

COMMAND_PACKAGE = "hoshi.command"

CommandCallable = Callable[[Sequence[str] | None], int | None]


@dataclass(frozen=True)
class CommandSpec:
    canonical_name: str
    module_name: str
    entrypoint: CommandCallable
    help_text: str | None


def _select_entrypoint(module: ModuleType) -> CommandCallable | None:
    for attr in ("cli", "main"):
        func = getattr(module, attr, None)
        if callable(func):
            return func
    return None


def _derive_help_text(module: ModuleType) -> str | None:
    doc = (module.__doc__ or "").strip()
    if doc:
        return doc.splitlines()[0]
    builder = getattr(module, "build_parser", None)
    if callable(builder):
        parser = builder(argparse.ArgumentParser(add_help=False))
        return parser.description
    return None


def _discover_commands() -> tuple[dict[str, CommandSpec], dict[str, CommandSpec]]:
    package = importlib.import_module(COMMAND_PACKAGE)
    canonical: dict[str, CommandSpec] = {}
    aliases: dict[str, CommandSpec] = {}
    for module_info in pkgutil.iter_modules(package.__path__):
        if module_info.ispkg:
            continue
        module_name = module_info.name
        if module_name.startswith("_"):
            continue
        module = importlib.import_module(f"{COMMAND_PACKAGE}.{module_name}")
        entrypoint = _select_entrypoint(module)
        if entrypoint is None:
            continue
        help_text = _derive_help_text(module)
        canonical_name = module_name.replace("_", "-")
        spec = CommandSpec(
            canonical_name=canonical_name,
            module_name=module_name,
            entrypoint=entrypoint,
            help_text=help_text,
        )
        canonical[canonical_name] = spec
        aliases[canonical_name] = spec
        aliases[module_name] = spec
        for alias in getattr(module, "ALIASES", ()):
            aliases[alias] = spec
    return canonical, aliases


def _build_description(commands: dict[str, CommandSpec]) -> str:
    if not commands:
        return "Hoshi command line interface."
    width = max(len(name) for name in commands)
    lines = ["Hoshi command line interface.", "", "Available commands:"]
    for name in sorted(commands):
        help_text = commands[name].help_text
        line = f"  {name.ljust(width)}"
        if help_text:
            line += f"  {help_text}"
        lines.append(line)
    return "\n".join(lines)


def _resolve_version() -> str:
    try:
        return importlib.metadata.version("hoshi")
    except importlib.metadata.PackageNotFoundError:
        project_root = Path(__file__).resolve().parents[2]
        pyproject = project_root / "pyproject.toml"
        if pyproject.exists():
            with pyproject.open("rb") as handle:
                data = tomllib.load(handle)
            version = data.get("project", {}).get("version")
            if version:
                return version
        return "unknown"


def main(argv: Sequence[str] | None = None) -> int:
    commands, aliases = _discover_commands()
    parser = argparse.ArgumentParser(
        prog="hoshi",
        description=_build_description(commands),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    parser.add_argument(
        "-h",
        "--help",
        action="store_true",
        dest="show_help",
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        dest="show_version",
        help="Show the installed hoshi version and exit.",
    )
    parser.add_argument("command", nargs="?", help="Command to execute.")
    args, remaining = parser.parse_known_args(argv)
    remaining = list(remaining)

    if args.show_version:
        print(f"{parser.prog} {_resolve_version()}")
        return 0

    if args.show_help and args.command is None:
        parser.print_help()
        return 0

    if args.command is None:
        parser.print_help()
        return 0

    command_key = args.command
    spec = aliases.get(command_key) or aliases.get(command_key.replace("-", "_"))
    if spec is None:
        parser.error(f"Unknown command: {command_key}")

    if args.show_help:
        remaining.insert(0, "--help")

    original_prog = sys.argv[0]
    sys.argv[0] = f"{parser.prog} {spec.canonical_name}".strip()
    try:
        result = spec.entrypoint(remaining)
    except SystemExit as exc:
        code = exc.code
        if isinstance(code, int):
            return code
        if code is None:
            return 0
        print(code, file=sys.stderr)
        return 1
    finally:
        sys.argv[0] = original_prog
    if isinstance(result, int):
        return result
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
