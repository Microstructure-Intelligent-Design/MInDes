#!/usr/bin/env python3
"""Expand MInDes INCLUDE directives into one input file."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import sys
import tempfile


INCLUDE_LINE = re.compile(r"^\s*INCLUDE(?:\s+(.*?))?\s*$")


class PreprocessError(Exception):
    pass


def parse_include(line: str, source: Path, line_number: int) -> str | None:
    match = INCLUDE_LINE.match(line)
    if match is None:
        return None

    value = match.group(1)
    if not value:
        raise PreprocessError(f"{source}:{line_number}: INCLUDE requires a file path")

    if value.startswith(('"', "'")):
        quote = value[0]
        if len(value) < 2 or value[-1] != quote:
            raise PreprocessError(f"{source}:{line_number}: unterminated quoted INCLUDE path")
        value = value[1:-1]
    elif value.endswith(('"', "'")):
        raise PreprocessError(f"{source}:{line_number}: malformed quoted INCLUDE path")

    if not value:
        raise PreprocessError(f"{source}:{line_number}: INCLUDE path cannot be empty")
    return value


def canonical_file(path: Path, included_from: Path | None = None) -> Path:
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        context = f" included from {included_from}" if included_from else ""
        raise PreprocessError(f"cannot read input file {path}{context}: {error}") from error
    if not resolved.is_file():
        raise PreprocessError(f"input path is not a file: {resolved}")
    return resolved


def expand_file(
    path: Path,
    output_lines: list[str],
    visited: set[Path],
    active_stack: list[Path],
) -> bool:
    resolved = canonical_file(path, active_stack[-1] if active_stack else None)
    if resolved in active_stack:
        cycle_start = active_stack.index(resolved)
        chain = active_stack[cycle_start:] + [resolved]
        raise PreprocessError("circular INCLUDE: " + " -> ".join(map(str, chain)))
    if resolved in visited:
        return False

    visited.add(resolved)
    active_stack.append(resolved)
    found_include = False
    try:
        try:
            lines = resolved.read_text(encoding="utf-8-sig").splitlines()
        except (OSError, UnicodeError) as error:
            raise PreprocessError(f"cannot decode {resolved} as UTF-8: {error}") from error

        for line_number, line in enumerate(lines, start=1):
            include_value = parse_include(line, resolved, line_number)
            if include_value is None:
                output_lines.append(line)
                continue
            found_include = True
            include_path = Path(include_value.replace("\\", "/"))
            if not include_path.is_absolute():
                include_path = resolved.parent / include_path
            found_include = expand_file(include_path, output_lines, visited, active_stack) or found_include
    finally:
        active_stack.pop()

    return found_include


def write_atomically(output: Path, content: str) -> None:
    output_parent = output.parent.resolve()
    if not output_parent.is_dir():
        raise PreprocessError(f"output directory does not exist: {output_parent}")

    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=output_parent,
            prefix=output.name + ".",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary.write(content)
            temporary_name = temporary.name
        os.replace(temporary_name, output)
    except OSError as error:
        if temporary_name:
            try:
                Path(temporary_name).unlink(missing_ok=True)
            except OSError:
                pass
        raise PreprocessError(f"cannot write output file {output}: {error}") from error


def preprocess(input_path: Path, output_path: Path) -> bool:
    output_lines: list[str] = []
    found_include = expand_file(input_path, output_lines, set(), [])
    if not found_include:
        return False
    content = "\n".join(output_lines)
    if output_lines:
        content += "\n"
    write_atomically(output_path, content)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    try:
        if not preprocess(arguments.input, arguments.output):
            return 10
    except PreprocessError as error:
        print(f"MInDes input preprocessing error: {error}", file=sys.stderr)
        return 20
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
