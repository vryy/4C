#!/usr/bin/env python3

# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Check changed files for non-ASCII characters."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    errors: list[str] = []

    for filename in sys.argv[1:]:
        path = Path(filename)

        # Skip Markdown files.
        if path.suffix.lower() == ".md":
            continue

        try:
            content = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            # Binary files are ignored, matching the intent of
            # the original grep --binary-files=without-match behavior.
            continue
        except OSError as exc:
            print(f"Could not read {filename}: {exc}")
            return 1

        # Find unique non-ASCII characters.
        non_ascii_chars = sorted({char for char in content if ord(char) > 127})

        for char in non_ascii_chars:
            errors.append(f"{filename} contains non-ascii character '{char}'")

    if errors:
        print()
        print("The following file(s) contain non-ascii characters:")
        print()

        for error in errors:
            print(error)

        print()
        print("--> Please remove those non-ascii characters")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
