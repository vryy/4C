#!/usr/bin/env python3

# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Check commit messages against the 4C commit message guidelines."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


IMPERATIVE_MOOD_BLACKLIST = {
    "added",
    "adds",
    "adding",
    "affixed",
    "affixes",
    "affixing",
    "adjusted",
    "adjusts",
    "adjusting",
    "amended",
    "amends",
    "amending",
    "avoided",
    "avoids",
    "avoiding",
    "bumped",
    "bumps",
    "bumping",
    "changed",
    "changes",
    "changing",
    "checked",
    "checks",
    "checking",
    "committed",
    "commits",
    "committing",
    "copied",
    "copies",
    "copying",
    "corrected",
    "corrects",
    "correcting",
    "created",
    "creates",
    "creating",
    "decreased",
    "decreases",
    "decreasing",
    "deleted",
    "deletes",
    "deleting",
    "disabled",
    "disables",
    "disabling",
    "dropped",
    "drops",
    "dropping",
    "duplicated",
    "duplicates",
    "duplicating",
    "enabled",
    "enables",
    "enabling",
    "enhanced",
    "enhances",
    "enhancing",
    "excluded",
    "excludes",
    "excluding",
    "extracted",
    "extracts",
    "extracting",
    "fixed",
    "fixes",
    "fixing",
    "handled",
    "handles",
    "handling",
    "implemented",
    "implements",
    "implementing",
    "improved",
    "improves",
    "improving",
    "included",
    "includes",
    "including",
    "increased",
    "increases",
    "increasing",
    "installed",
    "installs",
    "installing",
    "introduced",
    "introduces",
    "introducing",
    "leased",
    "leases",
    "leasing",
    "managed",
    "manages",
    "managing",
    "merged",
    "merges",
    "merging",
    "moved",
    "moves",
    "moving",
    "normalised",
    "normalises",
    "normalising",
    "normalized",
    "normalizes",
    "normalizing",
    "passed",
    "passes",
    "passing",
    "pointed",
    "points",
    "pointing",
    "pruned",
    "prunes",
    "pruning",
    "ran",
    "runs",
    "running",
    "refactored",
    "refactors",
    "refactoring",
    "released",
    "releases",
    "releasing",
    "removed",
    "removes",
    "removing",
    "renamed",
    "renames",
    "renaming",
    "replaced",
    "replaces",
    "replacing",
    "resolved",
    "resolves",
    "resolving",
    "reverted",
    "reverts",
    "reverting",
    "sets",
    "setting",
    "showed",
    "shows",
    "showing",
    "swapped",
    "swaps",
    "swapping",
    "tested",
    "tests",
    "testing",
    "tidied",
    "tidies",
    "tidying",
    "updated",
    "updates",
    "updating",
    "upped",
    "ups",
    "upping",
    "used",
    "uses",
    "using",
}


URL_REGEX = re.compile(
    r"^\s*"
    r"(https?|ftp|file|wss?|git|ssh|data|irc|dat)://"
    r"[-A-Za-z0-9+&@#/%?=~_|!:,.;]*"
    r"[-A-Za-z0-9+&@#/%=~_|]"
)


CUT_LINE = "# ------------------------ >8 ------------------------"


def get_editor() -> str:
    """Get the editor configured for Git."""

    commands = [
        ["git", "config", "--get", "core.editor"],
    ]

    for command in commands:
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError:
            continue

        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()

    return (
        # Git normally uses these environment variables as fallbacks.
        __import__("os").environ.get("VISUAL")
        or __import__("os").environ.get("EDITOR")
        or "vi"
    )


def read_commit_message(commit_msg_file: Path) -> list[str]:
    """Read and clean the commit message."""

    lines: list[str] = []

    # utf-8-sig also handles a UTF-8 BOM if one is present.
    content = commit_msg_file.read_text(
        encoding="utf-8-sig",
        errors="replace",
    )

    for line in content.splitlines():
        # Remove trailing spaces.
        line = line.rstrip(" ")

        # Ignore everything after Git's cut line.
        if line == CUT_LINE:
            break

        # Ignore comments.
        if line.startswith("#"):
            continue

        lines.append(line)

    return lines


def add_warning(
    warnings: dict[int, list[str]],
    line_number: int,
    message: str,
) -> None:
    """Add a warning to a line."""

    warnings.setdefault(line_number, []).append(message)


def validate_commit_message(
    lines: list[str],
) -> dict[int, list[str]]:
    """Validate the commit message."""

    warnings: dict[int, list[str]] = {}

    if not lines:
        return warnings

    # Capture the subject and remove special prefixes.
    subject = lines[0]

    for prefix in ("squash! ", "fixup! ", "amend! ", "reword! "):
        if subject.startswith(prefix):
            subject = subject[len(prefix) :]

    # If the commit is effectively empty, there is nothing to validate.
    if not any(line.strip() for line in lines):
        return warnings

    # ------------------------------------------------------------------
    # 1. Separate subject from body with a blank line
    # ------------------------------------------------------------------

    if len(lines) >= 2 and lines[1] != "":
        add_warning(
            warnings,
            2,
            "Separate subject from body with a blank line",
        )

    # ------------------------------------------------------------------
    # 2. Limit the subject line to 50 characters
    # ------------------------------------------------------------------

    if len(subject) > 50:
        add_warning(
            warnings,
            1,
            f"Limit the subject line to 50 characters ({len(subject)} chars)",
        )

    # ------------------------------------------------------------------
    # 3. Capitalize the subject line
    # ------------------------------------------------------------------

    if not re.match(r"^\s*[A-Z]", subject):
        add_warning(
            warnings,
            1,
            "Capitalize the subject line",
        )

    # ------------------------------------------------------------------
    # 4. Do not end the subject line with a period
    # ------------------------------------------------------------------

    if not re.search(r"[^.]$", subject):
        add_warning(
            warnings,
            1,
            "Do not end the subject line with a period",
        )

    # ------------------------------------------------------------------
    # 5. Use the imperative mood in the subject line
    # ------------------------------------------------------------------

    first_word_match = re.match(r"^\s*(\S+)", subject)

    if first_word_match:
        first_word = first_word_match.group(1).lower()

        if first_word in IMPERATIVE_MOOD_BLACKLIST:
            add_warning(
                warnings,
                1,
                "Use the imperative mood in the subject line, " "e.g 'fix' not 'fixes'",
            )

    # ------------------------------------------------------------------
    # 6. Wrap the body at 72 characters
    # ------------------------------------------------------------------

    for index, line in enumerate(lines):
        line_number = index + 1

        if len(line) > 72 and not URL_REGEX.match(line):
            add_warning(
                warnings,
                line_number,
                f"Wrap the body at 72 characters ({len(line)} chars)",
            )

    # ------------------------------------------------------------------
    # 7. Use the body to explain what and why vs. how
    # ------------------------------------------------------------------

    # Intentionally not implemented in the original script.

    # ------------------------------------------------------------------
    # 8. Do not write single-worded commits
    # ------------------------------------------------------------------

    subject_words = subject.split()

    if len(subject_words) <= 1:
        add_warning(
            warnings,
            1,
            "Do no write single worded commits",
        )

    # ------------------------------------------------------------------
    # 9. Do not start the subject line with whitespace
    # ------------------------------------------------------------------

    if re.match(r"^\s+", subject):
        add_warning(
            warnings,
            1,
            "Do not start the subject line with whitespace",
        )

    return warnings


def display_warnings(
    warnings: dict[int, list[str]],
    lines: list[str],
    commit_msg_file: Path,
) -> None:
    """Display commit message warnings."""

    print()
    print("Your commit message is not following the commit message guidelines:")
    print()

    for line_number in sorted(warnings):
        line_index = line_number - 1
        line = lines[line_index] if 0 <= line_index < len(lines) else ""

        print(f"{line:<74} [line {line_number}]")

        for warning in warnings[line_number]:
            print(f" - {warning}")

    print()
    print(f"Your commit message is saved in {commit_msg_file}")
    print(
        "To retry the commit and edit that message, append the following "
        "to the git commit command:"
    )
    print()
    print(f"  -F {commit_msg_file} --edit")
    print()


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: check_commit_msg.py <commit-msg-file>")
        return 1

    commit_msg_file = Path(sys.argv[1])

    if not commit_msg_file.is_file():
        print(f"Commit message file does not exist: {commit_msg_file}")
        return 1

    try:
        lines = read_commit_message(commit_msg_file)
    except OSError as exc:
        print(f"Could not read commit message: {exc}")
        return 1

    warnings = validate_commit_message(lines)

    if not warnings:
        return 0

    display_warnings(
        warnings,
        lines,
        commit_msg_file,
    )

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
