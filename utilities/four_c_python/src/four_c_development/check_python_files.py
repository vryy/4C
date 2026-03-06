#!/usr/bin/env python3
# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check Python file locations and naming conventions."""

import sys
from pathlib import Path
from four_c_common_utils import common_utils as utils


def valid_python_file_name(file: Path) -> bool:
    """Validate location and naming of Python file."""

    # check 1: all python files in utilities/four_c_python/ are allowed
    if "utilities/four_c_python/" in str(file):
        return True

    # check 2: all python files in /tests/input_files/ are allowed
    if "tests/input_files/" in str(file):
        return True

    # check 3: all python files in src/module/python/ are allowed (for py4C)
    if len(file.parts) > 3 and file.parts[0] == "src" and file.parts[2] == "python":
        module_name = utils.get_module_name(file)

        if module_name is None:
            return False

        expected_prefix = f"4C_{module_name}_"

        return file.name.startswith(expected_prefix) and file.name.endswith("_test.py")

    # check 4: all python files in src/ must be named 4C_<module>_*.py, where <module> is the name of the module they are located in
    if str(file).startswith("src/"):
        module_name = utils.get_module_name(file)

        if module_name is None:
            return False

        expected_prefix = f"4C_{module_name}_"

        return file.name.startswith(expected_prefix)

    return False


def main():
    """Check that all Python files adhere to our naming conventions."""

    wrong_files = []

    for file in sys.argv[1:]:
        if not valid_python_file_name(Path(file)):
            wrong_files.append(file)

    if wrong_files:
        print("The following Python files violate location or naming rules:\n")
        for file in wrong_files:
            print(f"    {file}")

        print(
            "\nThe following rules apply:\n"
            "    - Python files for development (pre-commit hook, ...) should be located in 'utilities/four_c_python' /\n"
            "    - Python files inside modules under src/ must be named according to their module name: 4C_<module>_*.py\n"
            "    - Python files inside modules under python/ must be named according to their module name: 4C_<module>_*_test.py\n"
            "    - Python files in /tests/input_files/ are allowed (e.g., for testing purposes) and must be utilized/typed within input files\n"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
