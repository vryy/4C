"""Check that files have correct filenames and includes use correct style"""

import argparse
from baci_utils import common_utils as utils
import os
import re

# The name of files that indicate that a directory should be treated as module root directoy.
module_root_maker_file_name = ".contains_modules"


def find_all_module_roots(folders):
    """
    Recursively step through all directories and return all paths that contain module root marker files.
    """

    def is_module_root(files):
        return any(f == module_root_maker_file_name for f in files)

    return [
        base
        for folder in folders
        for base, dirs, files in os.walk(os.path.abspath(folder))
        if is_module_root(files)
    ]


def is_prefixed_by_module(file, module):
    """
    A valid prefix is equal to the module name followed either by an underscore or a file extension.
    """
    return file.startswith("4C_" + module + "_") or file.startswith(
        "4C_" + module + "."
    )


def most_specific_module_root(file, module_roots):
    return next((m for m in module_roots if file.startswith(m)), None)


def has_valid_filename(path, module_roots):
    """
    Check that the file has a valid filename.

    A valid filename
    - is prefixed by '4C' and the module name
    - is postfixed by '_test' if it is a test file
    """
    abs_path, file = os.path.split(os.path.abspath(path))

    my_module_root = most_specific_module_root(abs_path, module_roots)

    if "tests" in abs_path.split("/") or "unittests" in abs_path.split("/"):
        return (
            file.endswith("_test.cpp") or file.endswith("_test.hpp")
        ) and file.startswith("4C_")

    # If the file is not in a module, we do not require a specific prefix and return true.
    if my_module_root is None:
        return True

    assert abs_path.startswith(os.path.abspath(my_module_root))
    # The module name is the first directory in the file that can be reached from the module_root.
    module = os.path.relpath(abs_path, os.path.abspath(my_module_root)).split(
        os.path.sep
    )[0]
    return is_prefixed_by_module(file, module)


def find_header_include(line):
    """
    Return the included file's name and True if it is included with quotes and False otherwise.
    """
    # short-circuit expensive regex checks if line is clearly not an include directive
    if not line.startswith("#include"):
        return None, None

    match = re.search(r'^#include "(\S+)"', line)
    if match:
        return match[1], True
    match = re.search(r"^#include <(\S+)>", line)
    if match:
        return match[1], False
    else:
        return None, None


def invalid_includes(path, module_names):
    """
    Check that own files are included with quotes, others are included with angle brackets.
    No relative directory navigation (using '.' or '..') is allowed in includes.
    """
    invalid_include_lines = []
    for line_number, line in enumerate(utils.file_contents(path)):
        included_file, included_with_quotes = find_header_include(line)
        if included_file is not None:
            # a file that is part of a module's interface
            is_own_module_file = any(
                is_prefixed_by_module(included_file, m) for m in module_names
            )

            # a file residing in a subdirectory and thus an internal implementation detail
            is_internal_include = not is_own_module_file and os.path.exists(
                os.path.join(os.path.dirname(path), included_file)
            )

            # none of the above, thus an external file
            is_external_file = not is_own_module_file and not is_internal_include

            reason = None
            if "./" in included_file:
                reason = "no relative directories allowed"
            elif is_own_module_file:
                if not included_with_quotes:
                    reason = 'use quotes ""'
            elif is_internal_include:
                if not included_with_quotes:
                    reason = 'use quotes ""'
            elif is_external_file:
                if included_with_quotes:
                    reason = "use angle brackets <>"

            if reason is not None:
                invalid_include_lines.append(
                    "  l "
                    + str(line_number + 1)
                    + ": "
                    + line.strip()
                    + " ("
                    + reason
                    + ")"
                )

    if invalid_include_lines:
        invalid_include_lines.insert(0, "In '" + path + "':")
    return invalid_include_lines


def check_cpp_files_filename(look_cmd, module_roots):
    files_with_wrong_filename = [
        ff
        for ff in utils.files_changed(look_cmd)
        if utils.is_source_file(ff) and not has_valid_filename(ff, module_roots)
    ]
    return files_with_wrong_filename


def check_include_style(look_cmd, module_roots):
    def get_module_names(module_root):
        module_names = [
            os.path.basename(f.path) for f in os.scandir(module_root) if f.is_dir()
        ]
        return module_names

    all_module_names = [i for m in module_roots for i in get_module_names(m)]

    files_with_wrong_include_style = [
        output
        for ff in utils.files_changed(look_cmd)
        if utils.is_source_file(ff)
        for output in invalid_includes(ff, all_module_names)
    ]
    return files_with_wrong_include_style


def main():
    # build command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "paths",
        nargs="+",
        help="The paths in which to check source files for the correct filename.",
    )

    parser.add_argument(
        "--diff_only",
        action="store_true",
        help="Add this tag if only the difference to HEAD should be analyzed. This flag should be used as a pre-commit hook. Otherwise all files are checked.",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Add this tag if the error message should be written to a file.",
    )
    args = parser.parse_args()

    # error file (None for sys.stderr)
    errfile = args.out
    if args.diff_only:
        look_cmd = "git diff --name-only --cached --diff-filter=MRAC -- " + " ".join(
            args.paths
        )
    else:
        look_cmd = "git ls-files " + " ".join(args.paths)

    # Get all the module roots and sort them in reverse.
    # This has the effect that longer paths come before shorter paths. If we match files to the module roots in this
    # sorted order, we can abort at the first match and be sure to have the most specific match.
    module_roots = find_all_module_roots(args.paths)
    module_roots.sort(reverse=True)

    errors_filename = check_cpp_files_filename(look_cmd, module_roots)

    errors_include_style = check_include_style(look_cmd, module_roots)

    utils.pretty_print_error_report(
        "A valid filename looks like this: baci_<module_name>_<detailed_name>[_test].(.hpp|cpp). "
        + "The following files have incorrect filenames:",
        errors_filename,
        errfile,
    )
    utils.pretty_print_error_report(
        "The following files use an invalid style for the #include directive:",
        errors_include_style,
        errfile,
    )

    return len(errors_filename) + len(errors_include_style)


if __name__ == "__main__":
    import sys

    sys.exit(main())
