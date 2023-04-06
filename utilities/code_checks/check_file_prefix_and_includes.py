"Check that files are prefixed with module names and includes use correct style"

import argparse
import common_utils as utils
import os
import re


def is_prefixed_by_module(file, module):
    """
    A valid prefix is equal to the module name followed either by an underscore or a file extension.
    """
    return file.startswith(module + "_") or file.startswith(module + ".")


def has_valid_prefix(path, module_root):
    abs_path, file = os.path.split(os.path.abspath(path))

    assert abs_path.startswith(os.path.abspath(module_root))
    # The module name is the first directory in the file that can be reached from the module_root.
    module = os.path.relpath(abs_path, os.path.abspath(module_root)).split(os.path.sep)[
        0
    ]
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


def check_cpp_files_for_prefix(look_cmd, module_root):
    wrongly_prefixed_files = [
        ff
        for ff in utils.files_changed(look_cmd)
        if utils.is_source_file(ff) and not has_valid_prefix(ff, module_root)
    ]
    return wrongly_prefixed_files


def check_include_style(look_cmd, module_root):
    module_names = [os.path.basename(x[0]) for x in os.walk(module_root)]
    # first entry is the module directory itself
    module_names.pop(0)

    files_with_wrong_include_style = [
        output
        for ff in utils.files_changed(look_cmd)
        if utils.is_source_file(ff)
        for output in invalid_includes(ff, module_names)
    ]
    return files_with_wrong_include_style


def main():
    # build command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "path", help="The path in which to check source files for the correct prefix."
    )
    parser.add_argument(
        "module_root",
        help="The root folder that contains all the module folders. Every folder contained inside"
        "module_root will be treated as a valid module name",
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
        look_cmd = "git diff --name-only --cached --diff-filter=MRAC -- " + args.path
    else:
        look_cmd = "git ls-files " + args.path

    allerrors = check_cpp_files_for_prefix(look_cmd, args.module_root)

    allerrors_include = check_include_style(look_cmd, args.module_root)

    utils.pretty_print_error_report(
        "The following files are not prefixed by the module they reside in:",
        allerrors,
        errfile,
    )
    utils.pretty_print_error_report(
        "The following files use an invalid style for the #include directive:",
        allerrors_include,
        errfile,
    )

    return len(allerrors) + len(allerrors_include)


if __name__ == "__main__":
    import sys

    sys.exit(main())
