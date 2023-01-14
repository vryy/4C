"Check that files are prefixed with module names"

import argparse
import common_utils as utils
import os


def has_valid_prefix(path, module_root):
    """
    A valid prefix is equal to the module name followed either by an underscore or a file extension.
    """
    abs_path, file = os.path.split(os.path.abspath(path))

    assert abs_path.startswith(os.path.abspath(module_root))
    # The module name is the first directory in the file that can be reached from the module_root.
    module = os.path.relpath(abs_path, os.path.abspath(module_root)).split(os.path.sep)[
        0
    ]
    return file.startswith(module + "_") or file.startswith(module + ".")


def check_cpp_files_for_prefix(look_cmd, module_root):
    wrongly_prefixed_files = [
        ff
        for ff in utils.files_changed(look_cmd)
        if utils.is_source_file(ff) and not has_valid_prefix(ff, module_root)
    ]
    return wrongly_prefixed_files


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

    if len(allerrors) > 0:

        allerrors = [
            "The following files are not prefixed by the module they reside in:",
            "",
        ] + allerrors

        if errfile is None:
            utils.pretty_print_error_stderr(allerrors)
        else:
            utils.pretty_print_error_file(allerrors, errfile)
    return len(allerrors)


if __name__ == "__main__":
    import sys

    sys.exit(main())
