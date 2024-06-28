"""4C header check"""

import argparse
from four_c_utils import file_header as bh
from four_c_utils import common_utils as utils


# CHECK FOR TABS
def contains_tabs(filename):
    "Return True if this version of the file contains tabs."
    return "\t" in utils.file_contents(filename)


def check_support_files_for_tabs(file_names, allerrors):
    "Check support files in this transaction are tab-free."

    support_files_with_tabs = [ff for ff in file_names if contains_tabs(ff)]
    if len(support_files_with_tabs) > 0:
        if len(allerrors) > 0:
            allerrors.append("")
        allerrors.append("The following support files contain tabs:")
        allerrors += support_files_with_tabs
    return len(support_files_with_tabs)


# CHECK HEADER
def check_cpp_files_for_header(file_names, allerrors):
    "Check C/C++ files in this transaction."
    headers = dict([(ff, bh.Header(utils.file_contents(ff))) for ff in file_names])
    # \brief tag
    cpp_files_wo_brief = []
    # check for correct start of header
    cpp_files_wrong_start = [
        ff for ff, hdr in headers.items() if len(hdr.get_start()) < 1
    ]
    if len(cpp_files_wrong_start) > 0:
        if len(allerrors) > 0:
            allerrors.append("")
        allerrors.append(
            "The following files do not start with '/*! \\file' or '/** \\file' as an appropriate header marker:"
        )
        allerrors += cpp_files_wrong_start
    # \level tag
    cpp_files_wo_lvl = [
        ff for ff, hdr in headers.items() if not (0 <= hdr.get_level() <= 3)
    ]
    if len(cpp_files_wo_lvl) > 0:
        if len(allerrors) > 0:
            allerrors.append("")
        allerrors.append("The following files are missing a \\level tag:")
        allerrors += cpp_files_wo_lvl

    # print example header
    if (
        len(cpp_files_wo_brief) > 0
        or len(cpp_files_wrong_start) > 0
        or len(cpp_files_wo_lvl) > 0
    ):
        allerrors += bh.Header.get_example()

    return len(cpp_files_wo_brief) + len(cpp_files_wo_lvl) + len(cpp_files_wrong_start)


#######################################################################################################################


def main():
    # build command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="*")
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Add this tag if the error message should be written to a file.",
    )
    args = parser.parse_args()

    # error file (None for sys.stderr)
    errfile = args.out
    errors = 0
    allerrors = []
    try:
        errors += check_cpp_files_for_header(args.filenames, allerrors)
        errors += check_support_files_for_tabs(args.filenames, allerrors)
    except ValueError:
        print("Something went wrong! Check the error functions in this script again!")
        errors += 1
    utils.pretty_print_error_report("", allerrors, errfile)
    return errors


if __name__ == "__main__":
    import sys

    sys.exit(main())
