"4C CODEOWNERS check"

import argparse
from four_c_utils import common_utils as utils
from four_c_utils import patched_code_owners


def check_codeowners_rules_are_used(filenames, allerrors):
    codeowners_file = ".gitlab/CODEOWNERS"
    with open(codeowners_file, "r") as file:
        codeowners_content = file.read()
    owner = patched_code_owners.patched_code_owners(codeowners_content)
    pathlist = list()
    for i in owner.paths:
        pathlist.append(i[1])
    for ff in filenames:
        for i in owner.matching_lines(ff):
            try:
                pathlist.remove(i[2])
            except ValueError:
                ()

    if len(pathlist) > 0:
        if len(allerrors) > 0:
            allerrors.append("")
        allerrors.append(
            "The following CODEOWNER rules do not represent files or directories:"
        )
        allerrors += pathlist
    return len(pathlist)


def check_files_are_owned(filenames, allerrors):
    codeowners_file = ".gitlab/CODEOWNERS"
    with open(codeowners_file, "r") as file:
        codeowners_content = file.read()
    owner = patched_code_owners.patched_code_owners(codeowners_content)
    files_not_in_CO = [
        ff for ff in filenames if owner.matching_line(ff)[2] == "*" and ff != ""
    ]
    if len(files_not_in_CO) > 0:
        if len(allerrors) > 0:
            allerrors.append("")
        allerrors.append("The following files do not exist in the CODEOWNERS file:")
        allerrors += files_not_in_CO
    return len(files_not_in_CO)


def main():
    # build command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="*")
    parser.add_argument(
        "--check-rules",
        action="store_true",
        help="Check whether all rules in CODEOWNERS are used at least once with the set of filenames.",
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
    errors = 0
    allerrors = []
    if args.check_rules:
        errors += check_codeowners_rules_are_used(args.filenames, allerrors)
    errors += check_files_are_owned(args.filenames, allerrors)

    utils.pretty_print_error_report("", allerrors, errfile)
    return errors


if __name__ == "__main__":
    import sys

    sys.exit(main())
