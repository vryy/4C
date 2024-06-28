from codecs import namereplace_errors
import os
import sys
import argparse
import re
import json
from four_c_utils import common_utils as utils
from typing import Optional, List


# CHECK INPUT FILE TESTS
def check_inputtests(filenames, allerrors):
    errors = 0

    # read TestingFrameworkListOfTests.cmake
    with open("TestingFrameworkListOfTests.cmake", "r") as cmakefile:
        all_lines = "\n".join(cmakefile.readlines())

    # check if some input tests are missing
    missing_input_tests = []
    for input_test in filenames:
        # check, whether this input file is in TestingFrameworkListOfTests.cmake
        expected_test_name = os.path.splitext(os.path.basename(input_test))[0]
        if re.search(r"\b" + re.escape(expected_test_name) + r"\b", all_lines) is None:
            missing_input_tests.append(input_test)

    # Currently disabled since this needs a few fixes
    if False and len(missing_input_tests) > 0:
        errors += 1
        allerrors.append(
            "The following input files are missing in TestingFrameworkListOfTests.cmake:"
        )
        allerrors.append("")
        allerrors.extend(missing_input_tests)

    # check if input tests have empty sections
    tests_empty_sections = []

    for input_test in filenames:
        with open(input_test, "r") as f:
            num_current_section_non_empty_lines = None

            for line in f:
                if line.startswith("--"):
                    if num_current_section_non_empty_lines == 0:
                        tests_empty_sections.append(input_test)
                        break
                    else:
                        num_current_section_non_empty_lines = 0

                elif num_current_section_non_empty_lines is None:
                    # No section title until now
                    continue

                elif line.strip() != "":
                    num_current_section_non_empty_lines += 1

    if len(tests_empty_sections) > 0:
        errors += 1
        allerrors.append(
            "The following input files have empty sections. Please delete them or correct your input file."
        )
        allerrors.append("")
        allerrors.extend(tests_empty_sections)

    return errors


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
    # check input file tests
    errors += check_inputtests(args.filenames, allerrors)

    utils.pretty_print_error_report("", allerrors, errfile)
    return errors


if __name__ == "__main__":
    import sys

    sys.exit(main())
