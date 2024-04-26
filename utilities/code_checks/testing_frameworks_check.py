from codecs import namereplace_errors
import dataclasses
import os, sys, argparse, re, json
from baci_utils import common_utils as utils
from typing import Optional, List


@dataclasses.dataclass(frozen=True)
class TestMacro:
    #  regular expression matching all test names and mpi ranks
    pattern: re.Pattern

    # group ids in the regular expression that matches a test name
    test_names_groups: List[int] = dataclasses.field(default_factory=lambda: [0])

    # group id in the regular expression that matches the mpi rank
    mpirank_group: Optional[int] = None


# CHECK INPUT FILE TESTS
def check_inputtests(look_cmd, allerrors):
    errors = 0
    input_tests = [
        str(ff) for ff in utils.files_changed(look_cmd)[:-1] if utils.is_input_file(ff)
    ]

    mpi_non_compliant_tests: List[str] = []
    list_of_all_testnames: List[str] = []

    # read TestingFrameworkListOfTests.cmake
    with open("TestingFrameworkListOfTests.cmake", "r") as cmakefile:

        all_lines = "\n".join(cmakefile.readlines())

        # Defining the list of test macros with regular expressions that recognize them in the cmake file
        test_macros = [
            TestMacro(
                re.compile(r"baci_test\s*\(*([a-zA-Z0-9_\.\-]+)\s+(\d+)"),
                mpirank_group=1,
            ),
            TestMacro(
                re.compile(
                    r"baci_test_extended_timeout\s*\(*([a-zA-Z0-9_\.\-]+)\s+(\d+)"
                ),
                mpirank_group=1,
            ),
            TestMacro(
                re.compile(r"baci_omp_test\s*\(*([a-zA-Z0-9_\.\-]+)\s+(\d+)"),
                mpirank_group=1,
            ),
            TestMacro(
                re.compile(
                    r"baci_test_and_post_ensight_test\s*\(*([a-zA-Z0-9_\.\-]+)\s+(\d+)"
                ),
                mpirank_group=1,
            ),
            TestMacro(
                re.compile(
                    r"baci_test_restartonly\s*\(*([a-zA-Z0-9_\.\-]+)\s+(?:[a-zA-Z0-9_\.\-]+)\s+(\d+)"
                ),
                mpirank_group=1,
            ),
            TestMacro(
                re.compile(
                    r"baci_test_Nested_Par\s*\(\s*([a-zA-Z0-9_\.\-]+)\s+([a-zA-Z0-9_\.\-]+)\s*"
                ),
                test_names_groups=[0, 1],
            ),
        ]

        for test_macro in test_macros:
            # search for tests of this macro
            for test in test_macro.pattern.findall(all_lines):
                # get list of test names defined in this macro (ignore empty test names)
                my_names = [
                    test[i]
                    for i in test_macro.test_names_groups
                    if len(test[i].strip()) > 0 and test[i].strip() != '""'
                ]
                list_of_all_testnames.extend(my_names)

                # check mpi rank
                if test_macro.mpirank_group is not None:
                    mpi_rank = int(test[test_macro.mpirank_group])

                    if mpi_rank > 3:
                        for name in my_names:
                            mpi_non_compliant_tests.append(name)

    if len(mpi_non_compliant_tests) > 0:
        errors += 1
        allerrors.append("The following tests use an unjustified high mpi rank:")
        allerrors.append("")
        allerrors.extend(mpi_non_compliant_tests)

    # check if some input tests are missing
    missing_input_tests = []
    for input_test in input_tests:
        # check, whether this input file is in TestingFrameworkListOfTests.cmake

        expected_test_name = os.path.splitext(os.path.basename(input_test))[0]
        if expected_test_name not in list_of_all_testnames:
            missing_input_tests.append(input_test)

    if len(missing_input_tests) > 0:
        errors += 1
        allerrors.append(
            "The following input files are missing in TestingFrameworkListOfTests.cmake:"
        )
        allerrors.append("")
        allerrors.extend(missing_input_tests)

    # check if input tests have empty sections
    tests_empty_sections = []

    for input_test in input_tests:
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

    # flag, whether only touched files should be checked
    diff_only = args.diff_only

    # error file (None for sys.stderr)
    errfile = args.out
    errors = 0
    allerrors = []
    try:
        if diff_only:
            look_cmd = "git diff --name-only --cached --diff-filter=MRAC"
        else:
            look_cmd = "git ls-files"

        # check input file tests
        errors += check_inputtests(look_cmd, allerrors)
    except ValueError:
        print("Something went wrong! Check the error functions in this script again!")
        errors += 1
    utils.pretty_print_error_report("", allerrors, errfile)
    return errors


if __name__ == "__main__":
    import sys

    sys.exit(main())
