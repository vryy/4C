import argparse
import re
import enum


class LogSection(enum.Enum):
    SETUP = 0
    TEST_OUTPUT = 1
    SUMMARY = 2


def write_test_output_lines(test_map, file_out):
    for test_id in sorted(test_map.keys()):
        write_lines(test_map[test_id], file_out)


def write_lines(list_of_lines, file_out):
    for line in list_of_lines:
        file_out.write(line)


def main():
    parser = argparse.ArgumentParser()

    # read commmand line arguments
    parser.add_argument(
        "log_file",
        help="Location to log file that was written during pipeline run.",
    )
    parser.add_argument(
        "log_file_sorted",
        help="Location to log file that is sorted by the number of tests.",
    )
    args = parser.parse_args()

    # setup lines before the actual output
    setup_output_lines = []

    # map between test number and text lines of this test
    test_output_lines_map = {}

    # lines could not be recognized (will be printed after the tests)
    unrecognized_output_lines = []

    # summary lines after test output
    summary_output_lines = []

    # define some regex to classify each line
    regex_new_test = re.compile(r"^test\s(\d+)$")
    regex_start_test = re.compile(r"^\s*Start\s+(\d+):")
    regex_test_output_line = re.compile(r"^(\d+):")
    regex_test_finish = re.compile(r"^\s*\d+\/\d+\s+Test\s+#(\d+):\s")

    # loop over all lines of output
    with open(args.log_file, "r") as file_in:

        log_section = LogSection.SETUP
        for line in file_in:
            if log_section == LogSection.SETUP:
                # check whether this line indicates the start of the test suite
                if match_new_test := regex_new_test.match(line):
                    test_id = int(match_new_test.group(1))
                    test_output_lines_map[test_id] = []
                    test_output_lines_map[test_id].append(line)
                    log_section = LogSection.TEST_OUTPUT
                else:
                    setup_output_lines.append(line)

            elif log_section == LogSection.TEST_OUTPUT:
                if (
                    (match_test_line := regex_new_test.match(line)) != None
                    or (match_test_line := regex_start_test.match(line)) != None
                    or (match_test_line := regex_test_output_line.match(line)) != None
                    or (match_test_line := regex_test_finish.match(line)) != None
                ):
                    test_id = int(match_test_line.group(1))
                    if test_id not in test_output_lines_map:
                        test_output_lines_map[test_id] = []

                    test_output_lines_map[test_id].append(line)
                elif line == "\n":
                    continue
                elif line == "The following tests passed:\n":
                    summary_output_lines.append(line)
                    log_section = LogSection.SUMMARY
                else:
                    unrecognized_output_lines.append(line)
                    print(f"Warning: Unrecognized line `{line.strip()}`")

            elif log_section == LogSection.SUMMARY:
                summary_output_lines.append(line)

    with open(args.log_file_sorted, "w") as file_out:
        write_lines(setup_output_lines, file_out)
        write_test_output_lines(test_output_lines_map, file_out)
        write_lines(unrecognized_output_lines, file_out)
        write_lines(summary_output_lines, file_out)


if __name__ == "__main__":
    main()
