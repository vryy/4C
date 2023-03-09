import argparse
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    # read_input
    parser.add_argument(
        "number_of_test",
        help="Total number of tests.",
    )
    parser.add_argument(
        "log_file",
        help="Location to log file that was written during pipeline run.",
    )
    parser.add_argument(
        "log_file_sorted",
        help="Location to log file that is sorted by the number of tests.",
    )
    args = parser.parse_args()

    # flag, whether only touched files should be checked
    num_tests = int(args.number_of_test)
    log_file = os.path.abspath(args.log_file)
    log_file_sorted = os.path.abspath(args.log_file_sorted)

    # map between test number (index+1) and text lines of this test
    test_map = [[] for x in range(num_tests)]

    # lines with text from setup (i.e. configuration, compilation, etc.)
    setup_lines = []

    # lines of summary after test
    summary_lines = []

    # first test output. Everything before this line is setup output.
    first_test_line = 10000000000000

    # loop over all lines of output
    with open(log_file, "r") as file_in:
        for line_num, line in enumerate(file_in):
            # check if first characters of line contain ": ". This indicates a test output
            check_string = line[0:6]
            if check_string.find(": ") != -1:
                check_string_split = line[0:6].split(": ")[0]
                if check_string_split.isnumeric():
                    test_map[int(check_string_split) - 1].append(line)
                    first_test_line = line_num
            # special cases
            else:
                found_test_output = False
                # check for "test xxx". This indicates a new test
                for test_num in range(1, num_tests + 1):
                    if line.find("test " + str(test_num) + "\n") != -1:
                        test_map[test_num - 1].append(line)
                        found_test_output = True
                        first_test_line = line_num
                        break
                # check for "Start". This indicates the start of a test
                if line.find("      Start ") != -1:
                    check_string_2 = line.split(":")[0].split()[1]
                    if check_string_2.isnumeric():
                        test_map[int(check_string_2) - 1].append(line)
                        found_test_output = True
                        first_test_line = line_num
                # check for "/num_tests Test". This indicates the end of a test
                elif line.find(" Test ") != -1 and (
                    line.find("   Passed   ") != -1
                    or line.find("***Failed") != -1
                    or line.find("***Not Run") != -1
                    or line.find("***Timeout") != -1
                ):
                    check_string_3 = line.split("#")[1].split(":")[0]
                    if check_string_3.isnumeric():
                        test_map[int(check_string_3) - 1].append(line)
                        found_test_output = True
                        first_test_line = line_num

                # if output does not belong to a test, it is either from setup of from summary
                if found_test_output == False:
                    if line_num < first_test_line:
                        setup_lines.append(line)
                    elif line != "\n":
                        summary_lines.append(line)

    # write sorted file
    with open(log_file_sorted, "w") as file_out:
        for line in setup_lines:
            file_out.write(line)

        for test in test_map:
            for line in test:
                file_out.write(line)

        file_out.write("\n")

        for line in summary_lines:
            file_out.write(line)
