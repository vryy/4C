# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import re

"""
This script analyzes the log files of a test suite run and determines contiguous chunks of tests such that all chunks
have roughly the same run time.
"""


def find_run_time(log_file):
    """
    Regex needs to match lines like:
     Test #585: elch_3D_tet4_s2i_butlervolmer_mortar_standard-p3 ........................................................   Passed    3.65 sec
    """
    # Find all matches of the regex. Save in a dict with the test number as the key and the run time as the value.
    regex = re.compile(r"Test\s+#(\d+).*Passed\s+(\d+\.\d+)\s+sec")
    run_times = {}
    with open(log_file, "r") as f:
        for line in f:
            match = regex.search(line)
            if match:
                test_number = int(match.group(1))
                run_time = float(match.group(2))
                run_times[test_number] = run_time

    return run_times


def main():
    parser = argparse.ArgumentParser(description="Chunk Test Suite")
    parser.add_argument(
        "--chunks",
        "-c",
        type=int,
        help="Number of chunks to divide the test suite into",
    )
    parser.add_argument("log_files", type=str, nargs="+", help="Log file(s) to analyze")

    args = parser.parse_args()
    chunks = args.chunks
    log_files = args.log_files

    all_run_times = {}
    for log_file in log_files:
        all_run_times.update(find_run_time(log_file))

    # sort the dictionary by keys
    all_run_times = dict(sorted(all_run_times.items()))
    if list(all_run_times.keys()) != list(range(1, len(all_run_times) + 1)):
        raise ValueError("The test numbers in the log files are not contiguous.")

    average_run_time = sum(all_run_times.values()) / chunks

    chunks = []
    current_chunk_run_time = 0
    chunk_start = 1
    for test_number, run_time in all_run_times.items():
        current_chunk_run_time += run_time
        if current_chunk_run_time >= average_run_time:
            chunks.append((chunk_start, test_number, current_chunk_run_time))
            chunk_start = test_number + 1
            current_chunk_run_time = 0

    chunks.append((chunk_start, len(all_run_times), current_chunk_run_time))

    for i, chunk in enumerate(chunks):
        print(
            f"Chunk {i + 1}: {chunk[0]}-{chunk[1]} (expected run time: {chunk[2]:.2f} sec)"
        )

    print()

    print("Settings for github workflow:")
    print("test-chunk: [", end=" ")
    for chunk in chunks[:-1]:
        print(f'"{chunk[0]},{chunk[1]}"', end=", ")
    print(f'"{chunks[-1][0]},"]')


if __name__ == "__main__":
    main()
