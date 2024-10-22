# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from four_c_utils import patched_code_owners
import argparse
import sys

"""
Support script for danger bot: determine all code owners for a given list of
files. All code owners are concatenated into a space-separated string.
"""


def all_code_owners(codeowner_file, file_list):
    with open(codeowner_file, "r") as cf:
        owners = patched_code_owners.patched_code_owners(cf.read())

    all_owners = set()

    for f in file_list:
        for owner_type, owner_tag in owners.of(f):
            assert owner_type == "TEAM"
            all_owners.add(owner_tag)

    return all_owners


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("codeowner_file", type=str)
    parser.add_argument("--exclude_owners", nargs="*", type=str)
    args = parser.parse_args()

    # read files from stdin
    files = []
    for line in sys.stdin:
        files.append(line.strip())

    exclude_owners = []
    if args.exclude_owners is not None:
        exclude_owners.extend(args.exclude_owners)

    print(
        " ".join(
            all_code_owners(args.codeowner_file, files).difference(set(exclude_owners))
        )
    )
