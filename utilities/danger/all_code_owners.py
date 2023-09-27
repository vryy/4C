from codeowners import CodeOwners
import argparse

"""
Support script for danger bot: determine all code owners for a given list of
files. All code owners are concatenated into a space-separated string.
"""


def all_code_owners(codeowner_file, file_list):
    with open(codeowner_file, "r") as cf:
        owners = CodeOwners(cf.read())

    all_owners = set()

    for f in file_list:
        for owner_type, owner_tag in owners.of(f):
            assert owner_type == "TEAM"
            all_owners.add(owner_tag)

    return all_owners


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("codeowner_file", type=str)
    parser.add_argument("--files", nargs="*", type=str)
    parser.add_argument("--exclude_owners", nargs="*", type=str)
    args = parser.parse_args()

    print(
        " ".join(
            all_code_owners(args.codeowner_file, args.files).difference(
                set(args.exclude_owners)
            )
        )
    )
