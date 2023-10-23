import os
import yaml
import attrs
import cattrs
from typing import List, Optional, Dict, Set
from multiprocessing import Pool
import tqdm
import json
import argparse
import hashlib
import struct

# choose fastest yaml reader
try:
    from yaml import CLoader

    yamlreader = lambda file: yaml.load(file, Loader=CLoader)
except ImportError:
    yamlreader = yaml.safe_load


@attrs.define(frozen=True)
class ReplacementItem:
    FilePath: str
    Offset: int
    Length: int
    ReplacementText: str


@attrs.define(frozen=True)
class DiagnosticMessage:
    Message: str
    FilePath: str
    FileOffset: int
    Replacements: List[ReplacementItem]


@attrs.define(frozen=True)
class NotesItem:
    Message: str
    FilePath: str
    FileOffset: int
    Replacements: List[ReplacementItem]


@attrs.define(frozen=True)
class DiagnosticItem:
    DiagnosticName: str
    DiagnosticMessage: DiagnosticMessage
    Level: str
    BuildDirectory: str
    Notes: Optional[List[NotesItem]] = None


def get_diagnostic_item_hash(item: DiagnosticItem) -> bytes:
    sha = hashlib.sha256()
    sha.update(item.DiagnosticName.encode())
    sha.update(item.DiagnosticMessage.FilePath.encode())
    sha.update(struct.pack("<i", item.DiagnosticMessage.FileOffset))

    return sha.digest()


def read(file: str) -> List[DiagnosticItem]:
    with open(file, "r") as f:
        content = yamlreader(f)
        if content is None:
            return []
        return cattrs.structure(content["Diagnostics"], List[DiagnosticItem])


def load_clang_tidy_issues(yaml_files: List[str], num_procs=10) -> List[DiagnosticItem]:
    all_hashes = set()
    unique_items: List[DiagnosticItem] = []
    with Pool(num_procs) as p:
        for parsed_items in tqdm.tqdm(
            p.imap_unordered(read, yaml_files), total=len(yaml_files)
        ):
            for item in parsed_items:
                h = get_diagnostic_item_hash(item)
                if h not in all_hashes:
                    all_hashes.add(h)
                    unique_items.append(item)

    return unique_items


def get_file_issues_dict(
    issues: List[DiagnosticItem], src_dir: str
) -> Dict[str, List[DiagnosticItem]]:
    issues_dict = {}

    for item in issues:
        if not item.DiagnosticMessage.FilePath:
            # this error is not associated to a file
            continue
        relative_path = os.path.relpath(item.DiagnosticMessage.FilePath, src_dir)
        if relative_path not in issues_dict:
            issues_dict[relative_path] = []

        issues_dict[relative_path].append(item)

    return issues_dict


def sort_for_offset(issues: List[DiagnosticItem]) -> List[DiagnosticItem]:
    return sorted(issues, key=lambda x: x.DiagnosticMessage.FileOffset)


def get_fingerprint(
    name: str,
    relative_path: str,
    line: str,
    previous_line: str,
    used_fingerprints: Set[str],
) -> str:
    fingerprint = hashlib.sha256(
        (name + relative_path + line + previous_line).encode()
    ).hexdigest()

    # we compute the fingerprint based on the name of the file, the name of the check and the line
    # with the previous line. For some tests, it happens that this is not unique, e.g. a redundnant
    # "return;" statement is often following an empty line. In that case, distinct errors would have
    # the same fingerprint. To prevent this, we append a number if the fingerprint has already been
    # used.

    i = 0
    base_fingerprint = fingerprint
    while fingerprint in used_fingerprints:
        fingerprint = base_fingerprint + str(i)
        i += 1

    used_fingerprints.add(fingerprint)
    return fingerprint


def write_code_climate_report_for_file(
    outstream, file: str, issues: List[DiagnosticItem], source_directory: str
):
    sorted_issues = sort_for_offset(issues)

    sum_offset = 0
    previous_line = None
    current_line = None
    current_line_number = 0

    first_item = True

    used_fingerprints = set()

    with open(file, "r") as source_file:
        for item in sorted_issues:
            if sum_offset <= item.DiagnosticMessage.FileOffset:
                for line in source_file:
                    previous_line = current_line
                    current_line = line
                    sum_offset += len(line.encode("utf-8"))
                    current_line_number += 1

                    if sum_offset > item.DiagnosticMessage.FileOffset:
                        break

            relative_path = os.path.relpath(
                item.DiagnosticMessage.FilePath, source_directory
            )
            fingerprint = get_fingerprint(
                item.DiagnosticName,
                relative_path,
                current_line if current_line is not None else "",
                previous_line if previous_line is not None else "",
                used_fingerprints,
            )

            climate_item = {
                "description": item.DiagnosticMessage.Message,
                "check_name": item.DiagnosticName,
                "fingerprint": fingerprint,
                "severity": "blocker" if item.Level == "Error" else "info",
                "location": {
                    "path": relative_path,
                    "lines": {"begin": current_line_number},
                },
            }

            if first_item:
                first_item = False
            else:
                outstream.write(",\n")
            outstream.write(json.dumps(climate_item))


def write_code_climate_report(
    filename: str, issues: List[DiagnosticItem], source_directory: str
):
    with open(filename, "w") as f:
        f.write("[\n")
        first_item = True
        for relpath, file_issues in get_file_issues_dict(
            issues, source_directory
        ).items():
            if first_item:
                first_item = False
            else:
                f.write(",\n")

            write_code_climate_report_for_file(
                f,
                os.path.join(source_directory, relpath),
                file_issues,
                source_directory,
            )

        f.write("\n]")


def main():
    parser = argparse.ArgumentParser(
        prog="Collect all clang-tidy issues and export them as a codeclimate report for Gitlab.",
    )
    parser.add_argument("yaml_files", metavar="yaml-files", nargs="*")
    parser.add_argument("--code-climate-file")
    parser.add_argument("--source-folder")
    parser.add_argument("-j", type=int, default=10)

    args = parser.parse_args()

    unique_items = load_clang_tidy_issues(args.yaml_files, args.j)

    # write as codeclimate report
    if args.code_climate_file is not None:
        if args.source_folder is None:
            raise RuntimeError(
                "If you want to export a code climate file, you need to specify the source folder with --source_folder [...]"
            )
        write_code_climate_report(
            args.code_climate_file, unique_items, args.source_folder
        )


if __name__ == "__main__":
    main()
