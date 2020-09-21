import glob
import hashlib
import json
import os
import time

import yaml
from joblib import Parallel, delayed

from diagnostic_item import DiagnosticItem


def read(yamlfiles, src):

    unique_items = set()

    items = []

    for yamlfile in yamlfiles:
        try:
            with open(yamlfile, "r") as fhandle:
                content = yaml.safe_load(fhandle)
                if not content:
                    continue  # Skip empty files.

                for item in content.get("Diagnostics", []):

                    diagnosticItem = DiagnosticItem.parse(item, src)

                    if not diagnosticItem.is_of_interest():
                        continue

                    hashstr = diagnosticItem.hash()

                    if hashstr in unique_items:
                        continue  # Skip duplicates

                    unique_items.add(hashstr)
                    items.append(diagnosticItem)
        except:
            print("error during processing {0}".format(yamlfile))

    print("processed {0} files".format(len(yamlfiles)))
    return items


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    if len(lst) <= n:
        yield lst
        return

    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def read_parallel(yamlfiles, src, chunksize=10, max_task=1):
    processed_list = Parallel(n_jobs=max_task)(
        delayed(read)(chunk, src) for chunk in chunks(yamlfiles, chunksize)
    )
    return merge_diagnostic_items(processed_list)


def merge_diagnostic_items(items):
    unique_list = []
    unique_items = set()

    for diagnosticList in items:
        for diagnosticItem in diagnosticList:
            hashstr = diagnosticItem.hash()

            if hashstr in unique_items:
                continue

            unique_items.add(hashstr)
            unique_list.append(diagnosticItem)

    return unique_list


def sort_issues_by_file(issues):
    # sort diagnostic items by file
    issues_dict = {}
    for diagnosticItem in issues:
        if diagnosticItem.get_relative_path() not in issues_dict:
            issues_dict[diagnosticItem.get_relative_path()] = []

        issues_dict[diagnosticItem.get_relative_path()].append(diagnosticItem)

    return issues_dict


def sort_file_issues_by_offset(issues):
    return sorted(issues, key=lambda x: x.offset)


def write_list_of_issues(fout, lst):
    issues_dict = sort_issues_by_file(lst)

    hash_set = set()
    fout.write("[\n")
    first = True
    # go through files
    for issues in issues_dict.values():
        # sort issues for their offset
        issues_sorted = sort_file_issues_by_offset(issues)

        # open source file
        write_file_issues(fout, issues[0].file, issues_sorted, hash_set, first)
    fout.write("\n]")


def write_file_issues(fout, filename, issues, hash_set, first):
    current_line_number = 0
    previous_line = None
    current_line = None
    sum_offset = 0
    last_successful = True
    with open(filename, "r") as source_file:

        for diagnosticItem in issues:

            if sum_offset <= diagnosticItem.offset:
                for line in source_file:
                    previous_line = current_line
                    current_line = line
                    current_line_number += 1
                    sum_offset += len(line.encode("utf-8"))

                    if sum_offset > diagnosticItem.offset:
                        break

            # found line and line content
            diagnosticItem.line_number = current_line_number
            diagnosticItem.line = current_line
            diagnosticItem.previous_line = previous_line

            if not first and last_successful:
                fout.write(",\n")
            try:
                diagnosticItem.write(fout, hash_set)
                last_successful = True
            except:
                last_successful = False

            first = False
