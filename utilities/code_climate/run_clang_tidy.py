#!/usr/bin/env python
#
# ===- run-clang-tidy.py - Parallel clang-tidy runner ---------*- python -*--===#
#
# Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
# See https://llvm.org/LICENSE.txt for license information.
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
# ===------------------------------------------------------------------------===#
# FIXME: Integrate with clang-tidy-diff.py

"""
Parallel clang-tidy runner
==========================

Runs clang-tidy over all files in a compilation database. Requires clang-tidy
and clang-apply-replacements in $PATH.

Example invocations.
- Run clang-tidy on all files in the current working directory with a default
  set of checks and show warnings in the cpp files and all project headers.
    run-clang-tidy.py $PWD

Compilation database setup:
http://clang.llvm.org/docs/HowToSetupToolingForLLVM.html
"""

from __future__ import print_function

import argparse
import glob
import multiprocessing
import os
import queue as queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import json
import traceback

import yaml

import fast_yaml_reader


def find_compilation_database(path):
    """Adjusts the directory until a compilation database is found."""
    result = "./"
    while not os.path.isfile(os.path.join(result, path)):
        if os.path.realpath(result) == "/":
            print("Error: could not find compilation database.")
            sys.exit(1)
        result += "../"
    return os.path.realpath(result)


def make_absolute(f, directory):
    if os.path.isabs(f):
        return f
    return os.path.normpath(os.path.join(directory, f))


def get_tidy_invocation(
    f,
    clang_tidy_binary,
    tmpdir,
    build_path,
    header_filter,
    extra_arg,
    extra_arg_before,
    quiet,
):
    """Gets a command line for clang-tidy."""
    start = [clang_tidy_binary]
    if header_filter is not None:
        start.append("-header-filter=" + header_filter)
    if tmpdir is not None:
        start.append("-export-fixes")
        # Get a temporary file. We immediately close the handle so clang-tidy can
        # overwrite it.
        (handle, name) = tempfile.mkstemp(suffix=".yaml", dir=tmpdir)
        os.close(handle)
        start.append(name)
    for arg in extra_arg:
        start.append("-extra-arg=%s" % arg)
    for arg in extra_arg_before:
        start.append("-extra-arg-before=%s" % arg)
    start.append("-p=" + build_path)
    if quiet:
        start.append("-quiet")
    start.append(f)
    return start


def run_tidy(args, tmpdir, build_path, queue, lock, failed_files):
    """Takes filenames out of queue and runs clang-tidy on them."""
    while True:
        name = queue.get()
        invocation = get_tidy_invocation(
            name,
            args.clang_tidy_binary,
            tmpdir,
            build_path,
            args.header_filter,
            args.extra_arg,
            args.extra_arg_before,
            args.quiet,
        )
        with lock:
            sys.stdout.write("Analyzing file {0}\n".format(name))
            sys.stdout.flush()

        process = subprocess.Popen(
            " ".join(invocation),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )
        _, stderr = process.communicate()
        with lock:
            if process.returncode != 0:
                failed_files.append(name)
                sys.stderr.buffer.write(stderr)
        queue.task_done()


def main():
    parser = argparse.ArgumentParser(
        description="Runs clang-tidy over all files "
        "in a compilation database. Requires "
        "clang-tidy in "
        "$PATH."
    )
    parser.add_argument(
        "-clang-tidy-binary",
        metavar="PATH",
        default="clang-tidy",
        help="path to clang-tidy binary",
    )
    parser.add_argument(
        "-header-filter",
        default=None,
        help="regular expression matching the names of the "
        "headers to output diagnostics from. Diagnostics from "
        "the main file of each translation unit are always "
        "displayed.",
    )
    if yaml:
        parser.add_argument(
            "-export",
            metavar="filename",
            dest="export",
            help="Create a codeclimate.json file to store issues in.",
        )
    parser.add_argument(
        "-j",
        type=int,
        default=0,
        help="number of tidy instances to be run in parallel.",
    )
    parser.add_argument(
        "files", nargs="*", default=[".*"], help="files to be processed (regex on path)"
    )
    parser.add_argument(
        "-p", dest="build_path", help="Path used to read a compile command database."
    )
    parser.add_argument(
        "-extra-arg",
        dest="extra_arg",
        action="append",
        default=[],
        help="Additional argument to append to the compiler " "command line.",
    )
    parser.add_argument(
        "-extra-arg-before",
        dest="extra_arg_before",
        action="append",
        default=[],
        help="Additional argument to prepend to the compiler " "command line.",
    )
    parser.add_argument(
        "-quiet", action="store_true", help="Run clang-tidy in quiet mode"
    )
    parser.add_argument("-src", help="Path to the source directory")
    args = parser.parse_args()

    db_path = "compile_commands.json"

    if args.build_path is not None:
        build_path = args.build_path
    else:
        # Find our database
        build_path = find_compilation_database(db_path)

    # Load the database and extract all files.
    database = json.load(open(os.path.join(build_path, db_path)))
    files = [make_absolute(entry["file"], entry["directory"]) for entry in database]

    try:
        invocation = [args.clang_tidy_binary, "-list-checks"]
        invocation.append("-p=" + build_path)

        invocation.append(files[0])
        print(invocation)
        if args.quiet:
            # Even with -quiet we still want to check if we can call clang-tidy.
            with open(os.devnull, "w") as dev_null:
                subprocess.check_call(" ".join(invocation), stdout=dev_null, shell=True)
        else:
            subprocess.check_call(" ".join(invocation), shell=True)
    except:
        print("Unable to run clang-tidy.", file=sys.stderr)
        sys.exit(1)

    max_task = args.j
    if max_task == 0:
        max_task = multiprocessing.cpu_count()

    tmpdir = tempfile.mkdtemp()

    # Build up a big regexy filter from all command line arguments.
    file_name_re = re.compile("|".join(args.files))

    return_code = 0
    try:
        # Spin up a bunch of tidy-launching threads.
        task_queue = queue.Queue(max_task)
        # List of files with a non-zero return code.
        failed_files = []
        lock = threading.Lock()
        for _ in range(max_task):
            t = threading.Thread(
                target=run_tidy,
                args=(args, tmpdir, build_path, task_queue, lock, failed_files),
            )
            t.daemon = True
            t.start()

        count_files = 0
        for name in files:
            if file_name_re.search(name):
                count_files += 1

        print("Analyzing {0} files".format(count_files))

        # Fill the queue with files.
        for name in sorted(files):
            if file_name_re.search(name):
                task_queue.put(name)

        # Wait for all threads to be done.
        task_queue.join()
        if len(failed_files):
            return_code = 1

    except KeyboardInterrupt:
        # This is a sad hack. Unfortunately subprocess goes
        # bonkers with ctrl-c and we start forking merrily.
        print("\nCtrl-C detected, goodbye.")
        if tmpdir:
            shutil.rmtree(tmpdir)
        os.kill(0, 9)

    if args.export:
        try:
            print("Parsing issues...")
            lst = list(glob.iglob(os.path.join(tmpdir, "*.yaml")))
            issues = fast_yaml_reader.read_parallel(
                lst, args.src, chunksize=10, max_task=max_task
            )

            print("Writing issues to {0}...".format(args.export))
            with open(args.export, "w") as fout:
                fast_yaml_reader.write_list_of_issues(fout, issues)
        except:
            print("Error exporting codeclimate report.\n", file=sys.stderr)
            traceback.print_exc()
            return_code = 1

    if tmpdir:
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    main()
