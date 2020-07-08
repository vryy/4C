#!/bin/env python2
" Baci header check"

import os
import argparse
import baciheader as bh
import subprocess
import common_utils as utils

# CHECK FOR TABS

def contains_tabs(filename):
  " Return True if this version of the file contains tabs. "
  return "\t" in utils.file_contents(filename)

def check_support_files_for_tabs(look_cmd, allerrors):
  " Check support files in this transaction are tab-free. "

  support_files_with_tabs = [ff for ff in utils.files_changed(look_cmd) if utils.is_support_file(ff) and contains_tabs(ff)]
  if len(support_files_with_tabs) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following support files contain tabs:")
    allerrors += support_files_with_tabs
  return len(support_files_with_tabs)

# CHECK FOR TRAILING WHITESPACES

def trailing_whitespace(filename):
  " Return True if this version of the file contains a trailing whitespace. "
  return "\t" in utils.file_contents(filename)

def check_support_files_for_trailing_whitespace(look_cmd, allerrors):
  " Check support files in this transaction are trailing whitespace-free. "
  support_files_with_trailing_whitespace = [ff for ff in utils.files_changed(look_cmd) if utils.is_support_file(ff) and trailing_whitespace(ff)]
  if len(support_files_with_trailing_whitespace) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following support files contain trailing whitespaces:")
    allerrors += support_files_with_trailing_whitespace
  return len(support_files_with_trailing_whitespace)



#CHECK HEADER
def check_cpp_files_for_header(look_cmd, allerrors):
  " Check C/C++ files in this transaction. "
  headers = dict([(ff,bh.Header(utils.file_contents(ff))) for ff in utils.files_changed(look_cmd)[:-1] if utils.is_source_file(ff)])
# \brief tag
  cpp_files_wo_brief = []
# check for correct start of header
  cpp_files_wrong_start = [ff for ff,hdr in headers.items() if len(hdr.get_start())<1]
  if len(cpp_files_wrong_start) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files do not start with '/*! \\file' or '/** \\file' as an appropriate header marker:")
    allerrors += cpp_files_wrong_start
# \level tag
  cpp_files_wo_lvl = [ff for ff,hdr in headers.items() if not (0 <= hdr.get_level() <= 3)]
  if len(cpp_files_wo_lvl) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files are missing a \\level tag:")
    allerrors += cpp_files_wo_lvl

#print example header
  if len(cpp_files_wo_brief) > 0 or len(cpp_files_wrong_start) > 0 or len(cpp_files_wo_lvl) > 0:
    allerrors += bh.Header.get_example()

  return len(cpp_files_wo_brief)+len(cpp_files_wo_lvl)+len(cpp_files_wrong_start)

#CHECK FOR.GITIGNORE FILES

def build_filter(look_cmd):
  " Build a regex from entries found in .gitignore. "
  import pathspec
  #text= utils.file_contents(".gitignore")
  with open('.gitignore', 'r') as text:
    spec = pathspec.PathSpec.from_lines(pathspec.GitIgnorePattern,text)
  return spec

def check_all_files_for_gitignore(look_cmd, allerrors):
  " Check that the files to be added or changed are not in .gitignore. "
  ignored_files = build_filter(look_cmd).match_files(utils.files_changed(look_cmd))
  if len(ignored_files) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files are on the ignore list and may not be commited:")
    allerrors += ignored_files
  return len(ignored_files)


#######################################################################################################################

def main():
  # build command line arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('--diff_only', action='store_true', help='Add this tag if only the difference to HEAD should be analyzed. This flag should be used as a pre-commit hook. Otherwise all files are checked.')
  parser.add_argument('--out', type=str, default=None, help='Add this tag if the error message should be written to a file.')
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
      look_cmd = "find ./ -type f"
    errors += check_cpp_files_for_header(look_cmd, allerrors)
    errors += check_support_files_for_tabs(look_cmd, allerrors)
    errors += check_support_files_for_trailing_whitespace(look_cmd, allerrors)
    #errors += check_all_files_for_gitignore(look_cmd, allerrors)  # Did not work in latest Python env.
  except ValueError:
    print("Something went wrong! Check the error functions in this script again!")
    errors += 1
  if errors > 0:
    if errfile is None:
      utils.pretty_print_error_stderr(allerrors)
    else:
      utils.pretty_print_error_file(allerrors, errfile)
  return errors

if __name__ == "__main__":
  import sys
  sys.exit(main())
