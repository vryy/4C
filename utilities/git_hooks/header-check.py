#!/bin/env python2
" Baci header check to be executed along with the Git pre-commit hook. "
import os
import baciheader as bh
import inputheader as ih
#UTILS

def command_output(cmd):
  " Capture a command's standard output. "
  import subprocess
  message = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE).communicate()[0]
  return message

def is_source_file(fname):
  return os.path.splitext(fname)[1] in ".c .cpp .cxx .h .H .hpp".split()

def is_support_file(fname):
  return os.path.splitext(fname)[1] in ".dat .cmake .config .md".split()

def path_contains(test, path):
  head, tail = os.path.split(path)
  if test == tail or len(head) <= 0:
    return (test == tail)
  return path_contains(test, head)

def is_input_file(fname):
  return path_contains("Input", fname) and os.path.splitext(fname)[1] == ".dat"

def is_checked_file(fname):
  return is_source_file(fname) or is_support_file(fname)

def files_changed(look_cmd):
  """ List the files added or updated by this transaction.
  """
  def filename(line):
    return line
  looked_files  = command_output(look_cmd).decode().split("\n")
  changed_files = [filename(line) for line in looked_files]
  return changed_files

def file_contents(filename):
  " Return a file's contents for this transaction. "
  output=command_output("cat %s" %filename)
  return output

def pretty_print_error(allerrors):
  max_width = 56
  if len(allerrors) > 0:
    max_width = max(max_width,max([len(line) for line in allerrors]))
#header
  sys.stderr.write("\n"+"E"*(max_width+4)+"\nE"+" "*(max_width+2)+"E\n")
  sys.stderr.write("E Your commit was rejected due to the following reason(s):"+" "*(max_width-55)+"E\nE"+" "*(max_width+2)+"E\n")
#body
  for line in allerrors:
    sys.stderr.write("E "+line+" "*(max_width-len(line))+" E\n")
#footer
  sys.stderr.write("E"+" "*(max_width+2)+"E\n"+"E"*(max_width+4)+"\n")
  return

# CHECK FOR TABS

def contains_tabs(filename):
  " Return True if this version of the file contains tabs. "
  return "\t" in file_contents(filename)

def check_support_files_for_tabs(look_cmd, allerrors):
  " Check support files in this transaction are tab-free. "
  support_files_with_tabs = [ff for ff in files_changed(look_cmd) if is_support_file(ff) and contains_tabs(ff)]
  if len(support_files_with_tabs) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following support files contain tabs:")
    allerrors += support_files_with_tabs
  return len(support_files_with_tabs)

# CHECK FOR TRAILING WHITESPACES

def trailing_whitespace(filename):
  " Return True if this version of the file contains a trailing whitespace. "
  return " \n" in file_contents(filename)

def check_support_files_for_trailing_whitespace(look_cmd, allerrors):
  " Check support files in this transaction are trailing whitespace-free. "
  support_files_with_trailing_whitespace = [ff for ff in files_changed(look_cmd) if is_support_file(ff) and trailing_whitespace(ff)]
  if len(support_files_with_trailing_whitespace) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following support files contain trailing whitespaces:")
    allerrors += support_files_with_trailing_whitespace
  return len(support_files_with_trailing_whitespace)


#CHECK HEADER
def check_cpp_files_for_header(look_cmd, allerrors):
  " Check C/C++ files in this transaction have a maintainer. "
  import os
  headers = dict([(ff,bh.Header(file_contents(ff))) for ff in files_changed(look_cmd)[:-1] if is_source_file(ff)])
# \file tag
  cpp_files_wo_file = [ff for ff,hdr in headers.items() if hdr.get_file() != os.path.basename(ff)]
  if len(cpp_files_wo_file) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files have an incorrect or are missing a \\file tag:")
    allerrors += cpp_files_wo_file
# \brief tag
  cpp_files_wo_brief = []
#cpp_files_wo_brief = [ff for ff, hdr in headers.items() if len(hdr.get_brief()) < 5]
#if len(cpp_files_wo_brief) > 0:
#if len(allerrors) > 0:
#allerrors.append("")
#allerrors.append("The following files are missing a \\brief tag:")
#allerrors += cpp_files_wo_brief
# \maintainer tag
  cpp_files_wo_maint = [ff for ff,hdr in headers.items() if len(hdr.get_maintainer()) < 5]
  if len(cpp_files_wo_maint) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files are missing a \\maintainer tag:")
    allerrors += cpp_files_wo_maint

# check for correct start of header
  cpp_files_wrong_start = [ff for ff,hdr in headers.items() if len(hdr.get_start())>0]
  if len(cpp_files_wrong_start) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files do not start with /*! as an appropriate header marker:")
    allerrors += cpp_files_wrong_start
# \level tag
  cpp_files_wo_lvl = [ff for ff,hdr in headers.items() if not (0 <= hdr.get_level() <= 3)]
  if len(cpp_files_wo_lvl) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files are missing a \\level tag:")
    allerrors += cpp_files_wo_lvl
#print example header
  if len(cpp_files_wo_file) > 0 or len(cpp_files_wo_brief) > 0 or len(cpp_files_wrong_start) > 0 or len(cpp_files_wo_maint) > 0 or len(cpp_files_wo_lvl) > 0:
    allerrors += bh.Header.get_example()

  return len(cpp_files_wo_file)+len(cpp_files_wo_brief)+len(cpp_files_wo_maint)+len(cpp_files_wo_lvl)+len(cpp_files_wrong_start)

#CHECK INPUT FILE HEADERS
def check_input_files_for_header(look_cmd, allerrors):
  " Check .dat file in the Input folder for a proper header. "
  headers = dict([(ff,ih.Header(file_contents(ff))) for ff in files_changed(look_cmd)[:-1] if is_input_file(ff)])
  datfiles_without_header = [ff for ff,hdr in headers.items() if len(hdr.get_maintainer()) < 5]
  if len(datfiles_without_header) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files are missing a maintainer:")
    allerrors += datfiles_without_header
#print example header
  if len(datfiles_without_header) > 0:
    allerrors += ih.Header.get_example()
  return len(datfiles_without_header)

#CHECK FOR.GITIGNORE FILES

def build_filter(look_cmd):
  " Build a regex from entries found in .gitignore. "
  import pathspec
  #text= file_contents(".gitignore")
  with open('.gitignore', 'r') as text:
    spec = pathspec.PathSpec.from_lines(pathspec.GitIgnorePattern,text)
  return spec

def check_all_files_for_gitignore(look_cmd, allerrors):
  " Check that the files to be added or changed are not in .gitignore. "
  ignored_files = build_filter(look_cmd).match_files(files_changed(look_cmd))
  if len(ignored_files) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files are on the ignore list and may not be commited:")
    allerrors += ignored_files
  return len(ignored_files)


#######################################################################################################################

def main():
  errors = 0
  allerrors = []
  try:
    look_cmd = "git diff --name-only --cached --diff-filter=MRAC"
    errors += check_cpp_files_for_header(look_cmd, allerrors)
    errors += check_input_files_for_header(look_cmd, allerrors)
    errors += check_support_files_for_tabs(look_cmd, allerrors)
    errors += check_support_files_for_trailing_whitespace(look_cmd, allerrors)
    #errors += check_all_files_for_gitignore(look_cmd, allerrors)  # Did not work in latest Python env.
  except ValueError:
    print("Something went wrong! Check the error functions in this script again!")
    errors += 1
  if errors > 0:
    pretty_print_error(allerrors)
  return errors

if __name__ == "__main__":
  import sys
  sys.exit(main())
