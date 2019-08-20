#!/bin/env python2
" Baci header check"

import os
import argparse
import baciheader as bh
import inputheader as ih
import subprocess
#UTILS

def command_output(cmd):
  " Capture a command's standard output. "
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
  looked_files  = command_output(look_cmd).split("\n")
  changed_files = [filename(line) for line in looked_files]
  return changed_files

def file_contents(filename):
  " Return a file's contents for this transaction. "
  with open (filename) as myfile:
    output=myfile.readlines()
  return output

def pretty_print_error(allerrors):
  max_width = 56
  if len(allerrors) > 0:
    max_width = max(max_width,max([len(line) for line in allerrors]))
  if test_flag:
  #header
    sys.stderr.write("\n"+"E"*(max_width+4)+"\nE"+" "*(max_width+2)+"E\n")
    sys.stderr.write("E Your commit was rejected due to the following reason(s):"+" "*(max_width-55)+"E\nE"+" "*(max_width+2)+"E\n")
  #body
    for line in allerrors:
      sys.stderr.write("E "+line+" "*(max_width-len(line))+" E\n")
#footer
    sys.stderr.write("E"+" "*(max_width+2)+"E\n"+"E"*(max_width+4)+"\n")
  else:
    f = open('wrong_format.txt', 'w')
  #header
    f.write("\n"+"E"*(max_width+4)+"\nE"+" "*(max_width+2)+"E\n")
  #body
    for line in allerrors:
      f.write("E "+line+" "*(max_width-len(line))+" E\n")
  #footer
    f.write("E"+" "*(max_width+2)+"E\n"+"E"*(max_width+4)+"\n")
    f.close()
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
  return "\t" in file_contents(filename)

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
  headers = dict([(ff,bh.Header(file_contents(ff))) for ff in files_changed(look_cmd)[:-1] if is_source_file(ff)])
# \brief tag
  cpp_files_wo_brief = []
# \maintainer tag
  cpp_files_wo_maint = [ff for ff,hdr in headers.items() if len(hdr.get_maintainer()) < 1]
  if len(cpp_files_wo_maint) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files have an incorrect or are missing a \\maintainer tag:")
    allerrors += cpp_files_wo_maint
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
# check compliance of maintainer with the baci_developers.json (List of active BACI developers)
  cpp_files_uncompliant_maintainer = [ff for ff,hdr in headers.items() if len(hdr.get_compliant_maintainer())<1]
  if len(cpp_files_uncompliant_maintainer) > 0:
    if len(allerrors) > 0:
        allerrors.append("")
    allerrors.append("The following files contain an uncompliant maintainer. Please check './utilities/code_checks/baci_developers.json' for a list of compliant maintainers. Several maintainers are not allowed!")
    allerrors += cpp_files_uncompliant_maintainer

#print example header
  if len(cpp_files_wo_brief) > 0 or len(cpp_files_wrong_start) > 0 or len(cpp_files_wo_maint) > 0 or len(cpp_files_wo_lvl) > 0 or len(cpp_files_uncompliant_maintainer):
    allerrors += bh.Header.get_example()

  return len(cpp_files_wo_brief)+len(cpp_files_wo_maint)+len(cpp_files_wo_lvl)+len(cpp_files_wrong_start)+len(cpp_files_uncompliant_maintainer)


#CHECK INPUT FILE HEADERS
def check_input_files_for_header(look_cmd, allerrors):
  " Check .dat file in the Input folder for a proper header. "
  headers = dict([(ff,ih.Header(file_contents(ff))) for ff in files_changed(look_cmd)[:-1] if is_input_file(ff)])
  datfiles_without_header = [ff for ff,hdr in headers.items() if len(hdr.get_maintainer()) < 1]
  if len(datfiles_without_header) > 0:
    if len(allerrors) > 0:
      allerrors.append("")
    allerrors.append("The following files are missing a maintainer:")
    allerrors += datfiles_without_header

# check compliance of maintainer in input files with the baci_developers.json (List of active BACI developers)
  dat_files_uncompliant_maintainer = [ff for ff,hdr in headers.items() if len(hdr.get_compliant_maintainer())<1]

  if len(dat_files_uncompliant_maintainer) > 0:
    if len(allerrors) > 0:
        allerrors.append("")
    allerrors.append("The following input files (.dat) contain an uncompliant maintainer. Please check './utilities/code_checks/baci_developers.json' for a list of compliant maintainers. Several maintainers are not allowed!")
    allerrors += dat_files_uncompliant_maintainer

#print example header
  if len(datfiles_without_header) > 0 or len(dat_files_uncompliant_maintainer)>0:
    allerrors += ih.Header.get_example()

  return len(datfiles_without_header)+len(dat_files_uncompliant_maintainer)

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
  parser = argparse.ArgumentParser()
  parser.add_argument("--test_flag", help="Please specify the test_flag to determine the type of the code check. A value of 1 means pre-commit hook and a value of 0 a nightly test case.", type=int)
  args = parser.parse_args()
  flag = args.test_flag
  global test_flag
  if flag==1:
    test_flag = True # pre-commit hook
  elif flag==0:
    test_flag = False # nightly test case
  else: raise ValueError(r'Invalid argument for code checks! Must be integer')
  errors = 0
  allerrors = []
  try:
    if test_flag:
      look_cmd = "git diff --name-only --cached --diff-filter=MRAC"
    else:
      look_cmd = "find ./ -type f"
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
