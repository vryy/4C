import os, sys, argparse, re
import common_utils as utils

# CHECK UNITTESTS
def check_unittests(look_cmd, allerrors):
  errors = 0

  # get list of all unit tests
  unit_tests = [str(ff) for ff in utils.files_changed(look_cmd)[:-1] if utils.is_unittest_file(ff)]

  for fname in unit_tests:
    # check for naming convention
    if not str.startswith(os.path.split(fname)[1], 'unit_'):
      errors += 1
      allerrors.append('Unit test has to start with unit_')

    # check for naming convention of file extension
    if os.path.splitext(fname)[1] != ".H":
      errors += 1
      allerrors.append('The file extension of unit tests have to be .H')

    file_parts = utils.path_split_all(fname)
    # check, whether there is a corresponding source file to the unit test

    # NOTE: We decided not to enforce this for now, as there could be a reason for unit_tests that do not
    # belong to one specific file.
    #
    # 
    # basename = os.path.splitext(file_parts[-1])[0][
    #            5:]  # remove extension and unit_
    # basepath = os.path.sep.join(file_parts[1:-1])
    # file_exts = ['.h', '.H', '.hpp', '.cpp', '.c', '.cxx']
    # file_exists = False
    # for ext in file_exts:
    #   if os.path.isfile(os.path.join(basepath, '{0}{1}'.format(basename, ext))):
    #     file_exists = True
    #     break

    # if not file_exists:
    #   errors += 1
    #   allerrors.append(
    #     'There is no corresponding source file for the unittest {0}'.format(
    #       fname))

    # check, whether there is a corresponding CMakeLists.txt
    unittest_package = os.path.sep.join(file_parts[0:-1])
    if not os.path.isfile(os.path.join(unittest_package, 'CMakeLists.txt')):
      errors += 1
      allerrors.append(
        'There is no corresponding CMakeLists.txt file for the unittest package {0}'.format(
          unittest_package))
    else:
      # check, whether the file is added to the corresponding CMakeLists.txt file on the same level
      # IN ALPHABETICAL ORDER
      with open(os.path.join(unittest_package, 'CMakeLists.txt'),
                'r') as cmakefile:
        entries = []
        found = False
        entry_regex = re.compile(
          r'\$\{CMAKE_CURRENT_SOURCE_DIR\}\/([a-zA-Z0-9_\-\.]*)')
        for line in cmakefile:
          line_entries = entry_regex.findall(line)
          for entry in line_entries:
            entries.append(entry)
            if entry == file_parts[-1]:
              found = True

      if not found:
        errors += 1
        allerrors.append(
          'The unittest {0} is missing in the corresponding CMakeLists.txt file {1}'.format(
            fname, os.path.join(unittest_package, 'CMakeLists.txt')))

      if not utils.is_alphabetically(entries):
        errors += 1
        allerrors.append(
          'Tests are not in alphabetical order in CMakeLists.txt file {0}'.format(
            os.path.join(unittest_package, 'CMakeLists.txt')))

    # check the content of the unittest file itself
    # the name of the class must be in the same line as ": public CxxTest::TestSuite"
    with open(fname, 'r') as file:
      # go through all lines
      testsuite_re = re.compile(r':\s*public\s*CxxTest::TestSuite')
      testsuite_test_re = re.compile(
        r'class *([a-zA-Z0-9\-_:]*) *: *public *CxxTest::TestSuite')
      for line in file:
        if testsuite_re.search(line) is not None:
          # this line contains a TestSuite class.
          # There is no line break allowed here!
          if testsuite_test_re.search(line) is None:
            errors += 1
            allerrors.extend(
                '''The unittesttest {0} contains a class definition with a line break.
  
  The framework does not allow to have a line break between class NAMESPACE::Class_TestSuite : public Cxx::TestSuite.
  
  If the line is too long to be conform with the code style restrictions, add a "// clang-format off"
  before the class definition and "// clang-format on" after it. Example:
  
  ```
  // clang-format off
  class NAMESPACE::VeryLongClassName_TestSuite : public Cxx::TestSuite
  // clang-format on
  ```
  
  This will avoid automatic added linebreaks in between.'''.format(fname).split('\n'))

    
  return errors

# CHECK INPUT FILE TESTS
def check_inputtests(look_cmd, allerrors):
  errors = 0
  input_tests = [str(ff) for ff in utils.files_changed(look_cmd)[:-1] if utils.is_input_file(ff)]

  # read TestingFramework.cmake
  entries = []
  with open('TestingFramework.cmake', 'r') as cmakefile:
    entry_regex = [
      re.compile(r'baci_test *\( *([a-zA-Z0-9_\.\-]+) *'),
      re.compile(r'baci_test_restartonly *\( *([a-zA-Z0-9_\.\-]+) *'),
      re.compile(r'baci_test_Nested_Par *\( *([a-zA-Z0-9_\.\-]+) *'),
      re.compile(r'baci_test_Nested_Par_MultipleInvana *\( *([a-zA-Z0-9_\.\-]+) *([a-zA-Z0-9_\.\-]+)*'),
      re.compile(r'baci_test_Nested_Par_CopyDat *\( *([a-zA-Z0-9_\.\-]+) *')
    ]
    for line in cmakefile:
      # split comments
      line = line.split('#', 1)[0]

      line_entries = []
      for regex in entry_regex:
        for item in regex.findall(line):
          if isinstance(item, tuple):
            line_entries.extend(list(item))
          else:
            line_entries.append(item)
      
      for entry in line_entries:
        entries.append(entry)

  for input_test in input_tests:
    # check, whether this input file is in TestingFramework.cmake
    found = False

    for entry in entries:
      if entry == os.path.splitext(os.path.split(input_test)[1])[0]:
        found = True

    if not found:
      errors += 1
      allerrors.append(
        'The input file {0} is missing in TestingFramework.cmake'.format(input_test))
  
  return errors


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
    
    # execute error checks

    # check unit tests
    errors += check_unittests(look_cmd, allerrors)
    if errors > 0:
      # add hint to unit test wiki page
      allerrors.extend(['Refer to our Wiki-page for good-style unittests:',
        'https://gitlab.lrz.de/baci/baci/wikis/Unit-testing:-good-practice-in-software-development'])

    # check input fiile tests
    errors += check_inputtests(look_cmd, allerrors)
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