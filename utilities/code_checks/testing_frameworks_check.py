import os, sys, argparse, re
import common_utils as utils

# CHECK UNITTESTS
def check_unittests(look_cmd, allerrors):
  errors = 0

  # get list of all unit tests
  unit_tests = [str(ff) for ff in utils.files_changed(look_cmd)[:-1] if utils.is_unittest_file(ff)]
  unit_tests_cmake = [str(ff) for ff in utils.files_changed(look_cmd)[:-1] if utils.is_unittest_cmakefile(ff)]

  naming_convention = set()
  cmakelists = set()
  alphabetical_order = set()
  content = set()

  # check, whether there is a corresponding CMakeLists.txt
  cmake_regex_entries = [
    (re.compile(r'\$\{CMAKE_CURRENT_SOURCE_DIR\}\/([a-zA-Z0-9_\-\.]*)'), True),
    (re.compile(r'add_subdirectory\s*\(\s*([a-zA-Z0-9_\-\.]*)\s*\)'), False)
  ]

  for fname in unit_tests:
    # check for naming convention
    if not str.startswith(os.path.split(fname)[1], 'unit_'):
      naming_convention.add(fname)

    # check for naming convention of file extension
    if os.path.splitext(fname)[1] != ".H":
      naming_convention.add(fname)

    file_parts = utils.path_split_all(fname)
    unittest_package = os.path.sep.join(file_parts[0:-1])
    if not os.path.isfile(os.path.join(unittest_package, 'CMakeLists.txt')):
      cmakelists.add(fname)
    else:
      # check, whether the file is added to the corresponding CMakeLists.txt file on the same level
      # IN ALPHABETICAL ORDER
      with open(os.path.join(unittest_package, 'CMakeLists.txt'),
                'r') as cmakefile:
        found = False

        for line in cmakefile:

          for entry_regex in cmake_regex_entries:
            # this regex does not describe a unit test, it's a directory. Skip it
            if entry_regex[1] == False:
              continue
            line_entries = entry_regex[0].findall(line)
            for entry in line_entries:
              if entry == file_parts[-1]:
                found = True

      if not found:
        cmakelists.add(fname)

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
            content.add(fname)

  for fname in unit_tests_cmake:
    # check whether all tests are added in alphabetical order
    with open(fname, 'r') as cmakefile:
      entries = []
      for _ in range(len(cmake_regex_entries)):
        entries.append([])

      for line in cmakefile:
        for i, entry_regex in enumerate(cmake_regex_entries):
          line_entries = entry_regex[0].findall(line)
          entries[i].extend(line_entries)

      for l in entries:
        if not utils.is_alphabetically(l):
          alphabetical_order.add(fname)


  if len(naming_convention) > 0:
    errors += 1
    allerrors.append('The following unit tests are not compliant with the naming conventions.')
    allerrors.append('Unit tests have to start with unit_ and must have a .H file extension.')
    allerrors.append('')
    allerrors.extend(list(naming_convention))
    allerrors.append('')

  if len(cmakelists) > 0:
    errors += 1
    allerrors.append('The following unit tests are not added to a local CMakeLists.txt:')
    allerrors.append('')
    allerrors.extend(list(cmakelists))
    allerrors.append('')

  if len(alphabetical_order) > 0:
    errors += 1
    allerrors.append('The items in the following local CMakeLists.txt are not in alphabetical order:')
    allerrors.append('')
    allerrors.extend(list(alphabetical_order))
    allerrors.append('')

  if len(content) > 0:
    errors += 1
    allerrors.append('The contents of the following unit tests are not correct.')
    allerrors.append('')
    allerrors.append('The framework does not allow to have a line break between class NAMESPACE::Class_TestSuite : public Cxx::TestSuite.')
    allerrors.append('If the line is too long to be conform with the code style restrictions, add a "// clang-format off"')
    allerrors.append('before the class definition and "// clang-format on" after it. Example:')
    allerrors.append('')
    allerrors.append('```')
    allerrors.append('// clang-format off')
    allerrors.append('class NAMESPACE::VeryLongClassName_TestSuite : public Cxx::TestSuite')
    allerrors.append('// clang-format on')
    allerrors.append('')
    allerrors.append('This will avoid automatic added linebreaks in between.')
    allerrors.append('')
    allerrors.extend(list(content))
    allerrors.append('')

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
      re.compile(r'baci_test_Nested_Par *\( *([a-zA-Z0-9_\.\-]+) +([a-zA-Z0-9_\.\-]+) *'),
      re.compile(r'baci_test_Nested_Par_MultipleInvana *\( *([a-zA-Z0-9_\.\-]+) +([a-zA-Z0-9_\.\-]+)*'),
      re.compile(r'baci_test_Nested_Par_CopyDat *\( *([a-zA-Z0-9_\.\-]+) *'),
      re.compile(r'baci_test_Nested_Par_CopyDat_prepost *\( *([a-zA-Z0-9_\.\-]+)\s+([a-zA-Z0-9_\.\-]+)\s+([a-zA-Z0-9_\.\-]+) *'),
      re.compile(r'baci_test_Nested_Par_CopyDat_prepost *\( *([a-zA-Z0-9_\.\-]+) *')
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

  missing_input_tests = []
  for input_test in input_tests:
    # check, whether this input file is in TestingFramework.cmake
    found = False

    for entry in entries:
      if entry == os.path.splitext(os.path.split(input_test)[1])[0]:
        found = True

    if not found:
      missing_input_tests.append(input_test)

  if len(missing_input_tests) > 0:
    errors += 1
    allerrors.append('The following input files are missing in TestingFramework.cmake:')
    allerrors.append('')
    allerrors.extend(missing_input_tests)

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
  if len(allerrors) > 0:
    if errfile is None:
      utils.pretty_print_error_stderr(allerrors)
    else:
      utils.pretty_print_error_file(allerrors, errfile)
  return errors

if __name__ == "__main__":
  import sys
  sys.exit(main())