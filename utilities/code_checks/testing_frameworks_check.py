from codecs import namereplace_errors
import dataclasses
import os, sys, argparse, re, json
import common_utils as utils
from typing import Optional, List


@dataclasses.dataclass(frozen=True)
class TestMacro:
  #  regular expression matching all test names and mpi ranks
  pattern: re.Pattern

  # category used during whitelisting of tests
  category: str

  # group ids in the regular expression that matches a test name
  test_names_groups: List[int] = dataclasses.field(default_factory=lambda: [0])
  
  # group id in the regular expression that matches the mpi rank
  mpirank_group: Optional[int] = None

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
    # the name of the class must be in the same line as ": public BACICxxTestWrapper"
    with open(fname, 'r') as file:
      # go through all lines
      testsuite_re = re.compile(r':\s*public\s*BACICxxTestWrapper')
      testsuite_test_re = re.compile(
        r'class *([a-zA-Z0-9\-_:]*) *: *public *BACICxxTestWrapper')
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

def check_mpirank_against_whitelist(name: str, category: str, mpi_rank: int, whitelist: List) -> bool:
  for item in whitelist:
    if item["test_name"] == name and item["test_category"] == category and item["mpi_rank"] == mpi_rank:
      return True
  return False

# CHECK INPUT FILE TESTS
def check_inputtests(look_cmd, allerrors):
  errors = 0
  input_tests = [str(ff) for ff in utils.files_changed(look_cmd)[:-1] if utils.is_input_file(ff)]

  # load mpi rank whitelist
  with open(os.path.join(sys.path[0],'whitelist_mpi_ranks.json')) as whitelist_file:
      mpirank_whitelist = json.load(whitelist_file)
  
  mpi_non_compliant_tests: List[str] = []
  list_of_all_testnames: List[str] = []

  # read TestingFrameworkListOfTests.cmake
  with open('TestingFrameworkListOfTests.cmake', 'r') as cmakefile:

    all_lines = "\n".join(cmakefile.readlines())

    # Defining the list of test macros with regular expressions that recognize them in the cmake file
    test_macros = [
      TestMacro(re.compile(r'baci_test\s*\(*([a-zA-Z0-9_\.\-]+)\s+(\d+)'), "baci_test", mpirank_group=1),
      TestMacro(re.compile(r'baci_test_extended_timeout\s*\(*([a-zA-Z0-9_\.\-]+)\s+(\d+)'), "baci_test_extended_timeout", mpirank_group=1),
      TestMacro(re.compile(r'baci_test_and_post_ensight_test\s*\(*([a-zA-Z0-9_\.\-]+)\s+(\d+)'), "baci_test_and_post_ensight_test", mpirank_group=1),
      TestMacro(re.compile(r'baci_test_restartonly\s*\(*([a-zA-Z0-9_\.\-]+)\s+(?:[a-zA-Z0-9_\.\-]+)\s+(\d+)'), "baci_test_restartonly", mpirank_group=1),
      TestMacro(re.compile(r'baci_test_Nested_Par\s*\(\s*([a-zA-Z0-9_\.\-]+)\s+([a-zA-Z0-9_\.\-]+)\s*'), "Nested_Par", test_names_groups=[0, 1]),
      TestMacro(re.compile(r'baci_test_Nested_Par_CopyDat\s*\(\s*([a-zA-Z0-9_\.\-]+)\s+(\d+)'), "Nested_Par_CopyDat", mpirank_group=1),
      TestMacro(re.compile(r'baci_test_Nested_Par_CopyDat_prepost\s*\(\s*([a-zA-Z0-9_\.\-]+)\s+([a-zA-Z0-9_\.\-]+)\s+([a-zA-Z0-9_\.\-\"\"]+)\s+(\d+)'), "Nested_Par_CopyDat_prepost", test_names_groups=[0, 1, 2], mpirank_group=3),
    ]


    for test_macro in test_macros:
      # search for tests of this macro
      for test in test_macro.pattern.findall(all_lines):
        # get list of test names defined in this macro (ignore empty test names)
        my_names = [test[i] for i in test_macro.test_names_groups if len(test[i].strip()) > 0 and test[i].strip() != '""']
        list_of_all_testnames.extend(my_names)

        # check mpi rank
        if test_macro.mpirank_group is not None:
          mpi_rank = int(test[test_macro.mpirank_group])

          if mpi_rank > 3:
            # check if tests are allowed to run with more than 4 procs
            for name in my_names:
              if not check_mpirank_against_whitelist(name, test_macro.category, mpi_rank, mpirank_whitelist):
                mpi_non_compliant_tests.append(name)

  if len(mpi_non_compliant_tests) > 0:
    errors += 1
    allerrors.append('The following tests use an unjustified high mpi rank:')
    allerrors.append('')
    allerrors.extend(mpi_non_compliant_tests)

  # check if some input tests are missing
  missing_input_tests = []
  for input_test in input_tests:
    # check, whether this input file is in TestingFrameworkListOfTests.cmake

    expected_test_name = os.path.splitext(os.path.basename(input_test))[0]
    if expected_test_name not in list_of_all_testnames:
      missing_input_tests.append(input_test)

  if len(missing_input_tests) > 0:
    errors += 1
    allerrors.append('The following input files are missing in TestingFrameworkListOfTests.cmake:')
    allerrors.append('')
    allerrors.extend(missing_input_tests)
  
  # check if input tests have empty sections
  tests_empty_sections = []

  for input_test in input_tests:
    with open(input_test, "r") as f:
      num_current_section_non_empty_lines = None
      
      for line in f:
        if line.startswith("--"):
          if num_current_section_non_empty_lines == 0:
            tests_empty_sections.append(input_test)
            break
          else:
            num_current_section_non_empty_lines = 0

        elif num_current_section_non_empty_lines is None:
          # No section title until now
          continue

        elif line.strip() != "":
          num_current_section_non_empty_lines += 1
  
  if len(tests_empty_sections) > 0:
    errors += 1
    allerrors.append('The following input files have empty sections. Please delete them or correct your input file.')
    allerrors.append('')
    allerrors.extend(tests_empty_sections)

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
      look_cmd = "git ls-files"

    # execute error checks

    # check unit tests
    errors += check_unittests(look_cmd, allerrors)
    if errors > 0:
      # add hint to unit test wiki page
      allerrors.extend(['Refer to our Wiki-page for good-style unittests:',
        'https://gitlab.lrz.de/baci/baci/wikis/Unit-testing:-good-practice-in-software-development'])

    # check input file tests
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
