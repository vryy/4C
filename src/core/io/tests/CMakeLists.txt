# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_core_io)

set(SOURCE_LIST
    4C_io_file_reader_test.cpp
    4C_io_inputreader_test.cpp
    4C_io_linedefinition_test.cpp
    4C_io_value_parser_test.cpp
    4C_io_pstream_test.cpp
    4C_io_string_converter_test.cpp
    )

file(GLOB_RECURSE SUPPORT_FILES CONFIGURE_DEPENDS test_files/*)

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  1
  SOURCE
  ${SOURCE_LIST}
  SUPPORT_FILES
  ${SUPPORT_FILES}
  )
