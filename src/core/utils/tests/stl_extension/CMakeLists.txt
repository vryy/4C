set(TESTNAME unittests_stl_extension)

set(SOURCE_LIST 4C_utils_lazy_ptr_test.cpp 4C_utils_std_cxx20_ranges_test.cpp)

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  1
  SOURCE
  ${SOURCE_LIST}
  )
