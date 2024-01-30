set(TESTNAME unittests_stl_extension)

set(SOURCE_LIST baci_utils_std_cxx20_ranges_test.cpp)

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  1
  SOURCE
  ${SOURCE_LIST}
  )

baci_link_google_test_necessary_libraries(${TESTNAME} core)
