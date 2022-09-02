set(TESTNAME unittests_drt_lib)

set(SOURCE_LIST
    unit_cubic_spline_interpolation.cpp
    unit_demangle.cpp
    unit_domainreader.cpp
    unit_filesystem_utils.cpp
    unit_function_library.cpp
    unit_functionvariables.cpp
    unit_linedefinition.cpp
    unit_tensor_transformation.cpp
    unit_utils_reference_configuration.cpp
    unit_voigt_notation.cpp
    )

baci_add_google_test_executable(${TESTNAME} ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME}
  -Wl,--start-group
  ${BACI_LIBRARIES}
  -Wl,--end-group
  ${LIBRARIES}
  )
