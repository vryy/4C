set(TESTNAME unittests_lib)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_cubic_spline_interpolation.cpp
    unit_demangle.cpp
    unit_discretization_nodal_coordinates.cpp
    unit_domainreader.cpp
    unit_function_library.cpp
    unit_functionvariables.cpp
    unit_gridgenerator.cpp
    unit_linedefinition.cpp
    unit_symbolic_expression.cpp
    unit_tensor_transformation.cpp
    unit_utils_reference_configuration.cpp
    unit_voigt_notation.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_link_executable_with_all(${TESTNAME})

set(TESTNAME unittests_lib_np_3)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_discretization_nodal_coordinates_np_3.cpp
    unit_gid_vector.cpp
    unit_gridgenerator_np_3.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_link_executable_with_all(${TESTNAME})
