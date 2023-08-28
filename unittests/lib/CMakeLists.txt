set(TESTNAME unittests_lib)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_cubic_spline_interpolation_test.cpp
    baci_demangle_test.cpp
    baci_discretization_nodal_coordinates_test.cpp
    baci_domainreader_test.cpp
    baci_function_library_test.cpp
    baci_functionvariables_test.cpp
    baci_gridgenerator_test.cpp
    baci_inputreader_test.cpp
    baci_linedefinition_test.cpp
    baci_symbolic_expression_test.cpp
    baci_tensor_transformation_test.cpp
    baci_utils_reference_configuration_test.cpp
    baci_voigt_notation_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_target_requires_cyclic_legacy_linking(${TESTNAME})

set(TESTNAME unittests_lib_np_3)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_discretization_nodal_coordinates_np_3_test.cpp
    baci_gid_vector_test.cpp
    baci_gridgenerator_np_3_test.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_target_requires_cyclic_legacy_linking(${TESTNAME})
