# Define a special version of drt_error which will be linked to unit test executables
add_library(drt_error_for_testing dserror_for_testing.cpp)

set(LIBS_DRT_ERROR_FOR_TESTING
    # cmake-format: sortable
    drt_lib_error_interface
    )

target_link_libraries(drt_error_for_testing PUBLIC ${LIBS_DRT_ERROR_FOR_TESTING})
