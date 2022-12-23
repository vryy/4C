# Define a special version of error which will be linked to unit test executables
add_library(error_for_testing dserror_for_testing.cpp)

set(LIBS_ERROR_FOR_TESTING
    # cmake-format: sortable
    lib_error_interface
    )

target_link_libraries(error_for_testing PUBLIC ${LIBS_ERROR_FOR_TESTING})

add_library(unittests_common INTERFACE)
target_include_directories(unittests_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unittests_common INTERFACE linalg)
