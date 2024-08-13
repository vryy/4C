find_package(CLN REQUIRED)

if(CLN_FOUND)
  message(STATUS "CLN include directory: ${CLN_INCLUDE_DIR}")
  message(STATUS "CLN library directory: ${CLN_LIBRARY}")

  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE cln::cln)
endif()
