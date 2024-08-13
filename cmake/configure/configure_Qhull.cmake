find_package(Qhull REQUIRED)

if(QHULL_FOUND)
  message(STATUS "QHULL include directory: ${QHULL_INCLUDE_DIR}")
  message(STATUS "QHULL library directory: ${QHULL_LIBRARY}")

  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE qhull::qhull)
endif()
