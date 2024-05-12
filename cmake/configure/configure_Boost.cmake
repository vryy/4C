find_package(
  Boost
  COMPONENTS graph system
  REQUIRED
  )

# post-process found targets
if(Boost_FOUND)
  message(STATUS "Boost component libraries: ${Boost_LIBRARIES}")
  message(STATUS "Boost libraries directory: ${Boost_LIBRARY_DIRS}")

  target_compile_definitions(
    Boost::graph
    INTERFACE "-DBOOST_MAJOR_VERSION=${Boost_MAJOR_VERSION}"
              "-DBOOST_MINOR_VERSION=${Boost_MINOR_VERSION}"
    )
  target_compile_definitions(
    Boost::system
    INTERFACE "-DBOOST_MAJOR_VERSION=${Boost_MAJOR_VERSION}"
              "-DBOOST_MINOR_VERSION=${Boost_MINOR_VERSION}"
    )

  target_link_libraries(
    four_c_all_enabled_external_dependencies INTERFACE Boost::system Boost::graph
    )
endif()
