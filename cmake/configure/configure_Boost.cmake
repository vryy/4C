# This is needed on kaiser to enforce the usage of the non-system boost libraries
option(BOOST_EXCLUDE_SYSTEM_PATHS "Avoid boost libraries in system paths" OFF)
if(BOOST_EXCLUDE_SYSTEM_PATHS)
  set(Boost_NO_SYSTEM_PATHS ON)
endif(BOOST_EXCLUDE_SYSTEM_PATHS)

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

  target_link_libraries(baci_all_enabled_external_dependencies INTERFACE Boost::system Boost::graph)
endif()
