# temporarily disable our own finders for module lookup
list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# This is needed on kaiser to enforce the usage of the non-system boost libraries
option(USE_NO_SYSTEM_PATHS "Avoid boost libraries in system paths" OFF)
if(USE_NO_SYSTEM_PATHS)
  set(Boost_NO_SYSTEM_PATHS ON)
endif(USE_NO_SYSTEM_PATHS)

# call the built-in finder with the specifications set outside
find_package(Boost)

# post-process found targets
if(Boost_FOUND)
  message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS}")
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

  list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS Boost::system Boost::graph)
endif()

# re-enable our own finders
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
