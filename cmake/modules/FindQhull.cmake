# Finder for Qhull
# Exports qhull::qhull as an imported target
# Note: The Qhull_ROOT variable is automatically considered by the find_ calls below.

find_path(QHULL_INCLUDE_DIR libqhull/libqhull.h)

find_library(QHULL_LIBRARY NAMES qhull)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qhull DEFAULT_MSG QHULL_LIBRARY QHULL_INCLUDE_DIR)

if(QHULL_FOUND AND NOT TARGET qhull::qhull)
  add_library(qhull::qhull UNKNOWN IMPORTED)
  set_target_properties(
    qhull::qhull
    PROPERTIES IMPORTED_LOCATION "${QHULL_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${QHULL_INCLUDE_DIR}"
    )
endif()

if(QHULL_FOUND)
  message(STATUS "QHULL include directory: ${QHULL_INCLUDE_DIR}")
  message(STATUS "QHULL library directory: ${QHULL_LIBRARY}")
  baci_add_dependency(baci_all_enabled_external_dependencies qhull::qhull)
endif(QHULL_FOUND)
