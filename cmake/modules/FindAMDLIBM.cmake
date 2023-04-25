# include(FindLibraryWithDebug)

if(AMDLIBM_INCLUDES AND AMDLIBM_LIBRARIES)
  set(AMDLIBM_FIND_QUIETLY TRUE)
endif(AMDLIBM_INCLUDES AND AMDLIBM_LIBRARIES)
find_path(
  AMDLIBM_INCLUDES
  NAMES amdlibm.h
  PATHS ${AMDLIBMDIR}/include ${AMDLIBM_DIR}/include ${INCLUDE_INSTALL_DIR}
  )

find_library(
  AMDLIBM_LIBRARIES
  NAMES amdlibm
  PATHS ${AMDLIBMDIR}/lib/static ${AMDLIBM_DIR}/lib/static ${LIB_INSTALL_DIR}
  )

find_file(
  AMDLIBM_LIBRARIES
  NAMES libamdlibm.so
  PATHS /usr/lib ${AMDLIBMDIR}/lib/dynamic ${LIB_INSTALL_DIR}
  )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(AMDLIBM DEFAULT_MSG AMDLIBM_INCLUDES AMDLIBM_LIBRARIES)
mark_as_advanced(AMDLIBM_INCLUDES AMDLIBM_LIBRARIES)

if(AMDLIBM_FOUND AND NOT TARGET amdlibm::amdlibm)
  add_library(amdlibm::amdlibm UNKNOWN IMPORTED)
  set_target_properties(
    amdlibm::amdlibm
    PROPERTIES IMPORTED_LOCATION "${AMDLIBM_LIBRARIES}"
               INTERFACE_INCLUDE_DIRECTORIES "${AMDLIBM_INCLUDES}"
    )
endif()

if(AMDLIBM_FOUND)
  list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS amdlibm::amdlibm)
  message(STATUS "Found AMD libM: ${AMDLIBM_LIBRARIES}")
  message(STATUS "Found AMD libM: ${AMDLIBM_INCLUDES}")
endif()
