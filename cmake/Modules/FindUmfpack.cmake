if (Umfpack_INCLUDE_DIR AND Umfpack_LIBRARY)
  set(Umfpack_FIND_QUIETLY TRUE)
endif (Umfpack_INCLUDE_DIR AND Umfpack_LIBRARY)

if(CMAKE_Fortran_COMPILER_WORKS)

  find_package(BLAS)

  if(BLAS_FOUND)

    find_path(Umfpack_INCLUDE_DIR
      NAMES
      umfpack.h
      PATHS
      $ENV{UMFPACKDIR}
      ${INCLUDE_INSTALL_DIR}
      PATH_SUFFIXES
      suitesparse
      )

    find_library(Umfpack_LIBRARY umfpack PATHS $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})

    if(Umfpack_LIBRARY)

      get_filename_component(Umfpack_LIBDIR ${Umfpack_LIBRARY} PATH)

      find_library(AMD_LIBRARY amd PATHS ${Umfpack_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
      if (AMD_LIBRARY)
        set(Umfpack_LIBRARY ${Umfpack_LIBRARY} ${AMD_LIBRARY})
        #else (AMD_LIBRARY)
        #  set(Umfpack_LIBRARY FALSE)
      endif (AMD_LIBRARY)

    endif(Umfpack_LIBRARY)

#    if(Umfpack_LIBRARY)

#      find_library(COLAMD_LIBRARY colamd PATHS ${Umfpack_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
#      if (COLAMD_LIBRARY)
#        set(Umfpack_LIBRARY ${Umfpack_LIBRARY} ${COLAMD_LIBRARY})
#        #else (COLAMD_LIBRARY)
#        #  set(Umfpack_LIBRARY FALSE)
#      endif (COLAMD_LIBRARY)

#    endif(Umfpack_LIBRARY)

  endif(BLAS_FOUND)

endif(CMAKE_Fortran_COMPILER_WORKS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Umfpack DEFAULT_MSG Umfpack_INCLUDE_DIR Umfpack_LIBRARY)

IF(Umfpack_LIBRARY)
  set( Umfpack_FOUND "YES" )
  SET( Umfpack_LIBRARIES ${Umfpack_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY} ${BLAS_LIBRARIES})
ENDIF(Umfpack_LIBRARY)

mark_as_advanced(Umfpack_INCLUDE_DIR Umfpack_LIBRARY AMD_LIBRARY)
