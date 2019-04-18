#
# Find the PARMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of
# sparse matrices. It can be found at:
#       http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# PARMETIS_INCLUDE_DIR - where to find autopack.h
# PARMETIS_LIBRARIES   - List of fully qualified libraries to link against.
# PARMETIS_FOUND       - Do not attempt to use if "no" or undefined.

IF (PARMETIS_INCLUDE_DIR)
  # Already in cache, be silent
  SET(PARMETIS_FIND_QUIETLY TRUE)
ENDIF (PARMETIS_INCLUDE_DIR)

FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h
  ${INCLUDE_INSTALL_DIR}
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(PARMETIS_LIBRARY parmetis
  ${LIB_INSTALL_DIR}
  /usr/local/lib
  /usr/lib
)

FIND_LIBRARY(METIS_LIBRARY metis
  ${LIB_INSTALL_DIR}
  /usr/local/lib
  /usr/lib
)

#IF(PARMETIS_INCLUDE_DIR)
#  IF(PARMETIS_LIBRARY)
#
#    SET( PARMETIS_FOUND "YES" )
#  ENDIF(PARMETIS_LIBRARY)
#ENDIF(PARMETIS_INCLUDE_DIR)

# Per-recommendation
SET( PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR})
SET( PARMETIS_LIBRARIES    ${PARMETIS_LIBRARY} ${METIS_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PARMETIS DEFAULT_MSG PARMETIS_LIBRARIES PARMETIS_INCLUDE_DIRS)

MARK_AS_ADVANCED(PARMETIS_INCLUDE_DIRS PARMETIS_LIBRARIES)
