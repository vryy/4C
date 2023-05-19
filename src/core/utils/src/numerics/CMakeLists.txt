set(OBJS_UTILS_NUMERICS
    # cmake-format: sortable
    utils_local_numeric_methods.cpp
    )

set(LIBS_UTILS_NUMERICS
    # cmake-format: sortable
    linalg
    )

baci_add_library_to_module(
  MODULE
  utils
  NAME
  numerics
  SOURCES
  ${OBJS_UTILS_NUMERICS}
  DEPENDENCIES
  ${LIBS_UTILS_NUMERICS}
  )
