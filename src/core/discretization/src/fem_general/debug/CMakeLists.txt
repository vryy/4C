set(OBJS_FEM_GENERAL_DEBUG
    # cmake-format: sortable
    discretization_fem_general_debug_nan.cpp
    )

set(LIBS_FEM_GENERAL_DEBUG
    # cmake-format: sortable
    geometry
    lib
    lib_cut
    )

baci_add_library_to_module(
  MODULE
  discretization
  NAME
  fem_general_debug
  SOURCES
  ${OBJS_FEM_GENERAL_DEBUG}
  DEPENDENCIES
  ${LIBS_FEM_GENERAL_DEBUG}
  )
