set(OBJS_FEM_GENERAL_NURBS
    # cmake-format: sortable
    discretization_fem_general_utils_bspline.cpp
    discretization_fem_general_utils_nurbs_shapefunctions.cpp
    )

set(LIBS_FEM_GENERAL_NURBS
    # cmake-format: sortable
    geometry
    lib
    lib_cut
    )

baci_add_library_to_module(
  MODULE
  discretization
  NAME
  fem_general_nurbs
  SOURCES
  ${OBJS_FEM_GENERAL_NURBS}
  DEPENDENCIES
  ${LIBS_FEM_GENERAL_NURBS}
  )
