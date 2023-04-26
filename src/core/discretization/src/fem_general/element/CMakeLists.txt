set(OBJS_FEM_GENERAL_ELEMENT
    # cmake-format: sortable
    discretization_fem_general_utils_gauss_point_postprocess.cpp
    )

set(LIBS_FEM_GENERAL_ELEMENT
    # cmake-format: sortable
    geometry
    lib
    lib_cut
    )

baci_add_library_to_module(
  MODULE
  discretization
  NAME
  fem_general_element
  SOURCES
  ${OBJS_FEM_GENERAL_ELEMENT}
  DEPENDENCIES
  ${LIBS_FEM_GENERAL_ELEMENT}
  )
