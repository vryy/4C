set(OBJS_FEM_GENERAL_UTILS
    # cmake-format: sortable
    discretization_fem_general_largerotations.cpp
    discretization_fem_general_utils_boundary_integration.cpp
    discretization_fem_general_utils_fem_shapefunctions.cpp
    discretization_fem_general_utils_fem_shapefunctions.cpp
    discretization_fem_general_utils_gauss_point_extrapolation.cpp
    discretization_fem_general_utils_gausspoints.cpp
    discretization_fem_general_utils_gausspoints.cpp
    discretization_fem_general_utils_integration.cpp
    discretization_fem_general_utils_integration.cpp
    discretization_fem_general_utils_local_connectivity_matrices.cpp
    discretization_fem_general_utils_polynomial.cpp
    discretization_fem_general_utils_shapevalues_hdg.cpp
    )

set(LIBS_FEM_GENERAL_UTILS
    # cmake-format: sortable
    geometry
    lib
    lib_cut
    )

baci_add_library_to_module(
  MODULE
  discretization
  NAME
  fem_general_utils
  SOURCES
  ${OBJS_FEM_GENERAL_UTILS}
  DEPENDENCIES
  ${LIBS_FEM_GENERAL_UTILS}
  )
