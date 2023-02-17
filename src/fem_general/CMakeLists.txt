set(OBJS_FEM_GENERAL
    # cmake-format: sortable
    fem_general_debug_nan.cpp
    fem_general_largerotations.cpp
    fem_general_utils_boundary_integration.cpp
    fem_general_utils_bspline.cpp
    fem_general_utils_fem_shapefunctions.cpp
    fem_general_utils_fem_shapefunctions.cpp
    fem_general_utils_gauss_point_extrapolation.cpp
    fem_general_utils_gauss_point_postprocess.cpp
    fem_general_utils_gausspoints.cpp
    fem_general_utils_gausspoints.cpp
    fem_general_utils_integration.cpp
    fem_general_utils_integration.cpp
    fem_general_utils_local_connectivity_matrices.cpp
    fem_general_utils_nurbs_shapefunctions.cpp
    fem_general_utils_polynomial.cpp
    fem_general_utils_shapevalues_hdg.cpp
    )

baci_add_library(fem_general ${OBJS_FEM_GENERAL})

set(LIBS_FEM_GENERAL
    # cmake-format: sortable
    geometry
    lib
    lib_cut
    )

baci_add_dependency(fem_general ${LIBS_FEM_GENERAL})
