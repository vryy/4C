set(OBJS_DRT_FEM_GENERAL
    # cmake-format: sortable
    debug_nan.cpp
    drt_utils_boundary_integration.cpp
    drt_utils_bspline.cpp
    drt_utils_fem_shapefunctions.cpp
    drt_utils_fem_shapefunctions.cpp
    drt_utils_gausspoints.cpp
    drt_utils_gausspoints.cpp
    drt_utils_integration.cpp
    drt_utils_integration.cpp
    drt_utils_local_connectivity_matrices.cpp
    drt_utils_nurbs_shapefunctions.cpp
    drt_utils_polynomial.cpp
    drt_utils_shapevalues_hdg.cpp
    largerotations.cpp
    )

baci_add_library(drt_fem_general ${OBJS_DRT_FEM_GENERAL})

set(LIBS_DRT_FEM_GENERAL
    # cmake-format: sortable
    drt_geometry
    drt_lib
    drt_lib_cut
    )

target_link_libraries(drt_fem_general PUBLIC ${LIBS_DRT_FEM_GENERAL})
