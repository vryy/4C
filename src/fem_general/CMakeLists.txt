set(OBJS_FEM_GENERAL
    # cmake-format: sortable
    debug_nan.cpp
    largerotations.cpp
    utils_boundary_integration.cpp
    utils_bspline.cpp
    utils_fem_shapefunctions.cpp
    utils_fem_shapefunctions.cpp
    utils_gausspoints.cpp
    utils_gausspoints.cpp
    utils_integration.cpp
    utils_integration.cpp
    utils_local_connectivity_matrices.cpp
    utils_nurbs_shapefunctions.cpp
    utils_polynomial.cpp
    utils_shapevalues_hdg.cpp
    )

baci_add_library(fem_general ${OBJS_FEM_GENERAL})

set(LIBS_FEM_GENERAL
    # cmake-format: sortable
    geometry
    lib
    lib_cut
    )

target_link_libraries(fem_general PUBLIC ${LIBS_FEM_GENERAL})
