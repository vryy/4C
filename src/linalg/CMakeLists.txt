set(OBJS_LINALG
    linalg_ana.cpp
    linalg_blocksparsematrix.cpp
    linalg_downwindmatrix.cpp
    linalg_equilibrate.cpp
    linalg_krylov_projector.cpp
    linalg_mapextractor.cpp
    linalg_matrixtransform.cpp
    linalg_mlapi_operator.cpp
    linalg_multiply.cpp
    linalg_nullspace.cpp
    linalg_precond.cpp
    linalg_projected_operator.cpp
    linalg_projected_precond.cpp
    linalg_serialdensematrix.cpp
    linalg_serialdensevector.cpp
    linalg_sparsematrix.cpp
    linalg_sparsematrixbase.cpp
    linalg_gauss.cpp
    linalg_utils_nullspace.cpp
    linalg_utils_densematrix_communication.cpp
    linalg_utils_densematrix_determinant.cpp
    linalg_utils_densematrix_eigen.cpp
    linalg_utils_densematrix_inverse.cpp
    linalg_utils_densematrix_print.cpp
    linalg_utils_densematrix_svd.cpp
    linalg_utils_sparse_algebra_assemble.cpp
    linalg_utils_sparse_algebra_create.cpp
    linalg_utils_sparse_algebra_manipulation.cpp
    linalg_utils_sparse_algebra_math.cpp
    linalg_utils_sparse_algebra_print.cpp
    )

baci_add_library(linalg ${OBJS_LINALG})
target_link_libraries(linalg PUBLIC drt_inpar drt_lib_utils)
