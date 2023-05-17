set(OBJS
    # cmake-format: sortable
    linalg_blocksparsematrix.cpp
    linalg_downwindmatrix.cpp
    linalg_equilibrate.cpp
    linalg_krylov_projector.cpp
    linalg_mapextractor.cpp
    linalg_matrixtransform.cpp
    linalg_mlapi_operator.cpp
    linalg_multiply.cpp
    linalg_projected_operator.cpp
    linalg_projected_precond.cpp
    linalg_sparsematrix.cpp
    linalg_sparsematrixbase.cpp
    linalg_utils_sparse_algebra_assemble.cpp
    linalg_utils_sparse_algebra_create.cpp
    linalg_utils_sparse_algebra_manipulation.cpp
    linalg_utils_sparse_algebra_math.cpp
    linalg_utils_sparse_algebra_print.cpp
    )

baci_add_library(linalg_sparse ${OBJS})

set(LIBS
    # cmake-format: sortable
    lib
    )

baci_add_dependency(linalg_sparse ${LIBS})
