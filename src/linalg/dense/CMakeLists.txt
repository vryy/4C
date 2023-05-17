set(OBJS
    # cmake-format: sortable
    linalg_gauss.cpp
    linalg_serialdensematrix.cpp
    linalg_serialdensevector.cpp
    linalg_utils_densematrix_communication.cpp
    linalg_utils_densematrix_determinant.cpp
    linalg_utils_densematrix_eigen.cpp
    linalg_utils_densematrix_inverse.cpp
    linalg_utils_densematrix_print.cpp
    linalg_utils_densematrix_svd.cpp
    )

baci_add_library(linalg_dense ${OBJS})

set(LIBS
    # cmake-format: sortable
    lib
    )

baci_add_dependency(linalg_dense ${LIBS})
