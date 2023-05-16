set(OBJS_LINALG
    # cmake-format: sortable
    linalg_nullspace.cpp
    linalg_utils_nullspace.cpp
    )

baci_add_library(linalg ${OBJS_LINALG})

set(LIBS_LINALG
    # cmake-format: sortable
    adapter
    lib
    )

baci_add_dependency(linalg ${LIBS_LINALG})

add_subdirectory(dense)
target_link_libraries(linalg PUBLIC linalg_dense)

add_subdirectory(sparse)
target_link_libraries(linalg PUBLIC linalg_sparse)
