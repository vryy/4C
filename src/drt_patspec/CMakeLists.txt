set(OBJS_DRT_PATSPEC
    # cmake-format: sortable
    patspec.cpp
    )

baci_add_library(drt_patspec ${OBJS_DRT_PATSPEC})

set(LIBS_DRT_PATSPEC
    # cmake-format: sortable
    drt_mat
    )

target_link_libraries(drt_patspec PUBLIC ${LIBS_DRT_PATSPEC})
